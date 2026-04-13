import torch
import torch.nn as nn  # 已修复：补充导入 nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, random_split
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import json
from collections import defaultdict

# 配置 matplotlib 全局字体为新罗马字体 (Times New Roman)
plt.rcParams['font.family'] = 'serif'  # 将默认字体族更改为衬线字体
plt.rcParams['font.serif'] = ['Times New Roman']  # 指定衬线字体为 Times New Roman
plt.rcParams['mathtext.fontset'] = 'stix' # 可选：让数学符号(如果有)也呈现类似 Times New Roman 的风格
plt.rcParams['axes.unicode_minus'] = False

# 导入自定义模块
from dataset import OTFSDataset
from model import Universal_OTFS_Detector

# =========================================================================
# 1. 专属损失函数定义
# =========================================================================
class ProbMapFocalLoss(nn.Module):
    """用于概率图的 Focal Loss，严厉惩罚虚警(多检测)"""
    def __init__(self, alpha=0.9, gamma=2.0):
        super().__init__()
        self.alpha = alpha 
        self.gamma = gamma

    def forward(self, pred_prob, target_prob):
        pred_prob = torch.clamp(pred_prob, min=1e-6, max=1.0-1e-6)
        
        pt_1 = torch.where(target_prob == 1.0, pred_prob, torch.ones_like(pred_prob))
        loss_pos = -self.alpha * torch.pow(1.0 - pt_1, self.gamma) * torch.log(pt_1)
        
        pt_0 = torch.where(target_prob == 0.0, 1.0 - pred_prob, torch.ones_like(pred_prob))
        loss_neg = -(1.0 - self.alpha) * torch.pow(1.0 - pt_0, self.gamma) * torch.log(pt_0)
        
        fp_penalty = 0.5 * F.relu(pred_prob - target_prob).mean()
        
        return torch.mean(loss_pos + loss_neg) + fp_penalty

class MaskedOffsetLoss(nn.Module):
    """掩膜 MSE 损失：只在有真实目标(GT)的网格计算多普勒/时延偏移误差"""
    def __init__(self):
        super().__init__()

    def forward(self, pred_offset, target_offset, target_prob_mask):
        diff = (pred_offset - target_offset) * target_prob_mask
        num_targets = torch.sum(target_prob_mask) + 1e-5 
        loss = torch.sum(diff ** 2) / num_targets
        return loss

# =========================================================================
# 2. 核心工具：全新 3 通道标签生成器
# =========================================================================
def generate_3channel_labels(labels_tensor, M=64, N=32):
    """【极速全并行版】将连续坐标转化为网络需要的图像级标签，包含冲突覆盖逻辑"""
    B, P_max, _ = labels_tensor.shape
    device = labels_tensor.device
    
    prob_gt = torch.zeros((B, 1, M, N), device=device)
    d_nu_gt = torch.zeros((B, 1, M, N), device=device)
    d_tau_gt = torch.zeros((B, 1, M, N), device=device)
    mask = torch.zeros((B, 1, M, N), device=device)
    
    valid_mask = labels_tensor[:, :, 2] > 1e-6
    if not valid_mask.any():
        return prob_gt, d_nu_gt, d_tau_gt, mask
        
    nu = labels_tensor[:, :, 0] % N
    tau = labels_tensor[:, :, 1] % M
    
    tau_int = torch.round(tau).long() % M
    nu_int = torch.round(nu).long() % N
    
    delta_tau = tau - tau_int.float()
    delta_nu = nu - nu_int.float()
    
    delta_nu = torch.where(delta_nu > N/2, delta_nu - N, delta_nu)
    delta_nu = torch.where(delta_nu < -N/2, delta_nu + N, delta_nu)
    delta_tau = torch.where(delta_tau > M/2, delta_tau - M, delta_tau)
    delta_tau = torch.where(delta_tau < -M/2, delta_tau + M, delta_tau)
    
    cond_nu_pos = delta_nu > 0.5
    delta_nu[cond_nu_pos] -= 1.0
    nu_int[cond_nu_pos] = (nu_int[cond_nu_pos] + 1) % N
    
    cond_nu_neg = delta_nu < -0.5
    delta_nu[cond_nu_neg] += 1.0
    nu_int[cond_nu_neg] = (nu_int[cond_nu_neg] - 1) % N
    
    cond_tau_pos = delta_tau > 0.5
    delta_tau[cond_tau_pos] -= 1.0
    tau_int[cond_tau_pos] = (tau_int[cond_tau_pos] + 1) % M
    
    cond_tau_neg = delta_tau < -0.5
    delta_tau[cond_tau_neg] += 1.0
    tau_int[cond_tau_neg] = (tau_int[cond_tau_neg] - 1) % M
    
    b_indices = torch.arange(B, device=device).view(B, 1).expand(B, P_max)
    
    b_valid = b_indices[valid_mask]
    nu_valid = nu_int[valid_mask]
    tau_valid = tau_int[valid_mask]
    
    delta_nu_valid = delta_nu[valid_mask]
    delta_tau_valid = delta_tau[valid_mask]
    
    prob_gt[b_valid, 0, tau_valid, nu_valid] = 1.0
    d_nu_gt[b_valid, 0, tau_valid, nu_valid] = delta_nu_valid
    d_tau_gt[b_valid, 0, tau_valid, nu_valid] = delta_tau_valid
    mask[b_valid, 0, tau_valid, nu_valid] = 1.0
            
    return prob_gt, d_nu_gt, d_tau_gt, mask

# =========================================================================
# 3. 可视化工具
# =========================================================================
def visualize_dynamic_scatter(est_params, gt_labels, epoch, N=32, M=64):
    """ 
    est_params: [B, K, 4] -> (nu, tau, amp, phase) 
    带高级数值误差分析表格的可视化 (支持 GT/Est 自动配对)
    """
    if not os.path.exists('debug_pos'): os.makedirs('debug_pos')
    
    # 只可视化 batch 中的第一个样本
    est = est_params[0].detach().cpu().numpy()
    gt = gt_labels[0].detach().cpu().numpy()
    
    # 过滤无效点
    gt = gt[gt[:, 2] > 1e-6] 
    if len(est) > 0:
        est = est[est[:, 2] > 1e-4] 
    
    if len(gt) > 0:  gt = gt[np.argsort(gt[:, 1])]
    
    # ==========================================
    # 1. 画布与排版布局设定 (1行2列)
    # ==========================================
    fig = plt.figure(figsize=(16, 7))
    gs = GridSpec(1, 2, width_ratios=[1, 1.2], wspace=0.1)
    
    ax_plot = fig.add_subplot(gs[0])
    ax_table = fig.add_subplot(gs[1])
    ax_table.axis('off') # 隐藏表格区域的坐标轴
    
    # ==========================================
    # 2. 绘制左侧：检测星座图
    # ==========================================
    gt_nu = gt[:, 0] % N if len(gt) > 0 else []
    gt_tau = gt[:, 1] % M if len(gt) > 0 else []
    
    if len(gt) > 0:
        ax_plot.scatter(gt_nu, gt_tau, c='#77AC70', s=300, alpha=0.8, edgecolors='black', linewidth=1.5, label=f'Ground Truth (K={len(gt)})', zorder=3)
    
    if len(est) > 0:
        est_nu = est[:, 0] % N
        est_tau = est[:, 1] % M
        ax_plot.scatter(est_nu, est_tau, c='r', s=200, marker='x', linewidth=3, label=f'Estimated (K={len(est)})', zorder=4)
    
    ax_plot.set_xlim(-1, N)
    ax_plot.set_ylim(-1, 10)
    ax_plot.set_xlabel('Doppler Index')
    ax_plot.set_ylabel('Delay Index')
    ax_plot.set_title(f'Detection Constellation (Epoch {epoch})', pad=15, fontweight='bold')
    ax_plot.legend(loc='lower right')
    ax_plot.grid(True, linestyle='--', alpha=0.4)
    ax_plot.invert_yaxis()
    
    # ==========================================
    # 3. 自动配对算法 (GT 与 Est 距离最小匹配)
    # ==========================================
    matched_est = []
    unmatched_est = list(range(len(est)))
    pairs = [] # 格式: [(gt_idx, est_idx), ...]

    for i, g in enumerate(gt):
        if not unmatched_est:
            pairs.append((i, None)) # 漏检
            continue
            
        # 计算环形折叠距离
        dists = []
        for j in unmatched_est:
            d_nu = min(abs(g[0] - est[j][0]), N - abs(g[0] - est[j][0]))
            d_tau = min(abs(g[1] - est[j][1]), M - abs(g[1] - est[j][1]))
            dists.append(np.sqrt(d_nu**2 + d_tau**2))
            
        min_idx = np.argmin(dists)
        if dists[min_idx] < 3.0: # 匹配阈值
            pairs.append((i, unmatched_est[min_idx]))
            unmatched_est.pop(min_idx)
        else:
            pairs.append((i, None)) # 距离太远，视为漏检
            
    # 剩下的都是虚警 (False Alarm)
    for j in unmatched_est:
        pairs.append((None, j))

    # ==========================================
    # 4. 构建表格数据与色彩矩阵
    # ==========================================
    col_labels = ['Path', 'Type', 'Delay (Err)', 'Doppler (Err)', 'Complex Gain']
    cell_text = []
    cell_colors = []
    
    # 配色方案 (仿造截图高颜值配色)
    c_head = '#F4D03F' # 表头黄
    c_gt = '#E2EFDA'   # GT 浅绿
    c_est = '#FCE4D6'  # Est 浅橙
    c_path = '#FFF2CC' # 路径编号浅黄

    path_idx = 1
    for gt_idx, est_idx in pairs:
        # ---- 添加 GT 行 ----
        if gt_idx is not None:
            g = gt[gt_idx]
            g_comp = g[2] * np.exp(1j * g[3])
            cell_text.append([f"T{path_idx}", "GT", f"{g[1]:.4f}", f"{g[0]:.4f}", f"{g_comp.real:.4f}{g_comp.imag:+.4f}j"])
            cell_colors.append([c_path, c_gt, c_gt, c_gt, c_gt])
        else:
            # 虚警对应的空 GT
            cell_text.append([f"FA{path_idx}", "GT", "-", "-", "-"])
            cell_colors.append([c_path, c_gt, c_gt, c_gt, c_gt])

        # ---- 添加 Est 行 ----
        if est_idx is not None:
            e = est[est_idx]
            e_comp = e[2] * np.exp(1j * e[3])
            
            if gt_idx is not None:
                # 正常配对，计算误差
                err_tau = min(abs(e[1] - g[1]), M - abs(e[1] - g[1]))
                err_nu = min(abs(e[0] - g[0]), N - abs(e[0] - g[0]))
                cell_text.append(["", "Est", f"{e[1]:.4f} ({err_tau:.4f})", f"{e[0]:.4f} ({err_nu:.4f})", f"{e_comp.real:.4f}{e_comp.imag:+.4f}j"])
            else:
                # 虚警，无误差对比
                cell_text.append(["", "Est (FA)", f"{e[1]:.4f}", f"{e[0]:.4f}", f"{e_comp.real:.4f}{e_comp.imag:+.4f}j"])
                
            cell_colors.append([c_path, c_est, c_est, c_est, c_est])
        else:
            # 漏检
            cell_text.append(["", "Est", "Missed", "Missed", "Missed"])
            cell_colors.append([c_path, c_est, c_est, c_est, c_est])

        path_idx += 1

    # ==========================================
    # 5. 绘制精美表格
    # ==========================================
    if len(cell_text) > 0:
            # 1. 改回 loc='center'，让表格从正中心向上下展开，最稳定
            table = ax_table.table(cellText=cell_text, 
                                colLabels=col_labels, 
                                cellColours=cell_colors,
                                colColours=[c_head]*5,
                                loc='center', 
                                cellLoc='center')
            
            table.auto_set_font_size(False)
            table.set_fontsize(11)
            table.scale(1, 2.5) # 拉伸表格行高，更美观
            
            # 将表头文字加粗
            for (row, col), cell in table.get_celld().items():
                if row == 0:
                    cell.set_text_props(weight='bold')
                    
            # 2. 【核心修复】动态计算标题高度！
            # 每多一行，就把标题往上推 0.05 的相对高度
            num_rows = len(cell_text) + 1  # 数据行数 + 1行表头
            title_y = 0.5 + (num_rows * 0.035) + 0.05
            
            ax_table.set_title('Numerical Error Analysis', fontweight='bold', fontsize=14, y=title_y)
                
    else:
        ax_table.set_title('Numerical Error Analysis', pad=15, fontweight='bold', fontsize=14)
        ax_table.text(0.5, 0.5, "No Targets Detected", ha='center', va='center', fontsize=14)
    
    # 保存图片
    plt.tight_layout()
    save_path = f'debug_pos/analysis_epoch_{epoch}.png'
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

def save_loss_curve(train_losses, val_losses, filename='loss_curve.png', start_epoch=5):
    total_epochs = len(train_losses)
    if total_epochs < start_epoch:
        start_epoch = 1
        
    start_idx = start_epoch - 1 
    
    plot_train = train_losses[start_idx:]
    plot_val = val_losses[start_idx:]
    epochs = range(start_epoch, start_epoch + len(plot_train))
    
    plt.figure(figsize=(10, 5))
    plt.plot(epochs, plot_train, label='Train Total Loss', linewidth=2)
    plt.plot(epochs, plot_val, label='Val Total Loss', linewidth=2, linestyle='--')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title(f'Loss Curve (Starting from Epoch {start_epoch})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=100, bbox_inches='tight')
    plt.close()

# =========================================================================
# 4. 训练主流程 (纯净无噪版 + Teacher Forcing)
# =========================================================================
if __name__ == "__main__":
    BATCH_SIZE = 32
    LR = 1e-4  
    EPOCHS = 100  
    DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    NUM_WORKERS = min(4, os.cpu_count() or 1)
    PIN_MEMORY = DEVICE.type == 'cuda'
    PERSISTENT_WORKERS = NUM_WORKERS > 0
    USE_SEPARATE_VAL = False
    
    train_h5 = 'data_h5/OTFS_Train_Set_Physical.h5'
    val_h5 = ''
    
    train_dataset = OTFSDataset(train_h5, normalize=True)
    if USE_SEPARATE_VAL and os.path.exists(val_h5):
        val_dataset = OTFSDataset(val_h5, normalize=True)
    else:
        total = len(train_dataset)
        train_size = int(0.8 * total)
        val_size = total - train_size
        train_dataset, val_dataset = random_split(train_dataset, [train_size, val_size])

    train_loader = DataLoader(
        train_dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=NUM_WORKERS,
        pin_memory=PIN_MEMORY,
        persistent_workers=PERSISTENT_WORKERS
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=NUM_WORKERS,
        pin_memory=PIN_MEMORY,
        persistent_workers=PERSISTENT_WORKERS
    )
    
    model = Universal_OTFS_Detector(M=64, N=32).to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=LR, weight_decay=1e-4)  
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS, eta_min=1e-7)
    
    criterion_prob = ProbMapFocalLoss(alpha=0.9, gamma=2.0).to(DEVICE)
    criterion_offset = MaskedOffsetLoss().to(DEVICE)
    
    loss_history = {'train': [], 'val': []}
    
    best_val_loss = float('inf')
    patience_counter = 0
    PATIENCE = 20  
    
    print(f"Start Pure Physics-Driven Unfolded Training on {DEVICE}...")
    print(f"Early Stopping: Patience={PATIENCE} epochs\n")
    
    for epoch in range(EPOCHS):
        # --- Training ---
        model.train()
        train_loss_accum = 0.0
        

        for batch_idx, (y_in, x_in, labels, sy, sx, _) in enumerate(train_loader):
            y_in = y_in.to(DEVICE, non_blocking=PIN_MEMORY)
            labels = labels.to(DEVICE, non_blocking=PIN_MEMORY)
            B = y_in.shape[0]

            optimizer.zero_grad()
            
            # ==========================================================
            # 【极致清爽训练：纯净原图直出】
            # 因为贪心 SIC 完美处理了所有的迭代减法，
            # 网络现在只需要专心训练一个技能：“从当前图中找到最显著的波峰”
            # ==========================================================
            
            # 永远喂给网络干净的历史状态（假装这是第一轮寻优）
            y_recon_input = torch.zeros_like(y_in)
            prob_prev_input = torch.zeros((B, 1, 64, 32), device=DEVICE) 
            
            # 直接生成完整的上帝视角答案（网络只需尽力把图上的目标都亮起来）
            prob_gt, d_nu_gt, d_tau_gt, mask = generate_3channel_labels(labels, M=64, N=32)

            # --- 前向传播 ---
            prob_pred, d_nu_pred, d_tau_pred = model(y_in, y_recon_input, prob_prev_input)
            
            # --- 损失计算 ---
            loss_prob = criterion_prob(prob_pred, prob_gt)
            loss_nu = criterion_offset(d_nu_pred, d_nu_gt, mask)
            loss_tau = criterion_offset(d_tau_pred, d_tau_gt, mask)
            
            # 把亚像素回归的权重放大 10 倍，逼迫网络练好“火眼金睛”
            loss = loss_prob + 10.0 * (loss_nu + loss_tau)
            
            # --- 反向传播 ---
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()       
            
            train_loss_accum += loss.item()
            
            if batch_idx % 50 == 0:
                print(f"Ep {epoch+1} B {batch_idx} | Total: {loss.item():.4f} [Prob:{loss_prob.item():.4f}, dNu:{loss_nu.item():.4f}, dTau:{loss_tau.item():.4f}]")

        avg_train_loss = train_loss_accum / len(train_loader)
        loss_history['train'].append(avg_train_loss)
        
        # --- Validation ---
        model.eval()
        val_loss_accum = 0.0
        last_val_data = None 
        
        with torch.no_grad():
            for batch_idx, (y_in, x_in, labels, _, _, _) in enumerate(val_loader):
                y_in = y_in.to(DEVICE, non_blocking=PIN_MEMORY)
                labels = labels.to(DEVICE, non_blocking=PIN_MEMORY)
                
                # --- 验证集的损失计算保持不变，依然用 Teacher-Forcing 模式 ---
                # 这确保了我们衡量的是网络在理想单步预测下的性能，与训练过程一致
                y_recon_dummy = torch.zeros_like(y_in)
                prob_prev_dummy = torch.zeros((y_in.shape[0], 1, 64, 32), device=DEVICE)
                
                prob_pred, d_nu_pred, d_tau_pred = model(y_in, y_recon_dummy, prob_prev_dummy)
                prob_gt, d_nu_gt, d_tau_gt, mask = generate_3channel_labels(labels, M=64, N=32)
                
                loss_prob = criterion_prob(prob_pred, prob_gt)
                loss_nu = criterion_offset(d_nu_pred, d_nu_gt, mask)
                loss_tau = criterion_offset(d_tau_pred, d_tau_gt, mask)
                
                loss = loss_prob + 10.0 * (loss_nu + loss_tau)
                val_loss_accum += loss.item()
                
                if batch_idx == 0:
                    val_prob_loss = loss_prob.item()
                    val_nu_loss = loss_nu.item()
                    val_tau_loss = loss_tau.item()

                # --- 【重要】只在验证的最后一个 batch 且特定 epoch 运行新的、耗时较长的 inference ---
                if batch_idx == len(val_loader) - 1 and (epoch + 1) % 10 == 0:
                    print("  [Info] Running SIC-Adam inference for visualization...")
                    # 为了节省时间，只对验证集中的前4个样本进行完整推理
                    x_in = x_in.to(DEVICE, non_blocking=PIN_MEMORY)
                    est_params_inf, _, _ = model.inference(y_in[:4], x_in[:4], max_targets=8, prob_threshold=0.05, gain_ratio_stop=10.0)
                    last_val_data = (est_params_inf, labels[:4])#可视化的数据准备
                
        avg_val_loss = val_loss_accum / len(val_loader)
        loss_history['val'].append(avg_val_loss)
        
        print(f"  Val | Total: {avg_val_loss:.4f} [Prob:{val_prob_loss:.4f}, dNu:{val_nu_loss:.4f}, dTau:{val_tau_loss:.4f}]")
        
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            patience_counter = 0
            torch.save(model.state_dict(), 'universal_otfs_stripper_best.pth')
            print(f"  --> New best model saved! Val Loss: {best_val_loss:.5f}")
        else:
            patience_counter += 1
            if patience_counter >= PATIENCE:
                print(f"\n⚠ Early Stopping triggered! No improvement for {PATIENCE} epochs.")
                print(f"Best Val Loss: {best_val_loss:.5f}")
                break
        
        scheduler.step()  
        current_lr = optimizer.param_groups[0]['lr']
        print(f"Epoch {epoch+1}/{EPOCHS} | Train Loss: {avg_train_loss:.5f} | Val Loss: {avg_val_loss:.5f} | LR: {current_lr:.2e} | Patience: {patience_counter}/{PATIENCE}")
        
        if (epoch + 1) % 10 == 0 and last_val_data:
            visualize_dynamic_scatter(last_val_data[0], last_val_data[1], epoch+1)
                
    torch.save(model.state_dict(), 'universal_otfs_stripper.pth')
    loss_data = np.array([loss_history['train'], loss_history['val']])
    if not os.path.exists('data_npy'): os.makedirs('data_npy')
    np.save('data_npy/loss_history.npy', loss_data)
    save_loss_curve(loss_history['train'], loss_history['val'])
    
    print("Training finished.")
