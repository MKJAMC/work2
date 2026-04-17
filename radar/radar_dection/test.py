import torch
import torch.backends.cudnn as cudnn
cudnn.benchmark = True
import numpy as np
import os
import time
from torch.utils.data import DataLoader
import warnings  # <--- 新增
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib") # <--- 新增，让 matplotlib 闭嘴

# =========================================================================
# 【核心优化】：开启 Matplotlib 无头模式，极速后台画图，告别假死阻塞
# =========================================================================
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec  # 新增：用于表格与图像的并排布局

# 配置 matplotlib 支持中文显示及全局字体
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'Times New Roman'] 
plt.rcParams['mathtext.fontset'] = 'stix'  # 保证图中的公式和变量（如果有）依然是优美的类似 Times 的字体
plt.rcParams['axes.unicode_minus'] = False # 解决负号 '-' 显示为方块的问题

from dataset import OTFSDataset
from model import Universal_OTFS_Detector

# =========================================================================
# 配置参数
# =========================================================================
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
N = 32          # Doppler 维度
M = 64          # Delay 维度
BATCH_SIZE = 128
FAST_MODE = True
model_path = 'universal_otfs_stripper_best.pth'  # 建议测试 best 模型
test_h5 = 'data_h5/OTFS_Test_Set_T1T2_1000.h5'   
# test_h5 = 'data_h5/OTFS_Test_Set_Scatter_1000.h5'

# ================= 物理参数与分辨率计算 =================
fc = 64e9        # 载波频率 (Hz)
delta_f = 120e3  # 子载波间隔 (Hz)
c = 3e8          # 光速 (m/s)
distance_res = c / (2 * M * delta_f)   # 距离分辨率: m/bin
v_res = c * delta_f / (2 * fc * N)     # 速度分辨率: m/s/bin

# =========================================================================
# 5.5 按照精确SNR寻找最佳样本 (F1最高)
# =========================================================================
def update_best_by_exact_snr(best_dict, snr_db, sample_idx, metrics, est_params, gt_params):
    snr_key = int(round(snr_db))
    current = best_dict.get(snr_key)
    f1 = metrics.get('F1', 0.0)
    
    # 如果当前SNR还没记录过，或者遇到了F1更高的样本，则更新记录
    if current is None or f1 > current['f1']:
        best_dict[snr_key] = {
            'sample_idx': sample_idx, 'snr_db': snr_db, 'f1': f1,
            'precision': metrics.get('Precision', 0.0),
            'recall': metrics.get('Recall', 0.0),
            'tp': metrics.get('TP', 0), 'fp': metrics.get('FP', 0), 'fn': metrics.get('FN', 0),
            'est_params': est_params, 'gt_params': gt_params
        }

def collect_samples_by_exact_snr(samples_by_snr, snr_db, sample_idx, metrics, est_params, gt_params):
    snr_key = int(round(snr_db))
    if snr_key not in samples_by_snr:
        samples_by_snr[snr_key] = []

    samples_by_snr[snr_key].append({
        'sample_idx': sample_idx, 'snr_db': snr_db, 'f1': metrics.get('F1', 0.0),
        'precision': metrics.get('Precision', 0.0),
        'recall': metrics.get('Recall', 0.0),
        'tp': metrics.get('TP', 0), 'fp': metrics.get('FP', 0), 'fn': metrics.get('FN', 0),
        'est_params': est_params, 'gt_params': gt_params
    })

def get_middle_by_exact_snr(samples_by_snr):
    middle_by_exact_snr = {}
    for snr_key, sample_list in samples_by_snr.items():
        if len(sample_list) == 0:
            continue

        # 选择该SNR下F1中位数对应的样本作为“中间样本”
        sorted_samples = sorted(sample_list, key=lambda x: x['f1'])
        mid_idx = len(sorted_samples) // 2
        middle_by_exact_snr[snr_key] = sorted_samples[mid_idx]

    return middle_by_exact_snr



def _build_circular_dist_matrix(est_pos_np, gt_pos_np):
    # 向量化计算循环维度上的欧氏距离矩阵
    d_nu = np.abs(est_pos_np[:, None, 0] - gt_pos_np[None, :, 0])
    d_nu = np.minimum(d_nu, N - d_nu)
    d_tau = np.abs(est_pos_np[:, None, 1] - gt_pos_np[None, :, 1])
    d_tau = np.minimum(d_tau, M - d_tau)
    return np.sqrt(d_nu**2 + d_tau**2)

def _greedy_match_indices(dist_matrix, threshold):
    K_est, K_gt = dist_matrix.shape
    if K_est == 0 or K_gt == 0:
        return []

    working = dist_matrix.copy()
    matched_pairs = []

    for _ in range(min(K_est, K_gt)):
        flat_idx = np.argmin(working)
        min_dist = float(working.flat[flat_idx])
        if not np.isfinite(min_dist) or min_dist > threshold:
            break

        i, j = np.unravel_index(flat_idx, working.shape)
        matched_pairs.append((i, j))
        working[i, :] = np.inf
        working[:, j] = np.inf

    return matched_pairs

def compute_detection_and_rmse(est_params, gt_labels, det_threshold=1.0, rmse_threshold=1.5):
    est_pos = est_params[:, :2]
    gt_pos = gt_labels[:, :2]

    K_est = est_pos.shape[0]
    K_gt = gt_pos.shape[0]

    if K_est == 0 and K_gt == 0:
        metrics = {'TP': 0, 'FP': 0, 'FN': 0, 'Precision': 1.0, 'Recall': 1.0, 'F1': 1.0}
        return metrics, np.nan, np.nan, 0, np.array([]), np.array([])
    if K_est == 0:
        metrics = {'TP': 0, 'FP': 0, 'FN': K_gt, 'Precision': 0.0, 'Recall': 0.0, 'F1': 0.0}
        return metrics, np.nan, np.nan, 0, np.array([]), np.array([])
    if K_gt == 0:
        metrics = {'TP': 0, 'FP': K_est, 'FN': 0, 'Precision': 0.0, 'Recall': 0.0, 'F1': 0.0}
        return metrics, np.nan, np.nan, 0, np.array([]), np.array([])

    est_pos_np = est_pos.detach().cpu().numpy() if torch.is_tensor(est_pos) else est_pos
    gt_pos_np = gt_pos.detach().cpu().numpy() if torch.is_tensor(gt_pos) else gt_pos
    dist_matrix = _build_circular_dist_matrix(est_pos_np, gt_pos_np)

    det_pairs = _greedy_match_indices(dist_matrix, det_threshold)
    TP = len(det_pairs)
    FP = K_est - TP
    FN = K_gt - TP
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    metrics = {'TP': TP, 'FP': FP, 'FN': FN, 'Precision': precision, 'Recall': recall, 'F1': f1}

    rmse_pairs = _greedy_match_indices(dist_matrix, rmse_threshold)
    if len(rmse_pairs) == 0:
        return metrics, np.nan, np.nan, 0, np.array([]), np.array([])

    est_idx = np.array([p[0] for p in rmse_pairs], dtype=int)
    gt_idx = np.array([p[1] for p in rmse_pairs], dtype=int)

    nu_errors = est_pos_np[est_idx, 0] - gt_pos_np[gt_idx, 0]
    nu_errors = (nu_errors + N / 2) % N - N / 2
    tau_errors = est_pos_np[est_idx, 1] - gt_pos_np[gt_idx, 1]
    tau_errors = (tau_errors + M / 2) % M - M / 2

    rmse_tau = float(np.sqrt(np.mean(tau_errors**2)))
    rmse_nu = float(np.sqrt(np.mean(nu_errors**2)))
    return metrics, rmse_tau, rmse_nu, len(rmse_pairs), tau_errors, nu_errors

# =========================================================================
# 3. 绘制检测结果 (带数值误差分析表格) -> 融合了 Train 的代码
# =========================================================================
def plot_detection_result(est_params, gt_labels, sample_idx, title_suffix="", save_dir='test_results_scatter'):
    os.makedirs(save_dir, exist_ok=True)
    
    est = est_params.cpu().numpy() if torch.is_tensor(est_params) else est_params.copy()
    gt = gt_labels.cpu().numpy() if torch.is_tensor(gt_labels) else gt_labels.copy()
    
    # 过滤无效点
    if len(gt) > 0 and gt.shape[1] >= 3:
        gt = gt[gt[:, 2] > 1e-6] 
    if len(est) > 0 and est.shape[1] >= 3:
        est = est[est[:, 2] > 1e-5] 
    
    if len(gt) > 0:  
        gt = gt[np.argsort(gt[:, 1])] # 按 Delay 排序，方便标号
    
    # ==========================================
    # 1. 画布与排版布局设定 (1行2列)
    # ==========================================
    fig = plt.figure(figsize=(16, 7))
    gs = GridSpec(1, 2, width_ratios=[1, 1.2], wspace=0.1)
    
    ax_plot = fig.add_subplot(gs[0])
    ax_table = fig.add_subplot(gs[1])
    ax_table.axis('off') 
    
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
    ax_plot.set_xlabel('Doppler Index', fontsize=12, fontweight='bold')
    ax_plot.set_ylabel('Delay Index', fontsize=12, fontweight='bold')
    title = f'Detection Result - Sample {sample_idx}'
    if title_suffix:
        title += f'\n({title_suffix})'
    ax_plot.set_title(title, pad=15, fontsize=14, fontweight='bold')
    ax_plot.legend(loc='lower right')
    ax_plot.grid(True, linestyle='--', alpha=0.4)
    ax_plot.invert_yaxis()
    
    # ==========================================
    # 3. 自动配对算法 (GT 与 Est 距离最小匹配)
    # ==========================================
    matched_est = []
    unmatched_est = list(range(len(est)))
    pairs = [] 

    for i, g in enumerate(gt):
        if not unmatched_est:
            pairs.append((i, None)) # 漏检
            continue
            
        dists = []
        for j in unmatched_est:
            d_nu = min(abs(g[0] - est[j][0]), N - abs(g[0] - est[j][0]))
            d_tau = min(abs(g[1] - est[j][1]), M - abs(g[1] - est[j][1]))
            dists.append(np.sqrt(d_nu**2 + d_tau**2))
            
        min_idx = np.argmin(dists)
        if dists[min_idx] < 3.0: 
            pairs.append((i, unmatched_est[min_idx]))
            unmatched_est.pop(min_idx)
        else:
            pairs.append((i, None)) 
            
    for j in unmatched_est:
        pairs.append((None, j))

    # ==========================================
    # 4. 构建表格数据与色彩矩阵
    # ==========================================
    col_labels = ['Path', 'Type', 'Delay (Err)', 'Doppler (Err)', 'Complex Gain']
    cell_text = []
    cell_colors = []
    
    c_head = '#F4D03F' 
    c_gt = '#E2EFDA'   
    c_est = '#FCE4D6'  
    c_path = '#FFF2CC' 

    path_idx = 1
    for gt_idx, est_idx in pairs:
        # ---- 添加 GT 行 ----
        if gt_idx is not None:
            g = gt[gt_idx]
            g_comp = g[2] * np.exp(1j * g[3]) if len(g) >= 4 else 0j
            cell_text.append([f"T{path_idx}", "GT", f"{g[1]:.4f}", f"{g[0]:.4f}", f"{g_comp.real:.4f}{g_comp.imag:+.4f}j"])
            cell_colors.append([c_path, c_gt, c_gt, c_gt, c_gt])
        else:
            cell_text.append([f"FA{path_idx}", "GT", "-", "-", "-"])
            cell_colors.append([c_path, c_gt, c_gt, c_gt, c_gt])

        # ---- 添加 Est 行 ----
        if est_idx is not None:
            e = est[est_idx]
            e_comp = e[2] * np.exp(1j * e[3]) if len(e) >= 4 else 0j
            
            if gt_idx is not None:
                err_tau = min(abs(e[1] - g[1]), M - abs(e[1] - g[1]))
                err_nu = min(abs(e[0] - g[0]), N - abs(e[0] - g[0]))
                cell_text.append(["", "Est", f"{e[1]:.4f} ({err_tau:.4f})", f"{e[0]:.4f} ({err_nu:.4f})", f"{e_comp.real:.4f}{e_comp.imag:+.4f}j"])
            else:
                cell_text.append(["", "Est (FA)", f"{e[1]:.4f}", f"{e[0]:.4f}", f"{e_comp.real:.4f}{e_comp.imag:+.4f}j"])
                
            cell_colors.append([c_path, c_est, c_est, c_est, c_est])
        else:
            cell_text.append(["", "Est", "Missed", "Missed", "Missed"])
            cell_colors.append([c_path, c_est, c_est, c_est, c_est])

        path_idx += 1

    # ==========================================
    # 5. 绘制精美表格
    # ==========================================
    if len(cell_text) > 0:
        table = ax_table.table(cellText=cell_text, 
                               colLabels=col_labels, 
                               cellColours=cell_colors,
                               colColours=[c_head]*5,
                               loc='center', 
                               cellLoc='center')
        
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 2.5) 
        
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_text_props(weight='bold')
                
        num_rows = len(cell_text) + 1  
        title_y = 0.5 + (num_rows * 0.035) + 0.05
        ax_table.set_title('Numerical Error Analysis', fontweight='bold', fontsize=14, y=title_y)
            
    else:
        ax_table.set_title('Numerical Error Analysis', pad=15, fontweight='bold', fontsize=14)
        ax_table.text(0.5, 0.5, "No Targets Detected", ha='center', va='center', fontsize=14)
    
    save_path = os.path.join(save_dir, f'detection_sample_{sample_idx}.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    return save_path

# =========================================================================
# 4. 绘制RMSE vs SNR曲线
# =========================================================================
def plot_rmse_vs_snr(snr_rmse_data, save_dir='test_results'):
    os.makedirs(save_dir, exist_ok=True)
    if len(snr_rmse_data) == 0:
        return None
    
    snr_values = np.array([d['snr'] for d in snr_rmse_data])
    rmse_tau_values = np.array([d['rmse_tau'] for d in snr_rmse_data])
    rmse_nu_values = np.array([d['rmse_nu'] for d in snr_rmse_data])
    
    sort_idx = np.argsort(snr_values)
    snr_sorted = snr_values[sort_idx]
    rmse_tau_m = rmse_tau_values[sort_idx] * distance_res
    rmse_nu_mps = rmse_nu_values[sort_idx] * v_res
    
    unique_snrs = np.unique(np.round(snr_sorted))
    rmse_tau_binned, rmse_nu_binned = [], []
    
    for snr in unique_snrs:
        mask = np.abs(snr_sorted - snr) < 0.5
        rmse_tau_binned.append(np.mean(rmse_tau_m[mask]))
        rmse_nu_binned.append(np.mean(rmse_nu_mps[mask]))
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    ax1.scatter(snr_sorted, rmse_tau_m, alpha=0.5, s=30, label='样本数据', color='lightblue', edgecolors='blue')
    ax1.plot(unique_snrs, rmse_tau_binned, 'o-', linewidth=2.5, markersize=8, color='red', label='分组平均')
    ax1.set_xlabel('信噪比 (dB)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('均方根误差 (m)', fontsize=12, fontweight='bold')
    ax1.set_title('延迟RMSE随SNR变化', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.3)
    
    ax2.scatter(snr_sorted, rmse_nu_mps, alpha=0.5, s=30, label='样本数据', color='lightgreen', edgecolors='green')
    ax2.plot(unique_snrs, rmse_nu_binned, 'o-', linewidth=2.5, markersize=8, color='red', label='分组平均')
    ax2.set_xlabel('信噪比 (dB)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('均方根误差 (m/s)', fontsize=12, fontweight='bold')
    ax2.set_title('多普勒RMSE随SNR变化', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.3)
    
    save_path = os.path.join(save_dir, 'rmse_vs_snr.png')
    plt.savefig(save_path, dpi=100, bbox_inches='tight')
    plt.close()
    return save_path

# =========================================================================
# 4.1 终端打印RMSE vs SNR数值
# =========================================================================
def print_rmse_vs_snr(snr_rmse_data):
    if len(snr_rmse_data) == 0:
        print("\n[RMSE] 无可用匹配目标，无法统计RMSE随SNR的数值。")
        return

    snr_values = np.array([d['snr'] for d in snr_rmse_data])
    rmse_tau_values = np.array([d['rmse_tau'] for d in snr_rmse_data])
    rmse_nu_values = np.array([d['rmse_nu'] for d in snr_rmse_data])

    sort_idx = np.argsort(snr_values)
    snr_sorted = snr_values[sort_idx]
    rmse_tau_m = rmse_tau_values[sort_idx] * distance_res
    rmse_nu_mps = rmse_nu_values[sort_idx] * v_res

    unique_snrs = np.unique(np.round(snr_sorted).astype(int))

    print("\n[RMSE] 按SNR分组的平均RMSE数值：")
    print("  SNR(dB) | Delay RMSE(m) | Doppler RMSE(m/s) | Samples")
    print("  " + "-" * 55)

    for snr in unique_snrs:
        mask = np.abs(snr_sorted - snr) < 0.5
        count = int(np.sum(mask))
        if count == 0:
            continue

        tau_mean = float(np.mean(rmse_tau_m[mask]))
        nu_mean = float(np.mean(rmse_nu_mps[mask]))
        print(f"  {snr:7d} | {tau_mean:13.6f} | {nu_mean:17.6f} | {count:7d}")

# =========================================================================
# 5. 按照精确SNR寻找最差样本
# =========================================================================
def update_worst_by_exact_snr(worst_dict, snr_db, sample_idx, metrics, est_params, gt_params):
    # 修改：原先传的是 est_pos，现在统一传包含amp,phase的完整 params
    snr_key = int(round(snr_db))
    current = worst_dict.get(snr_key)
    f1 = metrics.get('F1', 0.0)
    
    if current is None or f1 < current['f1']:
        worst_dict[snr_key] = {
            'sample_idx': sample_idx, 'snr_db': snr_db, 'f1': f1,
            'precision': metrics.get('Precision', 0.0),
            'recall': metrics.get('Recall', 0.0),
            'tp': metrics.get('TP', 0), 'fp': metrics.get('FP', 0), 'fn': metrics.get('FN', 0),
            'est_params': est_params, 'gt_params': gt_params
        }

# =========================================================================
# 主程序
# =========================================================================
if __name__ == "__main__":
    print(f"[PARAM] Distance resolution: {distance_res:.6f} m/bin")
    print(f"[PARAM] Velocity resolution: {v_res:.6f} m/s/bin\n")
    print(f"Testing on {DEVICE}...")
    num_workers = min(4, os.cpu_count() or 1)
    pin_memory = DEVICE.type == 'cuda'
    persistent_workers = num_workers > 0
    prefetch_factor = 4 if num_workers > 0 else None
    
    model = Universal_OTFS_Detector(M=M, N=N).to(DEVICE)
    model.load_state_dict(torch.load(model_path, map_location=DEVICE), strict=False)
    model.eval()
    print(f"✓ 模型加载成功\n")
    
    test_dataset = OTFSDataset(test_h5, normalize=True)
    test_loader = DataLoader(
        test_dataset,
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=pin_memory,
        persistent_workers=persistent_workers,
        prefetch_factor=prefetch_factor
    )
     
    all_metrics = []
    snr_rmse_data = [] 
    worst_by_exact_snr = {}
    best_by_exact_snr = {}
    samples_by_exact_snr = {}

    # 仅统计 F1 != 1.0 的样本中，估计目标数与GT目标数的关系
    non_perfect_est_more_count = 0
    non_perfect_est_less_count = 0
    non_perfect_est_equal_count = 0

    # 记录 Delay / Doppler 误差最差样本（基于 RMSE 匹配结果）
    worst_delay_sample = None
    worst_doppler_sample = None
    
    sample_count = 0
    t_start = time.time()
    
    print("=" * 70)
    print("开始高速测试推理 (日志按Batch打印)...")
    print("=" * 70 + "\n")
    
    with torch.no_grad():
        for batch_idx, (y_in, x_in, labels, _, _, snr_values) in enumerate(test_loader):
            y_in = y_in.to(DEVICE, non_blocking=pin_memory)
            x_in = x_in.to(DEVICE, non_blocking=pin_memory)
            labels = labels.to(DEVICE, non_blocking=pin_memory)
            snr_values = snr_values.to(DEVICE, non_blocking=pin_memory)
            
            B = y_in.shape[0]
            
            est_params, _, _ = model.inference(
                y_in, x_in,
                max_targets=6,
                prob_threshold=0.04,
                close_bin_threshold=1.8,
                stop_gain_ratio=7.0,
            )
                        

            for b in range(B):
                snr_db = snr_values[b].item()
                
                # 提取完整的 4 维估计参数 (nu, tau, amp, phase)
                est_params_b = est_params[b] 
                valid_est_mask = est_params_b[:, 2] > 1e-5 
                est_params_b = est_params_b[valid_est_mask]
                est_pos_b = est_params_b[:, :2]
                
                # 提取完整的 4 维真实参数
                gt_params_b_all = labels[b]
                valid_gt_mask = gt_params_b_all[:, 2] > 1e-8 
                gt_pos_b = gt_params_b_all[valid_gt_mask, :2] 
                gt_labels_b = gt_params_b_all[valid_gt_mask, :]
                
                metrics, rmse_tau, rmse_nu, n_matched, _, _ = compute_detection_and_rmse(
                    est_params_b, gt_labels_b, det_threshold=1.0, rmse_threshold=1.5
                )
                all_metrics.append(metrics)

                # 统计：仅在 F1 != 1.0 时，比较估计样本数与GT样本数
                if metrics['F1'] != 1.0:
                    est_count_b = est_params_b.shape[0]
                    gt_count_b = gt_labels_b.shape[0]
                    if est_count_b > gt_count_b:
                        non_perfect_est_more_count += 1
                    elif est_count_b < gt_count_b:
                        non_perfect_est_less_count += 1
                    else:
                        non_perfect_est_equal_count += 1
                
                if n_matched > 0:
                    snr_rmse_data.append({
                        'snr': snr_db, 'rmse_tau': rmse_tau, 'rmse_nu': rmse_nu, 'n_matched': n_matched
                    })

                    # 更新最差 Delay RMSE 样本（bin 与 m）
                    if np.isfinite(rmse_tau):
                        rmse_tau_m = rmse_tau * distance_res
                        if (worst_delay_sample is None) or (rmse_tau > worst_delay_sample['rmse_tau_bin']):
                            worst_delay_sample = {
                                'sample_idx': sample_count,
                                'snr_db': snr_db,
                                'f1': metrics['F1'],
                                'tp': metrics['TP'],
                                'fp': metrics['FP'],
                                'fn': metrics['FN'],
                                'est_count': int(est_params_b.shape[0]),
                                'gt_count': int(gt_labels_b.shape[0]),
                                'n_matched': int(n_matched),
                                'rmse_tau_bin': float(rmse_tau),
                                'rmse_tau_m': float(rmse_tau_m),
                            }

                    # 更新最差 Doppler RMSE 样本（bin 与 m/s）
                    if np.isfinite(rmse_nu):
                        rmse_nu_mps = rmse_nu * v_res
                        if (worst_doppler_sample is None) or (rmse_nu > worst_doppler_sample['rmse_nu_bin']):
                            worst_doppler_sample = {
                                'sample_idx': sample_count,
                                'snr_db': snr_db,
                                'f1': metrics['F1'],
                                'tp': metrics['TP'],
                                'fp': metrics['FP'],
                                'fn': metrics['FN'],
                                'est_count': int(est_params_b.shape[0]),
                                'gt_count': int(gt_labels_b.shape[0]),
                                'n_matched': int(n_matched),
                                'rmse_nu_bin': float(rmse_nu),
                                'rmse_nu_mps': float(rmse_nu_mps),
                            }
                
                est_params_b_np = est_params_b.detach().cpu().numpy() if torch.is_tensor(est_params_b) else est_params_b
                gt_params_plot = gt_params_b_all.detach().cpu().numpy() if torch.is_tensor(gt_params_b_all) else gt_params_b_all.copy()
                    
                # 【修改核心】：将完整的 4维 参数传入，以便后续在表格绘制复数增益
                update_worst_by_exact_snr(
                    worst_by_exact_snr, snr_db, sample_count, metrics,
                    est_params_b_np,
                    gt_params_plot
                )
                update_best_by_exact_snr(
                    best_by_exact_snr, snr_db, sample_count, metrics,
                    est_params_b_np,
                    gt_params_plot
                )
                collect_samples_by_exact_snr(
                    samples_by_exact_snr, snr_db, sample_count, metrics,
                    est_params_b_np,
                    gt_params_plot
                )
                
                sample_count += 1
                
            if sample_count % BATCH_SIZE == 0 or sample_count == len(test_dataset):
                progress = (sample_count / len(test_dataset)) * 100.0
                print(f"[PROGRESS] {sample_count:5d}/{len(test_dataset):5d} ({progress:6.2f}%) - Time elapsed: {time.time()-t_start:.1f}s")
                
    # --- 全局统计代码 ---
    total_tp = sum([m['TP'] for m in all_metrics])
    total_fp = sum([m['FP'] for m in all_metrics])
    total_fn = sum([m['FN'] for m in all_metrics])
    
    overall_precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
    overall_recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0.0
    overall_f1 = 2 * overall_precision * overall_recall / (overall_precision + overall_recall) if (overall_precision + overall_recall) > 0 else 0.0
    
    print(f"\n📊 **整体性能指标**")
    print(f"  总检测数 (TP):       {total_tp}")
    print(f"  整体精度 (Precision): {overall_precision:.4f}")
    print(f"  整体召回 (Recall):    {overall_recall:.4f}")
    print(f"  整体 F1-Score:       {overall_f1:.4f}")
    
    perfect_count = len([m for m in all_metrics if m['F1'] == 1.0])
    print(f"\n🌟 **完美检测样本统计 (F1 = 1.0)**")
    print(f"  完全一致的样本数: {perfect_count} 个 (占比: {(perfect_count/sample_count)*100:.2f}%)")

    non_perfect_count = sample_count - perfect_count
    print(f"\n🔎 **非完美样本计数关系统计 (F1 != 1.0)**")
    print(f"  非完美样本总数: {non_perfect_count} 个")
    print(f"  估计数 > GT数 的样本数: {non_perfect_est_more_count} 个")
    print(f"  估计数 < GT数 的样本数: {non_perfect_est_less_count} 个")
    print(f"  估计数 = GT数 的样本数: {non_perfect_est_equal_count} 个")

    print(f"\n⚠️ **Delay / Doppler 最差样本**")
    if worst_delay_sample is not None:
        print(
            f"  [Worst Delay] sample={worst_delay_sample['sample_idx']}, "
            f"SNR={worst_delay_sample['snr_db']:.2f}dB, F1={worst_delay_sample['f1']:.4f}, "
            f"RMSE={worst_delay_sample['rmse_tau_bin']:.6f} bin ({worst_delay_sample['rmse_tau_m']:.6f} m), "
            f"TP/FP/FN={worst_delay_sample['tp']}/{worst_delay_sample['fp']}/{worst_delay_sample['fn']}, "
            f"K_est/K_gt={worst_delay_sample['est_count']}/{worst_delay_sample['gt_count']}, "
            f"matched={worst_delay_sample['n_matched']}"
        )
    else:
        print("  [Worst Delay] 无可用样本（可能没有任何匹配对）。")

    if worst_doppler_sample is not None:
        print(
            f"  [Worst Doppler] sample={worst_doppler_sample['sample_idx']}, "
            f"SNR={worst_doppler_sample['snr_db']:.2f}dB, F1={worst_doppler_sample['f1']:.4f}, "
            f"RMSE={worst_doppler_sample['rmse_nu_bin']:.6f} bin ({worst_doppler_sample['rmse_nu_mps']:.6f} m/s), "
            f"TP/FP/FN={worst_doppler_sample['tp']}/{worst_doppler_sample['fp']}/{worst_doppler_sample['fn']}, "
            f"K_est/K_gt={worst_doppler_sample['est_count']}/{worst_doppler_sample['gt_count']}, "
            f"matched={worst_doppler_sample['n_matched']}"
        )
    else:
        print("  [Worst Doppler] 无可用样本（可能没有任何匹配对）。")

    # =========================================================================
    # 延后统一画图阶段 (高速无头渲染)
    # =========================================================================
    print(f"\n🎨 开始后台渲染并保存分析图像，请稍候...")
    middle_by_exact_snr = get_middle_by_exact_snr(samples_by_exact_snr)
        
    if len(worst_by_exact_snr) > 0:
        for snr_key in sorted(worst_by_exact_snr.keys()):
            info = worst_by_exact_snr[snr_key]
            # 传入的已经是完整的 4 维参数列表
            plot_detection_result(
                info['est_params'], info['gt_params'], info['sample_idx'],
                title_suffix=f"SNR={info['snr_db']:.1f}dB, F1={info['f1']:.2f}",
                save_dir='test_results/worst_samples_by_snr'
            )
            
    # <--- 新增：绘制并保存最佳样本
    if len(best_by_exact_snr) > 0:
        for snr_key in sorted(best_by_exact_snr.keys()):
            info = best_by_exact_snr[snr_key]
            plot_detection_result(
                info['est_params'], info['gt_params'], info['sample_idx'],
                title_suffix=f"Perfect: SNR={info['snr_db']:.1f}dB, F1={info['f1']:.2f}",
                save_dir='test_results/best_samples_by_snr'
            )

    # 新增：绘制并保存中间样本（每个SNR下F1中位数）
    if len(middle_by_exact_snr) > 0:
        for snr_key in sorted(middle_by_exact_snr.keys()):
            info = middle_by_exact_snr[snr_key]
            plot_detection_result(
                info['est_params'], info['gt_params'], info['sample_idx'],
                title_suffix=f"Middle: SNR={info['snr_db']:.1f}dB, F1={info['f1']:.2f}",
                save_dir='test_results/middle_samples_by_snr'
            )

    if len(snr_rmse_data) > 0:
        print_rmse_vs_snr(snr_rmse_data)
        plot_rmse_vs_snr(snr_rmse_data, 'test_results')
        
    print(f"✅ 所有图像渲染完毕，保存在 test_results/ 目录下！")
    print(f"总耗时: {time.time() - t_start:.2f} 秒")
