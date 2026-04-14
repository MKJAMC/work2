import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np


# =========================================================================
# 1. 可微物理层 (Differentiable Physics Layer)
# =========================================================================
class DifferentiableOTFSGenerator(nn.Module):
    def __init__(self, N, M):
        super().__init__()
        self.N = N
        self.M = M

        nu_range = torch.arange(N).float()
        tau_range = torch.arange(M).float()
        grid_tau, grid_nu = torch.meshgrid(tau_range, nu_range, indexing='ij')
        self.register_buffer('grid_tau', grid_tau.unsqueeze(0).unsqueeze(0))
        self.register_buffer('grid_nu', grid_nu.unsqueeze(0).unsqueeze(0))

    def forward(self, params, x_tf_input):
        """
        params: [B, Targets, 4] -> (nu, tau, amp, phase)
        x_tf_input: [B, M, N]
        """
        B, Num_Targets, _ = params.shape

        p_nu = params[:, :, 0].view(B, Num_Targets, 1, 1)
        p_tau = params[:, :, 1].view(B, Num_Targets, 1, 1)
        p_amp = params[:, :, 2].view(B, Num_Targets, 1, 1)
        p_pha = params[:, :, 3].view(B, Num_Targets, 1, 1)

        # --- (A) Doppler Kernel ---
        delta_nu = self.grid_nu - p_nu
        pi_delta_nu = np.pi * delta_nu
        den_nu = self.N * torch.sin(pi_delta_nu / self.N)
        mask_nu = torch.abs(delta_nu) < 1e-4
        den_nu_safe = torch.where(mask_nu, torch.ones_like(den_nu), den_nu)
        mag_doppler = torch.sin(pi_delta_nu) / den_nu_safe
        mag_doppler = torch.where(mask_nu, torch.ones_like(mag_doppler), mag_doppler)
        phase_term_doppler = torch.exp(-1j * np.pi * (self.N - 1) * delta_nu / self.N)
        term_doppler = mag_doppler * phase_term_doppler

        # --- (B) Delay Kernel ---
        delta_tau = self.grid_tau - p_tau
        pi_delta_tau = np.pi * delta_tau
        den_tau = self.M * torch.sin(pi_delta_tau / self.M)
        mask_tau = torch.abs(delta_tau) < 1e-4
        den_tau_safe = torch.where(mask_tau, torch.ones_like(den_tau), den_tau)
        mag_delay = torch.sin(pi_delta_tau) / den_tau_safe
        mag_delay = torch.where(mask_tau, torch.ones_like(mag_delay), mag_delay)
        phase_term_delay = torch.exp(-1j * np.pi * (self.M - 1) * delta_tau / self.M)
        term_delay = mag_delay * phase_term_delay

        # (C) 额外相位
        phase_offset = torch.exp(-1j * 2 * np.pi * p_tau * p_nu / (self.M * self.N))

        # (D) 幅度 + 相位
        rot_vec = torch.exp(1j * p_pha)

        # 合成 H_DD
        h_comp = p_amp * rot_vec * term_doppler * term_delay * phase_offset
        H_DD = torch.sum(h_comp, dim=1)

        # (E) 卷积 (频域乘法)
        H_TF = torch.fft.fft2(H_DD)
        Y_TF = H_TF * x_tf_input
        Y_recon_complex = torch.fft.ifft2(Y_TF)

        return torch.stack([Y_recon_complex.real, Y_recon_complex.imag], dim=1).float()


# =========================================================================
# 2. 检测网络
# =========================================================================
class Unfolded_OTFS_Stage(nn.Module):
    def __init__(self):
        super().__init__()

        self.net = nn.Sequential(
            nn.Conv2d(5, 32, kernel_size=3, padding=1, padding_mode='circular'),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.Dropout2d(p=0.2),

            nn.Conv2d(32, 64, kernel_size=3, padding=2, dilation=2, padding_mode='circular'),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.Dropout2d(p=0.2),

            nn.Conv2d(64, 32, kernel_size=3, padding=1, padding_mode='circular'),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),

            nn.Conv2d(32, 3, kernel_size=3, padding=1)
        )

    def forward(self, y_orig, y_recon_prev, prob_prev):
        x = torch.cat([y_orig, y_recon_prev, prob_prev], dim=1)
        out = self.net(x)

        prob_new = torch.sigmoid(out[:, 0:1, :, :])
        delta_nu = torch.tanh(out[:, 1:2, :, :]) * 0.5
        delta_tau = torch.tanh(out[:, 2:3, :, :]) * 0.5

        return prob_new, delta_nu, delta_tau


# =========================================================================
# 3. 普适 OTFS 检测器
# =========================================================================
class Universal_OTFS_Detector(nn.Module):
    def __init__(self, N=32, M=64):
        super().__init__()
        self.N = N
        self.M = M

        self.physics_layer = DifferentiableOTFSGenerator(N, M)
        self.unfolded_stage = Unfolded_OTFS_Stage()

        tau_idx = torch.arange(M)
        nu_idx = torch.arange(N)
        tau_grid, nu_grid = torch.meshgrid(tau_idx, nu_idx, indexing='ij')
        self.register_buffer('tau_grid', tau_grid.view(1, 1, M, N))
        self.register_buffer('nu_grid', nu_grid.view(1, 1, M, N))

    def forward(self, y_orig, y_recon_prev, prob_prev):
        return self.unfolded_stage(y_orig, y_recon_prev, prob_prev)

    # ---------------------------------------------------------------------
    # 基础工具
    # ---------------------------------------------------------------------
    def _circular_distance(self, p1, p2):
        d_nu = abs(p1[0] - p2[0])
        d_nu = min(d_nu, self.N - d_nu)
        d_tau = abs(p1[1] - p2[1])
        d_tau = min(d_tau, self.M - d_tau)
        return float(np.sqrt(d_nu ** 2 + d_tau ** 2))


    #触发是否结束
    def _gain_ratio_exceeds(self, points, ratio_thresh=15.0):
        if len(points) <= 1:
            return False
        amps = [max(p[2], 1e-9) for p in points]
        return (max(amps) / min(amps)) > ratio_thresh

    def _collect_all_points(self, groups):
        all_points = []
        for g in groups:
            all_points.extend(g['points'])
        return all_points

    def _find_nearest_group(self, candidate_point, groups):
        if len(groups) == 0:
            return float('inf'), -1

        min_dist = float('inf')
        group_idx = -1
        for gi, g in enumerate(groups):
            for p in g['points']:
                dist = self._circular_distance(candidate_point, p)
                if dist < min_dist:
                    min_dist = dist
                    group_idx = gi
        return min_dist, group_idx

    def _is_physically_unreasonable_candidate(
        self,
        candidate_point,
        group_points,
        close_dist_thresh=1.8,
        weak_amp_ratio=0.28,
    ):
        if len(group_points) == 0:
            return False

        cand_amp = max(candidate_point[2], 1e-9)
        max_group_amp = max(max(p[2], 1e-9) for p in group_points)
        min_dist = min(self._circular_distance(candidate_point, p) for p in group_points)

        # =======================================================
        # 【新增：第一层防御】贴身极近区绞杀 (< 0.45 bin)
        # 这么近的距离，如果能量不到主峰的 50%，绝对是畸变过拟合！
        # =======================================================
        if min_dist < 0.45 and cand_amp < max_group_amp * 0.50:
            return True

        # =======================================================
        # 【原版：第二层防御】常规波纹拦截区 (< 1.8 bin)
        # =======================================================
        if min_dist < close_dist_thresh and cand_amp < max_group_amp * weak_amp_ratio:
            return True

        return False

    # ---------------------------------------------------------------------
    # 从残差中找一个候选峰
    # 将整数网格加上小数偏移量，得到高精度的连续坐标 rough_nu 和 rough_tau。这里的取余操作非常关键，它完美契合了 OTFS 信号在 Delay-Doppler 域特有的二维循环卷积物理边界条件（即目标从多普勒/时延一端出去，会从另一端绕回来）。
    # ---------------------------------------------------------------------
    def _find_peak_from_residual(self, y_residual, prob_threshold):
        device = y_residual.device
        y_zero_recon = torch.zeros_like(y_residual)
        prob_zero = torch.zeros((1, 1, self.M, self.N), device=device, dtype=y_residual.dtype)

        with torch.no_grad():
            prob_map, delta_nu, delta_tau = self.unfolded_stage(
                y_residual, y_zero_recon, prob_zero
            )

        peak_val = torch.max(prob_map).item()
        if peak_val < prob_threshold:
            return None

        flat_idx = torch.argmax(prob_map)
        tau_int = (flat_idx // self.N).item()
        nu_int = (flat_idx % self.N).item()

        nu_frac = delta_nu[0, 0, tau_int, nu_int].item()
        tau_frac = delta_tau[0, 0, tau_int, nu_int].item()

        rough_nu = (nu_int + nu_frac) % self.N
        rough_tau = (tau_int + tau_frac) % self.M

        return {
            'nu': rough_nu,
            'tau': rough_tau,
            'prob': peak_val
        }

    #从残差中估计最可能的点
    def _estimate_candidate_from_residual(self, y_residual, x_tf_input, prob_threshold):
        peak = self._find_peak_from_residual(y_residual, prob_threshold)
        if peak is None:
            return None

        rough_pos = y_residual.new_empty((1, 1, 2))
        rough_pos[0, 0, 0] = peak['nu']
        rough_pos[0, 0, 1] = peak['tau']

        h_opt, _ = self.solve_least_squares_gain(y_residual, x_tf_input, rough_pos)
        amp = h_opt.abs().item()
        phase = h_opt.angle().item()

        return [peak['nu'], peak['tau'], amp, phase]

    # ---------------------------------------------------------------------
    # LS 求解，在已知目标（或候选点）二维物理坐标（多普勒 $nu$、时延 $tau$）的前提下，利用发送信号 $X$ 和接收残差 $Y$，通过严格的物理方程推导出这几个目标到底具有多大的反射面（幅度 Amplitude）和经历了怎样的电磁波相移（相位 Phase）
    # ---------------------------------------------------------------------
    def solve_least_squares_gain(self, y_input, x_tf_input, est_pos, window_size=3):
        B, P, _ = est_pos.shape
        device = y_input.device
        y_complex = torch.complex(y_input[:, 0], y_input[:, 1])

        nu_c = torch.round(est_pos[:, :, 0]).view(B, P, 1, 1)
        tau_c = torch.round(est_pos[:, :, 1]).view(B, P, 1, 1)

        delta_nu = torch.abs(self.nu_grid - nu_c)
        delta_nu = torch.min(delta_nu, self.N - delta_nu)
        delta_tau = torch.abs(self.tau_grid - tau_c)
        delta_tau = torch.min(delta_tau, self.M - delta_tau)
        mask_bool = (delta_nu <= window_size) & (delta_tau <= window_size)
        mask = mask_bool.any(dim=1).float() # [B, M, N]

        ones_amp = torch.ones(B, P, 1, 1, device=device)
        zeros_ph = torch.zeros(B, P, 1, 1, device=device)
        pos_reshaped = est_pos.view(B, P, 2, 1)
        basis_params = torch.cat([pos_reshaped, ones_amp, zeros_ph], dim=2).squeeze(-1)

        basis_params_flat = basis_params.view(B * P, 1, 4)
        x_tf_flat = x_tf_input[:, None, :, :].expand(B, P, self.M, self.N).reshape(B * P, self.M, self.N)

        # 依然用 FFT 算出完整的波形 (因为频域乘法必须是完整维度)
        echo_split = self.physics_layer(basis_params_flat, x_tf_flat)
        echo_complex = torch.complex(echo_split[:, 0], echo_split[:, 1])
        echo_complex = echo_complex.view(B, P, self.M, self.N)

        # =================================================================
        # 🚀 物理先验截断 (Delay Truncation Optimization)
        # 因为目标 delay < 10，我们只截取前 15 个 bin (留 5 个 bin 给旁瓣)
        # =================================================================
        M_trunc = 15
        
        # 裁剪掩码、接收信号、和基向量
        mask_trunc = mask[:, :M_trunc, :]              # [B, 15, N]
        mask_vec = mask_trunc.reshape(B, -1, 1)        # [B, 480, 1]
        
        y_trunc = y_complex[:, :M_trunc, :]            # [B, 15, N]
        y_vec = y_trunc.reshape(B, -1, 1) * mask_vec   # [B, 480, 1]
        
        echo_trunc = echo_complex[:, :, :M_trunc, :]   # [B, P, 15, N]
        
        # 构建极速版 A 矩阵
        A = echo_trunc.reshape(B, P, -1).transpose(1, 2) * mask_vec # [B, 480, P]
        A_H = A.conj().transpose(1, 2)                              # [B, P, 480]
        Gram = torch.bmm(A_H, A)

        if P > 1:
            with torch.no_grad():
                # 使用 .detach() 获取不带梯度的位置张量
                pos_det = est_pos.detach() 
                
                d_nu = torch.abs(pos_det[:, :, 0:1] - pos_det[:, :, 0:1].transpose(1, 2))
                d_nu = torch.min(d_nu, self.N - d_nu)
                
                d_tau = torch.abs(pos_det[:, :, 1:2] - pos_det[:, :, 1:2].transpose(1, 2))
                d_tau = torch.min(d_tau, self.M - d_tau)
                
                dist_matrix = torch.sqrt(d_nu ** 2 + d_tau ** 2)
                
                # 【关键优化 2】预先在外层缓存，或使用对角线提取法，避免生成全新的 torch.eye
                dist_matrix.diagonal(dim1=1, dim2=2).fill_(float('inf'))
                
                min_dist = torch.min(dist_matrix.view(B, -1), dim=1)[0]  # shape: [B]
                
                base_lambda = 1e-4
                max_extra_lambda = 0.025  # 从 0.05 剧烈降低到 0.002
                critical_dist = 0.8      
                
                boost_factor = torch.relu(critical_dist - min_dist) / critical_dist
                dynamic_lambda = base_lambda + boost_factor * max_extra_lambda
        else:
            dynamic_lambda = torch.ones(B, device=device) * 1e-4

        damping = torch.eye(P, device=device).unsqueeze(0).expand(B, -1, -1)
        damping = damping * dynamic_lambda.view(B, 1, 1)

        Gram_Reg = Gram + damping
        RHS = torch.bmm(A_H, y_vec)
        h_opt = torch.linalg.solve(Gram_Reg, RHS)

        h_opt_detached = h_opt.detach()                      # shape: [B, P, 1]
        h_opt_4d = h_opt_detached.unsqueeze(-1)              # shape: [B, P, 1, 1]
        
        # 4D 完美相乘后沿 P(多径) 维度求和，得到全局重构回波
        y_recon_complex_4d = torch.sum(echo_complex * h_opt_4d, dim=1) # shape: [B, M, N]
        
        # 分离实部和虚部
        y_recon_opt = torch.stack(
            [y_recon_complex_4d.real, y_recon_complex_4d.imag], dim=1
        ) # shape: [B, 2, M, N]

        return h_opt.squeeze(-1), y_recon_opt

    # ---------------------------------------------------------------------
    # Adam 精修
    # ---------------------------------------------------------------------
    def _adam_optimize(
    self,
    initial_params,
    target_y,
    x_tf_input,
    num_steps=8,
    lr=0.01,
    ls_update_every=2,
    early_stop_after=2,
    stop_thresh=8e-4,
):
        if initial_params.shape[1] == 0:
            return initial_params

        # 关键：即使外层 inference 在 no_grad() 里，这里也强制恢复梯度
        with torch.enable_grad():
            params_pos = initial_params[:, :, 0:2].clone().detach()
            params_pos.requires_grad_(True)

            optimizer = torch.optim.Adam([params_pos], lr=lr)
            prev_pos = params_pos.clone().detach()
            cached_h_opt = None

            for step in range(num_steps):
                optimizer.zero_grad()

                need_refresh_ls = (
                    step == 0
                    or cached_h_opt is None
                    or (step % ls_update_every == 0)
                )

                if need_refresh_ls:
                    h_opt, y_recon_opt = self.solve_least_squares_gain(
                        target_y, x_tf_input, params_pos
                    )
                    cached_h_opt = h_opt.detach()
                else:
                    ones_amp = cached_h_opt.abs().unsqueeze(-1)
                    phases = cached_h_opt.angle().unsqueeze(-1)
                    full_params = torch.cat([params_pos, ones_amp, phases], dim=-1)
                    y_recon_opt = self.physics_layer(full_params, x_tf_input)

                loss = F.mse_loss(y_recon_opt, target_y)

                if params_pos.shape[1] > 1:
                    P = params_pos.shape[1]
                    repulsive_terms = []
                    for i in range(P):
                        for j in range(i + 1, P):
                            d_nu = torch.abs(params_pos[:, i, 0] - params_pos[:, j, 0])
                            d_nu = torch.min(d_nu, self.N - d_nu)
                            d_tau = torch.abs(params_pos[:, i, 1] - params_pos[:, j, 1])
                            d_tau = torch.min(d_tau, self.M - d_tau)
                            dist = torch.sqrt(d_nu ** 2 + d_tau ** 2 + 1e-8)
                            repulsive_terms.append(F.relu(0.15 - dist))

                    if len(repulsive_terms) > 0:
                        repulsive_penalty = torch.mean(torch.stack(repulsive_terms, dim=0))
                        loss = loss + 0.03 * repulsive_penalty

                loss.backward()
                optimizer.step()

                with torch.no_grad():
                    params_pos[:, :, 0] = params_pos[:, :, 0] % self.N
                    params_pos[:, :, 1] = params_pos[:, :, 1] % self.M

                    max_displacement = torch.max(
                        torch.sqrt(torch.sum((params_pos - prev_pos) ** 2, dim=-1))
                    )

                    if step > early_stop_after and max_displacement < stop_thresh:
                        break

                    prev_pos = params_pos.clone().detach()

            with torch.no_grad():
                h_opt_final, _ = self.solve_least_squares_gain(target_y, x_tf_input, params_pos)
                final_amp = h_opt_final.abs().unsqueeze(-1)
                final_phase = h_opt_final.angle().unsqueeze(-1)
                final_params = torch.cat([params_pos, final_amp, final_phase], dim=-1)

        return final_params.detach()

    def _single_refine(self, point, target_y, x_tf_input):
        init = target_y.new_tensor([[[point[0], point[1]]]])
        refined = self._adam_optimize(
            init, target_y, x_tf_input,
            num_steps=15,            # 从 8 提升到 15
            lr=0.02,                 # 从 0.01 提升到 0.02
            ls_update_every=2,
            early_stop_after=3,
            stop_thresh=5e-4,
        )
        return refined[0, 0].tolist()

    def _joint_refine(self, points, target_y, x_tf_input):
        init = target_y.new_tensor([[p[:2] for p in points]])
        refined = self._adam_optimize(
            init, target_y, x_tf_input,
            num_steps=30,            # 【核心修改 2】从 15 步暴增到 40 步！给它足够的时间汇聚
            lr=0.01,                 # 【核心修改 2】学习率翻 4 倍，赋予强大的拉扯力
            ls_update_every=2,
            early_stop_after=8,
            stop_thresh=5e-4,
        )
        return refined[0].tolist()
    
    def _active_fission_test(self, candidate_point, target_y, x_tf_input):
        """
        主动裂变探测 (极致灵敏 + 宽进严出版)
        核心逻辑：
        1. 绝不盲目启动 Adam。先用无梯度的 Least Squares (LS) 进行极速探路。
        2. 初审 (宽进)：允许 MSE 仅下降 5% 的双点模型通过，防止漏掉间距极近的真实目标。
        3. 终审 (严出)：利用 Adam 物理推演，严酷绞杀距离小于 0.3 或幅值被榨干的过拟合单点。
        """
        nu_c, tau_c = candidate_point[0], candidate_point[1]

        # 1. 正常单点精修作为 Baseline (仅调用一次 Adam)
        refined_single = self._single_refine(candidate_point, target_y, x_tf_input)
        nu_r, tau_r = refined_single[0], refined_single[1]

        # =====================================================================
        # 极速探路阶段 (无梯度，纯矩阵运算，速度极快)
        # 测试四种经典干涉裂变方向（主对角、副对角、水平、垂直）
        # =====================================================================
        # 探路偏移量设为 0.25，生成的种子间距约为 0.5~0.7 bin，完美契合极近多径盲区
        offset = 0.15
        
        # 打包送入 LS 评估
        pos_S = target_y.new_tensor([[[nu_r, tau_r]]])
        
        seeds = [
            [[nu_r + offset, tau_r + offset], [nu_r - offset, tau_r - offset]], # 主对角
            [[nu_r + offset, tau_r - offset], [nu_r - offset, tau_r + offset]], # 副对角
            [[nu_r + offset, tau_r],          [nu_r - offset, tau_r]],          # 水平
            [[nu_r, tau_r + offset],          [nu_r, tau_r - offset]]           # 垂直
        ]

        best_loss = float('inf')
        best_seeds = None

        with torch.no_grad():
            # 测单点 Baseline 残差
            _, y_recon_S = self.solve_least_squares_gain(target_y, x_tf_input, pos_S)
            loss_S = F.mse_loss(y_recon_S, target_y).item()
            best_loss = loss_S

            # 遍历四种双点方向，寻找 MSE 最小的物理路径
            for seed_pair in seeds:
                pos_pair = target_y.new_tensor([seed_pair])
                h_pair, y_recon_pair = self.solve_least_squares_gain(target_y, x_tf_input, pos_pair)
                loss_pair = F.mse_loss(y_recon_pair, target_y).item()

                # 初审阶段幅值健康度检查：放宽至 0.05，因为静态位置未优化，允许出现轻微畸形
                amp1, amp2 = h_pair[0, 0].abs().item(), h_pair[0, 1].abs().item()
                if min(amp1, amp2) > 0.05 * max(amp1, amp2):
                    if loss_pair < best_loss:
                        best_loss = loss_pair
                        best_seeds = seed_pair

        # 【终极奥卡姆剃刀】：门限从 0.05 提升到 0.12
        # 物理意义：想多用 4 个参数把一个峰劈成两个？可以！
        # 但你必须证明这对残差有实质性的、超过 12% 的巨大改善。
        # 否则你就是在“均分能量过拟合”，直接打回单点！
        if best_seeds is None or loss_S == 0 or (loss_S - best_loss) / loss_S < 0.12:
            return [refined_single]

        # =====================================================================
        # 确认疑似双径，启动唯一一次定向 Adam 精修
        # =====================================================================
        init_pair = [
            [best_seeds[0][0], best_seeds[0][1], 0, 0], 
            [best_seeds[1][0], best_seeds[1][1], 0, 0]
        ]
        
        # 只在这里对拿到“准考证”的双径跑一次 Adam
        refined_pair = self._joint_refine(init_pair, target_y, x_tf_input)

        # =====================================================================
        # 终审判决 (极其严格的物理绞杀)
        # =====================================================================
        amp1, amp2 = refined_pair[0][2], refined_pair[1][2]
        dist = self._circular_distance(refined_pair[0], refined_pair[1])
        
        min_a = min(amp1, amp2)
        max_a = max(amp1, amp2)

        # 【新增：阶梯式终审】
        if dist <= 0.25:
            return [refined_single] # 距离太近，强行打回
            
        elif dist < 0.45:
            # 贴身极近距：必须具备 50% 能量才能被承认为隐藏双径
            if min_a > 0.50 * max_a:
                return refined_pair
            else:
                return [refined_single]
                
        else:
            # 常规距离：维持 25% 的能量底线
            if min_a > 0.25 * max_a:
                return refined_pair
            else:
                return [refined_single]
        
        # 【核心逻辑】：如果是模型自由度过高导致的单点过拟合（幽灵峰），
        # Adam 在收敛后，要么把两个点吸回单点（dist < 0.3），
        # 要么通过反相抵消把其中一个点的有效幅值降到极低（< 15%）。
        if min(amp1, amp2) > 0.25 * max(amp1, amp2) and dist > 0.25:
            return refined_pair  # 经受住了 Adam 的终极考验，确认为隐藏双径！
        else:
            return [refined_single] # 伪装被物理法则识破，打回单点！
        
    # ---------------------------------------------------------------------
    # 用当前所有点做一次全局回代，更新每个点的 amp/phase，并得到新残差
    # ---------------------------------------------------------------------
    def _global_refit_and_update_residual(self, y_orig, x_tf_input, groups):
        all_points = self._collect_all_points(groups)

        if len(all_points) == 0:
            y_recon = torch.zeros_like(y_orig)
            y_residual = y_orig.clone()
            return y_recon, y_residual, []

        all_pos = y_orig.new_tensor([[p[:2] for p in all_points]])
        h_opt_global, y_recon_all = self.solve_least_squares_gain(y_orig, x_tf_input, all_pos)

        h_opt_np = h_opt_global[0].detach().cpu().numpy()

        idx = 0
        for g in groups:
            for p in g['points']:
                p[2] = float(np.abs(h_opt_np[idx]))
                p[3] = float(np.angle(h_opt_np[idx]))
                idx += 1

        y_residual = y_orig - y_recon_all
        updated_points = self._collect_all_points(groups)

        return y_recon_all, y_residual, updated_points

    # ---------------------------------------------------------------------
    # 主推理：按你描述的多目标流程
    # ---------------------------------------------------------------------
    # 【终极完美版】: VarPro 前置去干扰 + T3专属物理微调 + 解除位置枷锁
    def inference(
        self,
        y_input,
        x_input,
        max_targets=12,
        prob_threshold=0.03,
        close_bin_threshold=1.8,
        stop_gain_ratio=15.0,
    ):
        """
        多目标检测流程：
        1) 网络从残差中逐次找最强目标 T1, T2, ...
        2) 远目标：视为单目标，单独 Adam 精修（<10 次）
        3) 近目标：放入同一集合，先重构回波并更新残差
        4) 若新残差里出现的候选点与该集合“距离近但幅值弱”，视作非新目标，
           说明该近邻集合未估准，对该集合做联合精修
        5) 重复上述过程
        6) 若当前目标集合 max_gain / min_gain > 5，则停止
        """
        B = y_input.shape[0]
        device = y_input.device
        x_complex = torch.complex(x_input[:, 0], x_input[:, 1])
        X_TF_Global = torch.fft.fft2(x_complex)

        final_params_list = [[] for _ in range(B)]

        for b in range(B):#对于batchsize中的每个样本进行处理
            y_orig_b = y_input[b:b+1]
            x_tf_b = X_TF_Global[b:b+1]
            y_residual_b = y_orig_b.clone()

            groups = []
            accepted_count = 0

            while accepted_count < max_targets:
                # ----------------------------------------------------------
                # Step 1: 从当前残差中找最强候选目标
                # ----------------------------------------------------------
                candidate = self._estimate_candidate_from_residual(
                    y_residual_b, x_tf_b, prob_threshold
                )
                if candidate is None:
                    break

                # ----------------------------------------------------------
                # Step 2: 如果是第一个目标，直接单目标精修
                # ----------------------------------------------------------
                if len(groups) == 0:
                    fission_points = self._active_fission_test(candidate, y_residual_b, x_tf_b)# 主动裂变探测器
                    groups.append({
                        'points': fission_points, # 可能是1个点，也可能是裂变后的2个点
                        'is_close_group': False
                    })

                    _, y_residual_b, all_points = self._global_refit_and_update_residual(
                        y_orig_b, x_tf_b, groups
                    )

                    accepted_count += len(fission_points) # 动态增加计数值

                    if self._gain_ratio_exceeds(all_points, stop_gain_ratio):
                        break

                    continue


                # ----------------------------------------------------------
                # Step 3: 判断与已有目标/集合的距离
                # ----------------------------------------------------------
                min_dist, group_idx = self._find_nearest_group(candidate, groups)

                # 远：说明它是新单目标
                if min_dist > close_bin_threshold:
                    fission_points = self._active_fission_test(candidate, y_residual_b, x_tf_b)
                    groups.append({
                        'points': fission_points,
                        'is_close_group': False
                    })

                    _, y_residual_b, all_points = self._global_refit_and_update_residual(
                        y_orig_b, x_tf_b, groups
                    )
                    accepted_count += len(fission_points)

                    if self._gain_ratio_exceeds(all_points, stop_gain_ratio):
                        break

                    continue


                # ----------------------------------------------------------
                # Step 4: 近目标：拦截旁瓣 OR 放入集合
                # ----------------------------------------------------------
                # 【补上丢失的防御阵地：局部旁瓣拦截器】
                unreasonable_cand = self._is_physically_unreasonable_candidate(
                    candidate,
                    groups[group_idx]['points'],
                    close_dist_thresh=close_bin_threshold,
                    weak_amp_ratio=0.28, # FA4 只有 18.6%，会被这里直接斩杀！T2 有 54%，安全通过！
                )

                if unreasonable_cand:
                    # 发现旁瓣！拒绝加入集合。反手对现有集合做一次联合精修来消除波纹
                    groups[group_idx]['points'] = self._joint_refine(groups[group_idx]['points'], y_orig_b, x_tf_b)
                    _, y_residual_b, all_points = self._global_refit_and_update_residual(y_orig_b, x_tf_b, groups)
                    
                    if self._gain_ratio_exceeds(all_points, stop_gain_ratio):
                        break
                    continue 

                # 如果真的是合法近距目标 (比如 T2)，才允许加入集合
                groups[group_idx]['points'].append(candidate)
                groups[group_idx]['is_close_group'] = True

                # 先用当前集合重构一次，得到新残差
                _, y_residual_b, all_points = self._global_refit_and_update_residual(
                    y_orig_b, x_tf_b, groups
                )
                accepted_count += 1          

                
                #还是基于近目标的情况下
                # ----------------------------------------------------------
                # Step 5: 探针 (Probe) 测试与策略分流 (核心修复区)
                # ----------------------------------------------------------
                probe_candidate = self._estimate_candidate_from_residual(
                    y_residual_b, x_tf_b, prob_threshold
                )

                if probe_candidate is not None:
                    # 判断这个探针是当前 group 估偏产生的伪峰，还是真实存在的新目标
                    unreasonable = self._is_physically_unreasonable_candidate(
                        probe_candidate,
                        groups[group_idx]['points'],
                        close_dist_thresh=close_bin_threshold,
                        weak_amp_ratio=0.75,
                    )

                    if unreasonable:
                        # 【情况 A：发现伪峰】
                        # 说明刚才的 group 估偏了 -> 立即对该 group 联合精修
                        refined_points = self._joint_refine(groups[group_idx]['points'], y_orig_b, x_tf_b)
                        groups[group_idx]['points'] = refined_points
                        
                    else:
                        # 【情况 B：探针是有效的真目标】
                        min_dist_probe, _ = self._find_nearest_group(probe_candidate, groups)
                        
                        if min_dist_probe > close_bin_threshold:
                            # 按照你的思路：探针是离得远的真目标！
                            # 1. 先对远目标进行单目标精修，并加入组群
                            refined_far = self._single_refine(probe_candidate, y_residual_b, x_tf_b)
                            groups.append({'points': [refined_far], 'is_close_group': False})
                            accepted_count += 1
                            
                            # 2. 【核心回溯】再去对刚才那个没精修的 group 进行联合精修！
                            # 此时远目标的干扰已被剥离，group 精修极其精准
                            refined_group = self._joint_refine(groups[group_idx]['points'], y_orig_b, x_tf_b)
                            groups[group_idx]['points'] = refined_group
                        else:
                            # 【核心修改 3】补全分支：探针是真目标，且离得很近 (比如罕见的紧密三径)
                            groups[group_idx]['points'].append(probe_candidate)
                            refined_group = self._joint_refine(groups[group_idx]['points'], y_orig_b, x_tf_b)
                            groups[group_idx]['points'] = refined_group
                            accepted_count += 1
                            
                else:
                    # 【情况 C：没有新目标了 (残差干净)】
                    # 但刚才放入的 group 还没被联合精修，必须在结束前补上这一步！
                    refined_points = self._joint_refine(groups[group_idx]['points'], y_orig_b, x_tf_b)
                    groups[group_idx]['points'] = refined_points

                # 不管经历了哪种情况，最后统一做一次全局 LS 更新
                _, y_residual_b, all_points = self._global_refit_and_update_residual(
                    y_orig_b, x_tf_b, groups
                )

                if self._gain_ratio_exceeds(all_points, stop_gain_ratio):
                    break

            # --------------------------------------------------------------
            # 最终后处理：按 max/min <= 5 保留强点，再做去重
            # --------------------------------------------------------------
            all_points = self._collect_all_points(groups)

            if len(all_points) > 0:
                max_gain = max(max(p[2], 1e-9) for p in all_points)

                filtered_points = []
                for p in all_points:
                    if max_gain / max(p[2], 1e-9) <= stop_gain_ratio:
                        filtered_points.append(p)

                
                # =========================================================
                # 先按能量从大到小排序，确保强点优先保留，弱点必须接受强点的审视
                filtered_points = sorted(filtered_points, key=lambda x: x[2], reverse=True)
                final_kept = []
                
                for p in filtered_points:
                    is_fake = False
                    for kept_p in final_kept:
                        dist = self._circular_distance(p, kept_p)
                        
                        # 规则1：绝对重合 (< 0.25 bin)，无论多大能量，直接吞噬
                        if dist < 0.25:
                            is_fake = True
                            break
                        # 规则2：贴身畸变区 (0.25 ~ 0.50 bin)
                        # 如果你在强点旁边这么近，但能量连强点的 50% 都不到，你必是 Adam 拉扯出来的假波纹！
                        elif dist < 0.50:
                            if p[2] < 0.50 * kept_p[2]:
                                is_fake = True
                                break
                        # 规则3：常规旁瓣区 (0.50 ~ 0.85 bin)
                        # 如果距离稍远，但能量异常微弱（不到强点 25%），同样按旁瓣抹除
                        elif dist < 0.85:
                            if p[2] < 0.25 * kept_p[2]:
                                is_fake = True
                                break
                                
                    if not is_fake:
                        final_kept.append(p)

                final_params_list[b] = final_kept
            else:
                final_params_list[b] = []

        max_K = max([len(t) for t in final_params_list]) if any(final_params_list) else 0
        if max_K == 0:
            return torch.zeros((B, 1, 4), device=device), y_input, y_input

        output_tensor = torch.zeros((B, max_K, 4), device=device)
        for b, params in enumerate(final_params_list):
            if params:
                output_tensor[b, :len(params), :] = y_input.new_tensor(params)

        y_recon_final = self.physics_layer(output_tensor, X_TF_Global)
        y_residual_final = y_input - y_recon_final

        return output_tensor, y_recon_final, y_residual_final