import h5py
import torch
import numpy as np
from torch.utils.data import Dataset

class OTFSDataset(Dataset):
    def __init__(self, h5_file, M=64, N=32, normalize=True):
        """
        Args:
            h5_file: MATLAB 生成的 .h5 文件路径
            M: Delay grid size (64)
            N: Doppler grid size (32)
            normalize: 是否对 Y 和 X 进行归一化 (推荐 True)
        """
        self.file_path = h5_file
        self.M = M
        self.N = N
        self.normalize = normalize

        # 1. 读取数据
        with h5py.File(self.file_path, 'r') as f:
            # 读取原始数据
            # 注意：h5py 读取的数据通常是 float64，建议转为 float32 节省内存
            Y_real = f['/Y_real'][:].astype(np.float32)
            Y_imag = f['/Y_imag'][:].astype(np.float32)
            X_real = f['/X_real'][:].astype(np.float32)
            X_imag = f['/X_imag'][:].astype(np.float32)
            labels = f['/labels'][:].astype(np.float32)
            
            # 读取SNR（如果存在）
            if 'snr_index' in f.keys():
                snr_data = f['/snr_index'][:].astype(np.float32)
                # shape可能是(1, N)或(N, 1)，展平成一维数组
                self.snr = snr_data.flatten()
            else:
                # 如果没有SNR字段，设置为None或默认值
                self.snr = None

        # 2. 智能维度处理 (Auto-Permute)
        # 目标形状: 
        #   Y/X -> (Samples, M, N)  [Delay=M, Doppler=N]
        #   labels -> (Samples, P, 3)
        
        # 寻找哪个维度是 Samples (通常最大的那个)
        # 假设 Samples >> M, N
        shape_y = Y_real.shape
        sample_dim_idx = np.argmax(shape_y) 
        
        # 调整 Y 和 X
        self.Y_real = self._fix_dims(Y_real, sample_dim_idx, M, N)
        self.Y_imag = self._fix_dims(Y_imag, sample_dim_idx, M, N)
        self.X_real = self._fix_dims(X_real, sample_dim_idx, M, N)
        self.X_imag = self._fix_dims(X_imag, sample_dim_idx, M, N)
        
        # 调整 Labels
        # 我们寻找哪个维度等于 Samples
        if labels.shape[0] == self.Y_real.shape[0]: # Case A: (Samples, ...)
            # 检查剩余维度，哪个是 4 (参数维)
            if labels.shape[1] == 4: # (Samples, 4, P) -> 需要转置
                 self.labels = np.transpose(labels, (0, 2, 1))
            else: # (Samples, P, 4) -> 也就是 shape[2]==4，直接用
                 self.labels = labels
                 
        elif labels.shape[2] == self.Y_real.shape[0]: # Case B: (..., Samples)
            if labels.shape[0] == 4: # (4, P, Samples) -> 转置为 (Samples, P, 4)
                self.labels = np.transpose(labels, (2, 1, 0))
            else: # (P, 4, Samples) -> 转置为 (Samples, P, 4)
                self.labels = np.transpose(labels, (2, 0, 1))
        else:
            # 兜底：如果维度乱了，尝试强制 reshape 或者保持原样
            # 假设现在的 H5 只有 (Samples, P, 4) 这一种标准情况
            self.labels = labels

        # 3. 合成复数并转 Tensor
        self.Y = torch.from_numpy(self.Y_real + 1j * self.Y_imag).cfloat()
        self.X = torch.from_numpy(self.X_real + 1j * self.X_imag).cfloat()
        self.labels = torch.from_numpy(self.labels).float() # 现在是 [Samples, P, 4]
      
        print(f"Dataset loaded successfully.")
        print(f"  - Y Shape: {self.Y.shape} (Batch, Delay, Doppler)")
        print(f"  - X Shape: {self.X.shape}")
        print(f"  - Labels : {self.labels.shape} (Batch, Targets, Params)")

    def _fix_dims(self, arr, sample_idx, M, N):
        """辅助函数：将数组转置为 (Samples, M, N)"""
        if sample_idx == 0: # (Samples, ..., ...)
            # 此时需要判断后面两维哪个是 M
            if arr.shape[1] == M: return arr # 已经是 (Samples, M, N)
            else: return np.transpose(arr, (0, 2, 1)) # (Samples, N, M) -> (Samples, M, N)
            
        elif sample_idx == 2: # (..., ..., Samples)
            # 需要把 Samples 移到最前
            if arr.shape[0] == M: return np.transpose(arr, (2, 0, 1)) # (M, N, Samples) -> (Samples, M, N)
            else: return np.transpose(arr, (2, 1, 0)) # (N, M, Samples) -> (Samples, M, N)
            
        return arr # Fallback

    def __len__(self):
        return self.Y.shape[0]

    def __getitem__(self, idx):
        y = self.Y[idx] # (M, N)
        x = self.X[idx] # (M, N)
        
        # 1. 获取原始标签：Shape [P, 4] -> [tau, nu, amp, phase]
        raw_label = self.labels[idx]  # [3, 4] 三个目标
        
        # 2. 提取各参数（所有目标）
        gt_tau = raw_label[:, 0]      # [P] Delay
        gt_nu = raw_label[:, 1]       # [P] Doppler  
        gt_amp = raw_label[:, 2]      # [P] Amplitude
        gt_phase = raw_label[:, 3]    # [P] Phase (弧度)
        
        # 3. 归一化处理
        if self.normalize:
            scale_y = torch.max(torch.abs(y)) + 1e-8
            y = y / scale_y
            scale_x = torch.max(torch.abs(x)) + 1e-8
            x = x / scale_x
            
            # 幅度缩放（物理关系）
            gt_amp = gt_amp * (scale_x / scale_y)
            # 相位不需要缩放
        else:
            scale_y = torch.tensor(1.0)
            scale_x = torch.tensor(1.0)
        
        # 【关键修复】将物理有符号索引转换为周期索引
        # MATLAB生成的nu是有符号的（-8到+8），需要转为0-31
        gt_nu = gt_nu % self.N  # -2.5 → 29.5, 3.5 → 3.5
        gt_tau = gt_tau % self.M  # 同样处理tau
        
        # 4. 转换为模型期望格式：[nu, tau, amp, phase]
        final_label = torch.stack([
            gt_nu,      # Index 0: Doppler (nu) - 已周期化
            gt_tau,     # Index 1: Delay (tau) - 已周期化
            gt_amp,     # Index 2: Amplitude
            gt_phase    # Index 3: Phase
        ], dim=-1)  # Shape: [P, 4]

        y_input = torch.stack([y.real, y.imag], dim=0).float()
        x_input = torch.stack([x.real, x.imag], dim=0).float()
        
        # 获取SNR（如果有的话）
        if self.snr is not None:
            snr_value = torch.tensor(self.snr[idx], dtype=torch.float32)
        else:
            snr_value = torch.tensor(0.0, dtype=torch.float32)  # 默认值
        
        return y_input, x_input, final_label, scale_y, scale_x, snr_value
        