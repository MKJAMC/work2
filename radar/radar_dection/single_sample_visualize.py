import argparse
import os

import torch

from dataset import OTFSDataset
from model import Universal_OTFS_Detector
from test import compute_detection_and_rmse, plot_detection_result


def _format_complex_gain(amp, phase):
    comp = amp * complex(torch.cos(torch.tensor(phase)).item(), torch.sin(torch.tensor(phase)).item())
    return f"{comp.real:.6f}{comp.imag:+.6f}j"


def _print_params_table(title, params):
    print(f"\n[{title}] K={len(params)}")
    if len(params) == 0:
        print("  (empty)")
        return

    print("  idx |      nu |     tau |      amp |    phase | complex_gain")
    print("  " + "-" * 68)
    for i, p in enumerate(params):
        nu, tau, amp, phase = float(p[0]), float(p[1]), float(p[2]), float(p[3])
        cplx = _format_complex_gain(amp, phase)
        print(f"  {i:3d} | {nu:7.4f} | {tau:7.4f} | {amp:8.5f} | {phase:8.5f} | {cplx}")


def run_one_sample(
    sample_index,
    dataset_path,
    model_path,
    output_dir,
    max_targets,
    prob_threshold,
    close_bin_threshold,
    stop_gain_ratio,
):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    dataset = OTFSDataset(dataset_path, normalize=True)
    if sample_index < 0 or sample_index >= len(dataset):
        raise IndexError(f"sample_index={sample_index} 越界，合法范围 [0, {len(dataset)-1}]")

    y_in, x_in, gt_all, _, _, snr_value = dataset[sample_index]

    y_in = y_in.unsqueeze(0).to(device)
    x_in = x_in.unsqueeze(0).to(device)

    model = Universal_OTFS_Detector(M=64, N=32).to(device)
    model.load_state_dict(torch.load(model_path, map_location=device), strict=False)
    model.eval()

    with torch.no_grad():
        est_params, _, _ = model.inference(
            y_in,
            x_in,
            max_targets=max_targets,
            prob_threshold=prob_threshold,
            close_bin_threshold=close_bin_threshold,
            stop_gain_ratio=stop_gain_ratio,
        )

    est_params_b = est_params[0]
    est_params_b = est_params_b[est_params_b[:, 2] > 1e-5]

    gt_params_b = gt_all[gt_all[:, 2] > 1e-8]

    metrics, rmse_tau, rmse_nu, n_matched, _, _ = compute_detection_and_rmse(
        est_params_b,
        gt_params_b,
        det_threshold=1.0,
        rmse_threshold=1.5,
    )

    os.makedirs(output_dir, exist_ok=True)
    save_path = plot_detection_result(
        est_params_b.detach().cpu().numpy(),
        gt_all.detach().cpu().numpy(),
        sample_index,
        title_suffix=f"SNR={float(snr_value):.1f}dB, F1={metrics['F1']:.2f}",
        save_dir=output_dir,
    )

    print("\n" + "=" * 72)
    print("Single Sample Inference Finished")
    print("=" * 72)
    print(f"Sample Index: {sample_index}")
    print(f"SNR: {float(snr_value):.2f} dB")
    print(f"TP/FP/FN: {metrics['TP']}/{metrics['FP']}/{metrics['FN']}")
    print(f"Precision: {metrics['Precision']:.4f}")
    print(f"Recall: {metrics['Recall']:.4f}")
    print(f"F1: {metrics['F1']:.4f}")
    if n_matched > 0:
        print(f"Delay RMSE: {rmse_tau:.6f} bin")
        print(f"Doppler RMSE: {rmse_nu:.6f} bin")
    else:
        print("Delay/Doppler RMSE: 无匹配目标")
    print(f"Plot saved: {save_path}")

    _print_params_table("GT Params (valid)", gt_params_b.detach().cpu().numpy())
    _print_params_table("Estimated Params (valid)", est_params_b.detach().cpu().numpy())


def parse_args():
    parser = argparse.ArgumentParser(description="输入一个 sample 索引并可视化当前算法结果")
    parser.add_argument("--sample-index", type=int, default=None, help="样本索引，留空则交互输入")
    parser.add_argument("--dataset", type=str, default="data_h5/OTFS_Test_Set_T1T2_1000.h5", help="测试集路径")
    parser.add_argument("--model", type=str, default="universal_otfs_stripper_best.pth", help="模型权重路径")
    parser.add_argument("--output-dir", type=str, default="test_results/single_sample", help="图像输出目录")
    parser.add_argument("--max-targets", type=int, default=6, help="最大目标数")
    parser.add_argument("--prob-threshold", type=float, default=0.04, help="概率阈值")
    parser.add_argument("--close-bin-threshold", type=float, default=1.8, help="近邻距离阈值")
    parser.add_argument("--stop-gain-ratio", type=float, default=7.0, help="增益比停止阈值")
    return parser.parse_args()


def main():
    args = parse_args()

    sample_index = args.sample_index
    if sample_index is None:
        sample_index = int(input("请输入 sample 索引（从 0 开始）: ").strip())

    run_one_sample(
        sample_index=sample_index,
        dataset_path=args.dataset,
        model_path=args.model,
        output_dir=args.output_dir,
        max_targets=args.max_targets,
        prob_threshold=args.prob_threshold,
        close_bin_threshold=args.close_bin_threshold,
        stop_gain_ratio=args.stop_gain_ratio,
    )


if __name__ == "__main__":
    main()
