%% MMSE 检测
heff=H_eff';
Ie=eye(size(H_eff, 2));
H_lmmse=(heff*((H_eff*heff+(1/factor(i_snr))*Ie))^(-1));%平均snr每个符号的
dd_est_matrix=H_lmmse*y;
dd_est_matrix= reshape(dd_est_matrix, M, []);
estimated_bits_stream=[];
