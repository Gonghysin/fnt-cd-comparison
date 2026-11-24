function y = eq_fdc_fft(r, cfg, z, N)
%EQ_FDC_FFT  方案 B：FFT 频域色散补偿，overlap-save 分块（浮点，M2）。
%   y = eq_fdc_fft(r, cfg, z, N)
%   r  2 SPS 接收序列；z 传输距离（m）；N 块长（FFT 点数）
%
%   分块：块长 N，重叠 P = N*cfg.blk.overlap，每块新增 B = N-P 个样本。
%   每块取输入 [jB-P/2+1 ... jB-P/2+N]，FFT -> 乘 H -> IFFT，只保留中间 B 个
%   （两端各丢 P/2）。这是"零相位"的居中写法，对应 Xing 里丢弃 M_d/2 个样本。
%   H 是理想逆色散在 N 点网格上的采样（浮点，未量化）。
%
%   首尾块需要信号之外的样本，这里按**环绕**取（与 apply_cd 的圆卷积信道一致），
%   所以两端也不会失真；BER 统计仍然丢弃两端。

r = r(:);
L = numel(r);
P = round(N * cfg.blk.overlap);
B = N - P;
assert(mod(P, 2) == 0 && mod(L, B) == 0, '块参数与信号长度不匹配');
nBlk = L / B;

H = cd_response(cfg, -z, N);                         % N x 1

start = (0:nBlk-1) * B - P/2;                        % 每块的起点（0 基）
idx = mod(start + (0:N-1).', L) + 1;                 % N x nBlk 的输入下标，环绕
R = r(idx);                                          % N x nBlk，每列一块

Yb = ifft(fft(R, [], 1) .* H, [], 1);                % 逐列 FFT -> 乘 H -> IFFT
y = reshape(Yb(P/2+1 : P/2+B, :), [], 1);            % 每块保留中间 B 个，按顺序拼回
end
