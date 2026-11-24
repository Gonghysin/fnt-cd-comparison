function [y, h] = eq_fir(r, cfg, z, nTap)
%EQ_FIR  方案 A：时域 FIR 静态色散补偿（浮点，M2）。
% 时域 FIR 是把同一个 $H$ 换成时域的等价形式，直接做卷积

%   [y, h] = eq_fir(r, cfg, z, nTap)
%   r     2 SPS 接收序列；z 传输距离（m）
%   nTap  抽头数（缺省取 cfg.td.nTap_Rs 的公式：按符号率带宽估计）；偶数会加 1 变奇数
%   y     补偿后序列（与 r 等长，已对齐）；h 为实际使用的 FIR 系数
%
%   系数来源：理想逆色散的频率响应 exp(+j*beta2*z/2*w^2) 在 2^16 点网格上做
%   IFFT 得到时域冲激响应，取以零时刻为中心的 nTap 个样本（直接截断，不加窗）。
%   这与 Xing 的"IFFT -> 截位"预处理是同一条路，后面 FNT 方案的系数也从这里来。
%
%   滤波用线性卷积，'same' 取中间部分，抽头对称所以无延迟；两端约 nTap/2
%   个样本因缺少前后数据而不准，BER 统计时已丢弃。

if nargin < 4 || isempty(nTap) % 如果没有指定抽头数的话，就用参数里的默认值。
    dl = abs(cfg.fiber.D) * cfg.fiber.lambda^2 * abs(z) / cfg.fiber.c;   % s/Hz
    nTap = ceil(dl * cfg.sig.Rs * cfg.sig.fs) + 1;
end
nTap = nTap + (mod(nTap, 2) == 0);           % 保证奇数，中心抽头对准零时刻

nFFT = 2^16;
assert(nTap < nFFT/2, 'nTap 过大，超过 IFFT 网格');
Ginv = cd_response(cfg, -z, nFFT);           % 理想逆色散
hs = fftshift(ifft(Ginv));                   % 零时刻移到中间（第 nFFT/2+1 个）
c = nFFT/2 + 1;
half = (nTap - 1) / 2;
h = hs(c-half : c+half); % 截断

y = conv(r(:), h, 'same');
end
