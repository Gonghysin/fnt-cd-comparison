function [y, sigma2] = add_noise(x, cfg, EbN0_dB)
%ADD_NOISE  按 Eb/N0 加复高斯白噪声（M1）。
%   [y, sigma2] = add_noise(x, cfg, EbN0_dB)
%
%   定义：噪声是整个采样带宽 fs 上的白噪声，复噪声总方差为 sigma2。
%     符号能量   Es = Ps / Rs            （Ps 为信号平均功率）
%     噪声谱密度 N0 = sigma2 / fs
%     Es/N0 = sps * Ps / sigma2          （因为 fs/Rs = sps）
%     Eb/N0 = Es/N0 / bitsPerSym
%   所以   sigma2 = sps * Ps / (bitsPerSym * 10^(EbN0_dB/10))。
%   这样接收端用单位能量 RRC 匹配滤波后，QPSK 的 BER 应等于 Q(sqrt(2*Eb/N0))。

% 在我们的模型里，先色散、后在末端一次性加噪，和"边传边加噪"在统计上是等价的。

Ps = mean(abs(x).^2);
sigma2 = cfg.sig.sps * Ps / (cfg.sig.bitsPerSym * 10^(EbN0_dB/10));
n = sqrt(sigma2/2) * (randn(size(x)) + 1j*randn(size(x)));
y = x + n;
end
