function [rsym, hatBits] = rx_demod(r, h, cfg)
%RX_DEMOD  接收机：匹配滤波 -> 下采样到符号率 -> QPSK 硬判决（M1）。
%   [rsym, hatBits] = rx_demod(r, h, cfg)
%
%   输入只有三样，都是接收机真正能知道的：
%     r    （色散补偿后的）2 SPS 接收序列
%     h    RRC 滤波器系数（收发双方事先约定的参数，不是发送的数据）
%     cfg  系统参数（采样率、sps 等）
%   这里**看不到**发送的比特；评分在 score_ber 里做。
%
%   输出：rsym 符号率星座点；hatBits 判决出的比特 (nBits x 1)。
%   不做载波恢复和定时恢复（本阶段假设理想，采样相位固定取第 1 个样本起）。

sps = cfg.sig.sps;

y = conv(r(:), h(:), 'same');          % 匹配滤波；h 对称且长度为奇数，'same' 自动补偿群时延
rsym = y(1:sps:end);                   % 下采样：取符号中心样本

hatBits = zeros(cfg.sig.nBits, 1);
hatBits(1:2:end) = real(rsym) < 0;     % 实部为负 -> 比特 1
hatBits(2:2:end) = imag(rsym) < 0;     % 虚部为负 -> 比特 1
end
