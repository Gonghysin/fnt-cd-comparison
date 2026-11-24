function tx = gen_tx(cfg)
%GEN_TX  发射端：随机比特 -> Gray QPSK -> RRC 成形，直接输出 2 SPS 信号（M1）。
%   tx = gen_tx(cfg)
%   调用前由脚本入口 rng(cfg.seed) 设定随机种子。
%
%   输出 tx 的字段：
%     bits   发送比特 (nBits x 1)
%     sym    QPSK 符号 (nSym x 1)，单位平均功率
%     h      RRC 滤波器（单位能量），接收端匹配滤波也用它
%     x      发射信号 (nSym*sps x 1)，已归一化为单位平均功率；
%            第 k 个符号的中心在第 (k-1)*sps+1 个样本
%
%   说明：
%   - 采样率 fs = 2*Rs，信号带宽约 Rs*(1+rolloff) 的一半（单边），不会混叠，
%     所以不需要先升采样再抽取，成形时直接在 2 SPS 上做。
%   - Gray 映射：比特 0 -> +1，比特 1 -> -1；奇数位走 I，偶数位走 Q。

nBits = cfg.sig.nBits;
sps   = cfg.sig.sps;
nSym  = cfg.sig.nSym;

% 把随机比特变成 QPSK 符号

bits = randi([0 1], nBits, 1);
b1 = bits(1:2:end);
b2 = bits(2:2:end);
sym = ((1 - 2*b1) + 1j*(1 - 2*b2)) / sqrt(2); % Gray映射 相邻的星座点只差 1 个比特
% sym长度是 nBits/2

h = rcosdesign(cfg.pulse.rolloff, cfg.pulse.span, sps, 'sqrt');   % 单位能量 RRC ：根升余弦滤波器
h = h(:);
delay = (numel(h) - 1) / 2;                  % numel(h)=抽头数， 滤波器群时延（样本）

x = upfirdn(sym, h, sps);                    % 升采样（插入0） + 成形
% 升采样后长度是 nBits (2e17)
x = x(delay + 1 : delay + nSym*sps);         % 去掉滤波器时延，保持 nSym*sps 个样本
x = x / sqrt(mean(abs(x).^2));               % 归一化到单位平均功率

tx.bits = bits;
tx.sym  = sym;
tx.h    = h;
tx.x    = x;
end
