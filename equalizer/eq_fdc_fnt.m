function [y, info] = eq_fdc_fnt(r, cfg, z, N, opts)
%EQ_FDC_FNT  方案 C：FNT 变换域色散补偿（DFNT-TrDE），overlap-save 分块（M4）。
%   [y, info] = eq_fdc_fnt(r, cfg, z, N, opts)
%   r     2 SPS 接收序列；z 传输距离（m）；N 块长（FNT 点数）
%   opts  可选结构体（缺省值见下）：
%           b     16     位宽，q = 2^b+1
%           xbit  4      输入量化位宽（含符号位）
%           hbit  4      滤波器系数量化位宽（含符号位）
%           xmax  3      输入限幅：|实部|、|虚部| 超过 xmax 的截到 xmax（相当于 ADC 满量程）
%           L     0      0 表示 1D（要求 N <= 2^b）；>0 表示 2D，N = L*M，行 M 点、列 2L 点
%                        N > 2^b 时缺省自动取 L=16
%   y     补偿后序列（与 r 等长，已对齐）
%   info  诊断：sx, sh（两个量化比例）、maxAbs（有效区卷积结果的最大绝对值）、
%         limit（(q-1)/2，超过即溢出）、nOverflow（有效区溢出的数目）、
%         nMismatch（未溢出处 FNT 结果与精确整数卷积不一致的个数，应为 0）
%
%   流程（与 eq_fdc_fft 的分块完全相同，只是把"FFT·H·IFFT"换成"FNT·QRNS 乘·IFNT"）：
%     离线：理想逆色散在 N 点网格上的采样 -> N 点 IFFT（得到长 N、零时刻在下标 0 的时域系数）
%           -> 量化成整数（不截断抽头，与 eq_fdc_fft 用的是同一个滤波器）
%           -> FNT（实部、虚部各一次，2D 则打包成 2L x M）-> 得到 H 的 FNT 域系数
%     （这就是 Xing 的"IFFT -> 截位 -> FNT"预处理。若把抽头截到 P+1 个，误差会大很多。）
%     每块：量化成整数 -> mod q -> FNT（实部、虚部各一次）-> QRNS 逐点复数乘 -> IFNT
%           -> 提升成有符号数 -> 除以 sx*sh 还原 -> 保留中间 B 个
%   FNT 圆卷积在模 q 下是精确的；只要真实卷积结果的绝对值不超过 (q-1)/2 就不溢出，
%   所以误差只来自：输入量化、系数量化、系数截断，以及溢出（若发生）。

if nargin < 5, opts = struct(); end
o = fill_defaults(opts, N);
b = o.b;  q = 2^b + 1;

r = r(:);
Ls = numel(r);
P = round(N * cfg.blk.overlap);
B = N - P;
assert(mod(P, 2) == 0 && mod(Ls, B) == 0, '块参数与信号长度不匹配');
nBlk = Ls / B;
use2d = o.L > 0;

%% 离线：滤波器系数
h = ifft(cd_response(cfg, -z, N));                      % N 点，零时刻在下标 0，负时刻在末尾
sh = (2^(o.hbit-1) - 1) / max(abs([real(h); imag(h)]));
hN = round(h * sh);                                     % 整数（复数）

%% 分块
start = (0:nBlk-1) * B - P/2;
idx = mod(start + (0:N-1).', Ls) + 1;
R = r(idx);                                             % N x nBlk

sx = (2^(o.xbit-1) - 1) / o.xmax;
mx = 2^(o.xbit-1) - 1;
xr = round(max(min(real(R) * sx, mx), -mx));            % 限幅 + 量化
xi = round(max(min(imag(R) * sx, mx), -mx));

%% FNT 域处理（1D 或 2D，见 fnt_cconv）
Yint = fnt_cconv(xr, xi, hN, b, o.L * use2d);       % N x nBlk，有符号整数
Yb = Yint / (sx * sh);
y = reshape(Yb(P/2+1 : P/2+B, :), [], 1);

%% 诊断：与精确整数卷积对比（仅用于检查，不属于算法）
exact = round(ifft(fft(complex(xr, xi)) .* fft(hN)));   % 每列的圆卷积，真实整数值
lim = (q - 1) / 2; % 幅度限制
rows = P/2+1 : P/2+B; % 有效区
ev = exact(rows, :);
over = abs(real(ev)) > lim | abs(imag(ev)) > lim;
info.sx = sx;  info.sh = sh;  info.limit = lim;
info.maxAbs = max(abs([real(ev(:)); imag(ev(:))]));
info.nOverflow = nnz(over);
info.nMismatch = nnz(~over & (Yint(rows, :) ~= ev));
end

function o = fill_defaults(opts, N) % 用默认值填充缺省
o.b = 16;  o.xbit = 4;  o.hbit = 4;  o.xmax = 3;  o.L = 0;
f = fieldnames(opts);
for k = 1:numel(f), o.(f{k}) = opts.(f{k}); end
if o.L == 0 && N > 2^o.b, o.L = 16; end
end
