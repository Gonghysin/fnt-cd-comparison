function [y, info] = eq_fdc_fft_fixed(r, cfg, z, N, opts)
%EQ_FDC_FFT_FIXED  方案 B 的位精确（定点）版：FFT 频域色散补偿，overlap-save 分块（M5）。
%   [y, info] = eq_fdc_fft_fixed(r, cfg, z, N, opts)
%   r 2 SPS 接收序列；z 传输距离（m）；N 块长（FFT 点数）。分块与 eq_fdc_fft 完全相同。
%   opts 可选结构体（缺省值见下）：
%       W       16       数据通路字长（有符号，含符号位）。输入、每层蝶形输出都是 W 位整数
%       xmax    3        输入限幅（相当于 ADC 满量程），与 eq_fdc_fnt 一致
%       rnd     false    false：移位时直接截断（丢低位，向下取整）；true：四舍五入
%       scheme  'fwd'    缩放方案：
%                 'fwd'  正向 FFT 每层右移 1 位（总共除以 N），逆向 IFFT 不缩放。
%                        数学上总增益为 1，逆向各层的值单调增长到输出，不会溢出。
%                 'both' 正向、逆向每层都右移 1 位（Xing 第 VI 节的字面说法："每层截 1 位"）。
%                        这样总增益是 1/N，最后要左移 log2(N) 位恢复，丢掉 log2(N) 位精度。
%
%   定点格式（所有数都是整数，用 double 存放，乘积 < 2^53，所以精确）：
%       输入       x_int = round(x * sx)，sx = (2^(W-1)-1)/xmax
%       旋转因子   cos、sin 乘 2^F 取整，F = W-2 位小数（1.0 对应 2^F）
%       滤波器 H   同样乘 2^F 取整（复数，|H|=1）
%       复数乘     乘积右移 F 位回到原来的比例（旋转因子乘、逐点乘 H 都是）
%   每个蝶形：a ± (b*w >> F)，再按方案右移 1 位。写回前检查是否超出 W 位，
%   超出则饱和并计数（info.nSat）。
%
%   info：sx、nSat（饱和次数，应为 0）、maxAbs（所有中间值的最大绝对值）、limit（2^(W-1)-1）。

if nargin < 5, opts = struct(); end
o.W = 16;  o.xmax = 3;  o.rnd = false;  o.scheme = 'fwd';
f = fieldnames(opts);
for k = 1:numel(f), o.(f{k}) = opts.(f{k}); end
W = o.W;  F = W - 2;  lim = 2^(W-1) - 1;
assert(2*W - 3 <= 52, 'W 太大，乘积超出 double 的精确整数范围');

r = r(:);
Ls = numel(r);
P = round(N * cfg.blk.overlap);
B = N - P;
assert(mod(P, 2) == 0 && mod(Ls, B) == 0, '块参数与信号长度不匹配');
nBlk = Ls / B;
logN = log2(N);
assert(N == 2^logN, 'N 需为 2 的幂');

%% 常数：旋转因子、滤波器（离线，量化成 F 位小数）
k = (0:N/2-1).';
twc = round(cos(2*pi*k/N) * 2^F);
tws = round(sin(2*pi*k/N) * 2^F);
Hc = cd_response(cfg, -z, N);
Hr = round(real(Hc) * 2^F);
Hi = round(imag(Hc) * 2^F);

%% 分块 + 输入量化
start = (0:nBlk-1) * B - P/2;
idx = mod(start + (0:N-1).', Ls) + 1;
R = r(idx);
sx = lim / o.xmax;
xr = round(max(min(real(R) * sx, lim), -lim));
xi = round(max(min(imag(R) * sx, lim), -lim));

rev = bitrev(N);
stat.nSat = 0;  stat.maxAbs = 0;

%% 正向 FFT -> 逐点乘 H -> 逆向 IFFT
scaleFwd = true;
scaleInv = strcmp(o.scheme, 'both');
[Xr, Xi, stat] = fft_fx(xr, xi, N, +1, scaleFwd, twc, tws, F, lim, o.rnd, rev, stat);
Zr = shr(Xr .* Hr - Xi .* Hi, F, o.rnd);
Zi = shr(Xr .* Hi + Xi .* Hr, F, o.rnd);
[stat, Zr, Zi] = sat(stat, Zr, Zi, lim);
[Yr, Yi, stat] = fft_fx(Zr, Zi, N, -1, scaleInv, twc, tws, F, lim, o.rnd, rev, stat);
if scaleInv                                              % 'both'：左移 log2(N) 位恢复比例
    Yr = Yr * N;  Yi = Yi * N;
    [stat, Yr, Yi] = sat(stat, Yr, Yi, lim);
end

Yb = complex(Yr, Yi) / sx;
y = reshape(Yb(P/2+1 : P/2+B, :), [], 1);

info.sx = sx;  info.nSat = stat.nSat;  info.maxAbs = stat.maxAbs;  info.limit = lim;
end

function [Xr, Xi, stat] = fft_fx(Xr, Xi, N, dir, scale, twc, tws, F, lim, rnd, rev, stat)
%FFT_FX  定点 radix-2 DIT FFT（dir=+1 正向，-1 逆向），沿列方向。
Xr = Xr(rev + 1, :);  Xi = Xi(rev + 1, :);
M = size(Xr, 2);
for s = 1:log2(N)
    m = 2^s;  half = m/2;
    wr = twc(1 : N/m : N/2);
    wi = -dir * tws(1 : N/m : N/2);                     % 正向 e^{-j...}，逆向取共轭
    Ar = reshape(Xr, m, N/m, M);  Ai = reshape(Xi, m, N/m, M);
    ar = Ar(1:half, :, :);  ai = Ai(1:half, :, :);
    br = Ar(half+1:m, :, :);  bi = Ai(half+1:m, :, :);
    tr = shr(br .* wr - bi .* wi, F, rnd);               % b*w，右移 F 位回到原比例
    ti = shr(br .* wi + bi .* wr, F, rnd);
    u_r = ar + tr;  u_i = ai + ti;                       % a + b*w
    l_r = ar - tr;  l_i = ai - ti;                       % a - b*w
    if scale
        u_r = shr(u_r, 1, rnd);  u_i = shr(u_i, 1, rnd);
        l_r = shr(l_r, 1, rnd);  l_i = shr(l_i, 1, rnd);
    end
    [stat, u_r, u_i] = sat(stat, u_r, u_i, lim);
    [stat, l_r, l_i] = sat(stat, l_r, l_i, lim);
    Ar(1:half, :, :) = u_r;  Ar(half+1:m, :, :) = l_r;
    Ai(1:half, :, :) = u_i;  Ai(half+1:m, :, :) = l_i;
    Xr = reshape(Ar, N, M);  Xi = reshape(Ai, N, M);
end
end

function v = shr(v, s, rnd)
%SHR  算术右移 s 位：截断（向下取整）或四舍五入。
if s == 0, return; end
if rnd
    v = floor((v + 2^(s-1)) / 2^s);
else
    v = floor(v / 2^s);
end
end

function [stat, a, b] = sat(stat, a, b, lim)
%SAT  检查是否超出 W 位，超出则饱和并计数；同时记录最大绝对值。
m = max(max(abs(a(:))), max(abs(b(:))));
stat.maxAbs = max(stat.maxAbs, m);
over = abs(a) > lim | abs(b) > lim;
if any(over(:))
    stat.nSat = stat.nSat + nnz(over);
    a = max(min(a, lim), -lim);
    b = max(min(b, lim), -lim);
end
end

function rev = bitrev(N)
nb = log2(N);  rev = zeros(N, 1);
for i = 0:N-1
    v = i;  rr = 0;
    for t = 1:nb
        rr = bitor(bitshift(rr, 1), bitand(v, 1));
        v = bitshift(v, -1);
    end
    rev(i+1) = rr;
end
end
