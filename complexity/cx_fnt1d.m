function [T, Ropt] = cx_fnt1d(N, b, R)
%CX_FNT1D  N 点 1D-FNT 的成本，单位 b 位加法（M7）。
%   T = cx_fnt1d(N, b, R)    固定基 R（2 的幂）；R=2 即 1D-R2
%   [T, Ropt] = cx_fnt1d(N, b)  在 R = 2..min(N,128) 里取最小的（1D-HR）
%
%   算法：DIT。N = R*(N/R)：R 个 N/R 点子变换；级间旋转因子 w^(j*k)，j=1..R-1，
%   k=0..N/R-1，按 cx_twcost 分类计价；N/R 个 R 点变换，每个 R 点用 R2 递归计价
%   （R 小时它的核 w^(N/R) 是 2 的幂，乘法全免，只剩 R*log2R 次加）。
%   每个 R2 蝶形 a +- b*w 计 2 次加 + 乘 w 的成本。
%   不含模约简（Xing 表 IV 注：NTT 的模约简不计）。
if nargin < 3
    best = inf;
    for r = 2.^(1:min(log2(N), 7))
        t = cx_fnt1d(N, b, r);
        if t < best, best = t;  Ropt = r; end
    end
    T = best;  return
end
par = fnt_params(b, N);
T = cost(N, par.omega, R, b);
Ropt = R;
end

function T = cost(N, w, R, b)
q = 2^b + 1;
if N == 1, T = 0; return; end
R = min(R, N);
M = N / R;
wR = mpow(w, R, q);
T = R * cost(M, wR, R, b);
if M > 1 || true
    [j, k] = ndgrid(1:R-1, 0:M-1);
    v = mpowv(w, j .* k, q);
    T = T + sum(cx_twcost(v(:), b));
end
wd = mpow(w, M, q);
T = T + M * cost2(R, wd, b);
end

function T = cost2(N, w, b)
% R2 递归（N 小，w 是 N 次本原根）
q = 2^b + 1;
if N == 1, T = 0; return; end
v = mpowv(w, (0:N/2-1).', q);
T = 2 * cost2(N/2, mpow(w, 2, q), b) + sum(cx_twcost(v, b)) + N;
end

function r = mpow(w, e, q)
r = mpowv(w, e, q);
end

function r = mpowv(w, e, q)
r = ones(size(e));  base = w;  ee = e;
while any(ee(:) > 0)
    odd = mod(ee, 2) == 1;
    r(odd) = mod(r(odd) * base, q);
    base = mod(base * base, q);
    ee = floor(ee / 2);
end
end
