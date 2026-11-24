function Y = fnt_1d(X, par, dir)
%FNT_1D  沿列方向做一维 FNT / IFNT，radix-2 DIT 蝶形，精确整数运算（M3）。
%   Y = fnt_1d(X, par, dir)
%   X    N x M 整数矩阵，元素在 [0, q-1]（负数先 mod q）；每一列各做一次变换
%   par  fnt_params(b, N) 的输出
%   dir  +1 正变换，-1 逆变换（用 omega^-1，并乘 N^-1）
%
%   定义：  X[k] = sum_n x[n] * omega^(n*k)  mod q
%   和 FFT 一样的蝶形：  a' = a + w*b,  b' = a - w*b  （全部 mod q）。
%   这里的 w*b 直接用 mod 乘法实现；当 w 是 2 的幂时它在硬件里是移位
%   （见 fnt_mulpow2），乘法复杂度的计数留到 M7。

N = par.N;  q = par.q;
assert(size(X, 1) == N, '输入行数需等于 N');
if dir > 0, W = par.W; else, W = par.Winv; end

Y = mod(X(par.rev + 1, :), q);                 % 位反转重排
M = size(Y, 2);

% 蝶形运算 算FNT 省运算
for s = 1:par.logN
    m = 2^s;  half = m/2;
    tw = W(1 : N/m : N/2);                     % omega^(k*N/m), k = 0..half-1
    Yr = reshape(Y, m, N/m, M);                % 每个 m 行为一个蝶形组
    a  = Yr(1:half, :, :);
    bb = mod(Yr(half+1:m, :, :) .* tw, q);
    Yr(1:half, :, :)   = mod(a + bb, q);
    Yr(half+1:m, :, :) = mod(a - bb, q);
    Y = reshape(Yr, N, M);
end
if dir < 0
    Y = mod(Y * par.Ninv, q);
end
end
