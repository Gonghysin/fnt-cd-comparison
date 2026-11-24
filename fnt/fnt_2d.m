function Y = fnt_2d(X, parR, parC, dir)
%FNT_2D  二维 FNT / IFNT：先对每行做 M 点，再对每列做 2L 点（M3）。
%   Y = fnt_2d(X, parR, parC, dir)
%   X     2L x M 整数矩阵（元素在 [0, q-1]）
%   parR  行变换（长度 M）的参数；parC 列变换（长度 2L）的参数，均由 fnt_params 给出
%   dir   +1 正变换，-1 逆变换（顺序相反）
%
%   二维变换是可分离的，行列之间没有旋转因子，所以只需要 M 点和 2L 点的变换核。
%   这样 1D 长度受限于 2^b 的问题就被绕开：N = L*M 可以到 2^(2b-1) 量级
%   （Xing 第 V-C 节："最大圆卷积长度由 Nmax 扩展到 Nmax^2/2"）。

if dir > 0
    Y = fnt_1d(X.', parR, 1).';
    Y = fnt_1d(Y, parC, 1);
else
    Y = fnt_1d(X, parC, -1);
    Y = fnt_1d(Y.', parR, -1).';
end
end
