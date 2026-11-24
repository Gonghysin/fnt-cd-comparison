function Y = fnt_cconv(xr, xi, hN, b, L)
%FNT_CCONV  用 FNT 做一批复数圆卷积（每列一块），结果提升为有符号整数。
%   Y = fnt_cconv(xr, xi, hN, b, L)
%   xr, xi  N x nBlk 整数（输入实部、虚部）；hN  N x 1 复整数滤波器；q = 2^b+1
%   L       0：1D-FNT（要求 N <= 2^b）；>0：2D-FNT，N = L*M 排成 2L x M（行 M 点、列 2L 点）
%   Y       N x nBlk 复整数 = 圆卷积 x (*) hN；真实值超出 +-(q-1)/2 时会回绕（溢出），调用方负责检查。
%   流程：实部、虚部各一次 FNT -> QRNS 逐点复数乘 -> 逆 FNT -> 提升到 +-(q-1)/2（M3/M4 的组合）。
%   从 eq_fdc_fnt 抽出来，供单遍（eq_fdc_fnt）和多遍（eq_fdc_fnt_mp）共用。

q = 2^b + 1;
N = size(xr, 1);
if L == 0
    par = fnt_params(b, N);
    Xr = fnt_1d(mod(xr, q), par, 1);   Xi = fnt_1d(mod(xi, q), par, 1);
    Hr = fnt_1d(mod(real(hN), q), par, 1);  Hi = fnt_1d(mod(imag(hN), q), par, 1);
    [Yr, Yi] = qrns_pointmul(Xr, Xi, Hr, Hi, b);
    yr = fnt_lift(fnt_1d(Yr, par, -1), q);
    yi = fnt_lift(fnt_1d(Yi, par, -1), q);
else
    M = N / L;
    assert(2*L <= 2^b && M <= 2^b, '2D 的行列长度需不超过 2^b');
    parR = fnt_params(b, M);  parC = fnt_params(b, 2*L);
    pack = @(v) cat(1, reshape(v, L, M, []), zeros(L, M, size(v, 2)));   % -> 2L x M x nBlk
    Xr = fnt2d_batch(pack(mod(xr, q)), parR, parC, 1);
    Xi = fnt2d_batch(pack(mod(xi, q)), parR, parC, 1);
    Hr = fnt2d_batch(pack(mod(real(hN), q)), parR, parC, 1);
    Hi = fnt2d_batch(pack(mod(imag(hN), q)), parR, parC, 1);
    [Yr, Yi] = qrns_pointmul(Xr, Xi, Hr, Hi, b);
    zr = fnt2d_batch(Yr, parR, parC, -1);
    zi = fnt2d_batch(Yi, parR, parC, -1);
    yr = fnt_lift(fold(zr, L, M, q), q);
    yi = fnt_lift(fold(zi, L, M, q), q);
end
Y = complex(yr, yi);
end

function Y = fnt2d_batch(X, parR, parC, dir)
%FNT2D_BATCH  对 2L x M x nBlk 的每一页做 2D-FNT / IFNT（与 fnt_2d 同一算法，批量版）。
[L2, M, nb] = size(X);
if dir > 0
    A = fnt_1d(reshape(permute(X, [2 1 3]), M, L2*nb), parR, 1);     % 先行（M 点）
    A = reshape(permute(reshape(A, M, L2, nb), [2 1 3]), L2, M*nb);
    Y = reshape(fnt_1d(A, parC, 1), L2, M, nb);                       % 后列（2L 点）
else
    A = fnt_1d(reshape(X, L2, M*nb), parC, -1);
    A = reshape(permute(reshape(A, L2, M, nb), [2 1 3]), M, L2*nb);
    A = fnt_1d(A, parR, -1);
    Y = permute(reshape(A, M, L2, nb), [2 1 3]);
end
end

function y = fold(Z, L, M, q)
%FOLD  2D 卷积结果折回一维：y(l,m) = Z(l,m) + Z(l+L, m-1)（见 fnt_unpack2d）。
low = Z(1:L, :, :);
up  = Z(L+1:2*L, :, :);
y = mod(low + circshift(up, 1, 2), q);
y = reshape(y, L*M, []);
end
