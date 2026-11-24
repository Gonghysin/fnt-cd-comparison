function Xhat = fnt_pack2d(x, L, M)
%FNT_PACK2D  把长 N=L*M 的一维序列排成 2L x M 的二维阵列，下半部补零（M3）。
%   Xhat = fnt_pack2d(x, L, M)
%   下标对应  n = l + L*m （l 为行、m 为列，连续样本沿列向下排），
%   即 MATLAB 的列优先 reshape；下半 L 行补零，使行方向变成"线性卷积"。
%   与 Xing 图 3 一致：一半的行是 0，每列是"L 点序列后补 L 个零"。

x = x(:);
assert(numel(x) == L*M, 'x 的长度需等于 L*M');
Xhat = [reshape(x, L, M); zeros(L, M)];
end
