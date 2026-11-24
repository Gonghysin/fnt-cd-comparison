function y = fnt_unpack2d(Z, q)
%FNT_UNPACK2D  2D 卷积结果折回成一维 N 点圆卷积的结果（M3）。
%   y = fnt_unpack2d(Z, q)
%   Z 是 2L x M 的 2D 卷积结果（行方向线性、列方向圆周）。
%
%   折回公式：  y(l,m) = Z(l,m) + Z(l+L, m-1)，   l = 0..L-1
%   推导：n = l + L*m，n-n1 中 l-l1<0 的项要向上一列借位，这些项落在 Z 的下半 L 行。
%   所以下半部分不能丢，折回要 L*M = N 次加法——这正是 Xing 里 INTT2D 比 NTT2D
%   多出的 N 次加法。折回后与一维 N 点圆卷积逐点完全相等。

L = size(Z, 1) / 2;
low = Z(1:L, :);
up  = Z(L+1:end, :);
y = mod(low + circshift(up, 1, 2), q);
y = y(:);
end
