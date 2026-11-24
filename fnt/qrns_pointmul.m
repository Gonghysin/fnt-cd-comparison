function [Yre, Yim] = qrns_pointmul(Xre, Xim, Hre, Him, b)
%QRNS_POINTMUL  用 QRNS 做 FNT 域里的复数逐点相乘（M3）。
%   [Yre, Yim] = qrns_pointmul(Xre, Xim, Hre, Him, b)
%   输入输出都是模 q=2^b+1 的整数（FNT 变换之后的实部、虚部）。
%
%   在模 q 下 J = 2^(b/2) 满足 J^2 = 2^b = -1，所以 J 扮演虚数单位 j。QRNS：
%       X+ = Xre + J*Xim,   X- = Xre - J*Xim        （正向转换）
%   复数乘法 (X)(H) 变成两路独立的实数逐点乘：  Y+ = X+ * H+,  Y- = X- * H-
%       Yre = (Y+ + Y-) / 2,   Yim = (Y+ - Y-) / (2J)  （反向转换）
%   1/2 与 1/J 是模 q 下的常数：1/2 = (q+1)/2，1/J = -J。
%   这样一次复数乘只要 2 次实数乘（直接算要 4 次，常见的 3 乘法也要 3 次）。

q = 2^b + 1;
J = 2^(b/2);
inv2 = (q + 1) / 2; % inv2 就是模 $q$ 下的"除以 2"。
invJ = q - J;

Xp = mod(Xre + J*Xim, q);   Xm = mod(Xre - J*Xim, q);
Hp = mod(Hre + J*Him, q);   Hm = mod(Hre - J*Him, q);

Yp = mod(Xp .* Hp, q);      Ym = mod(Xm .* Hm, q);

Yre = mod((Yp + Ym) * inv2, q);
Yim = mod(mod(Yp - Ym, q) * inv2, q);
Yim = mod(Yim * invJ, q);
end
