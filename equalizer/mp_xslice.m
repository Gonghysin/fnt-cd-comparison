function [xs, wx] = mp_xslice(X, o)
%MP_XSLICE  输入量化 xq = round(X*sx)，再按基 S 切成 o.Xs 片：xq = sum wx(t)*xs{t}。
%   Xs=1：xs{1} = xq（限幅到 +-2^(b-1)）。Xs=2：xs{1} = round(xq/S)，xs{2} = xq - S*xs{1}（在 [-S/2, S/2]）。
lim = 2^o.b / 2;
c = lim * o.S^(o.Xs - 1);                                 % 限幅：保证最高片不超出 +-lim
xq = complex(max(min(round(real(X)*o.sx), c), -c), max(min(round(imag(X)*o.sx), c), -c));
xs = {};  wx = [];  rem = xq;
for t = 1:o.Xs
    Sw = o.S^(o.Xs - t);
    hi = complex(round(real(rem)/Sw), round(imag(rem)/Sw));
    hi = complex(max(min(real(hi), lim), -lim), max(min(imag(hi), lim), -lim));
    xs{end+1} = hi;  wx(end+1) = Sw;  %#ok<AGROW>
    rem = rem - Sw * hi;
end
end
