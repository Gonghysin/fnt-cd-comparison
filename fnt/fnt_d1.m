function varargout = fnt_d1(mode, a, z)
%FNT_D1  Diminished-1（D1）表示的编码与解码（M3）。
%   [d, z] = fnt_d1('enc', x)   x 在 [0, 2^b] 内（共 2^b+1 种取值，b 位装不下）
%   x      = fnt_d1('dec', d, z)
%
%   D1：非零数存 x-1（取值 0..2^b-1，刚好 b 位），0 另用一个标志位 z 表示。
%   这样 F_n 的所有残余类只需 b 位 + 1 个标志位。本阶段只验证编码解码互逆；
%   位精确的 D1 加法留到定点仿真（M5 之后）。

switch mode
    case 'enc'
        zf = (a == 0);
        d = a - 1;
        d(zf) = 0;
        varargout = {d, zf};
    case 'dec'
        x = a + 1;
        x(z) = 0;
        varargout = {x};
    otherwise
        error('mode 需为 enc 或 dec');
end
end
