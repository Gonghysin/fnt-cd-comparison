function [FA, info] = cx_fdc_fnt(N, b)
%CX_FDC_FNT  方案 C（FNT-TrDE）一个 overlap-save 块的复杂度，单位 1 位全加器（M7）。
%   [FA, info] = cx_fdc_fnt(N, b)
%   N <= 2^b 用 1D-HR；N > 2^b 用 2D-HR（L 在可行范围内取成本最小）。
%
%   一个块 = 实部、虚部各一次正向 FNT（2 次）+ QRNS + 逐点乘 + 2 次逆向 FNT。
%     QRNS 正向 X+ = Xr + J*Xi, X- = Xr - J*Xi：每点 2 次加（J 是 2 的幂，移位免费）
%     QRNS 逆向同样每点 2 次加；逆变换的 1/N 是移位，免费
%     逐点乘：滤波器预先算好，X+、X- 各乘一个常数 -> 每点 2 次 b 位常数乘 = 2*b/3 次加
%   1D：4 次变换 + N 点逐点乘 + QRNS 共 4N 次加。
%   2D(2L x M)：行 M 点变换只对 L 个非零行做（另 L 行全零）；列 2L 点变换对 M 列做；
%     逆向 行、列都满做；逐点乘 2N 点；QRNS 8N 次加；折叠 y(l,m)=Z(l,m)+Z(l+L,m-1) 每分量 N 次加。
%     没有计入 Xing 第 V-C 节说的"节省一层加法"的优化，所以 2D 的值是偏保守的上界。
%   b 位 D1 加法 = b 个全加器（Xing 第 IV 节）。
if N <= 2^b
    [T, R] = cx_fnt1d(N, b);
    adds = 4*T + N*(2*b/3) + 4*N;
    info = struct('mode', '1D-HR', 'R', R, 'L', 0, 'adds_b', adds);
else
    best = inf;
    for L = 2.^(0:log2(N)-1)
        M = N / L;
        if M < 2 || 2*L > 2^b || M > 2^b, continue; end
        [Tr, Rr] = cx_fnt1d(M, b);
        [Tc, Rc] = cx_fnt1d(2*L, b);
        % 正向：L 个非零行 + M 列；逆向：2L 行 + M 列；实/虚（或 +/-）各一份
        adds = 2*(L*Tr + M*Tc) + 2*(M*Tc + 2*L*Tr) + 2*N*(2*b/3) + 8*N + 2*N;
        if adds < best
            best = adds;  info = struct('mode', '2D-HR', 'R', [Rr Rc], 'L', L, 'adds_b', adds);
        end
    end
    adds = best;
end
FA = adds * b;
end
