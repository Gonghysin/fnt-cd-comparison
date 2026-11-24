function [FA, info] = cx_fdc_fnt_mp(N, b, Xs, nPass)
%CX_FDC_FNT_MP  多遍（multi-pass）FNT-TrDE 一个 overlap-save 块的复杂度，单位 1 位全加器。
%   [FA, info] = cx_fdc_fnt_mp(N, b, Xs, nPass)
%   Xs     输入切成几片（每片一次正向变换，实部、虚部各一次）
%   nPass  总遍数 = 输入片数 x 滤波器片数（每遍一次逐点乘 + 一次逆变换，实部、虚部各一次）
%   Xs = nPass = 1 时与 cx_fdc_fnt 完全相同。
%
%   计数（规则同 cx_fdc_fnt）：
%     正向变换  2*Xs 次；QRNS 正向 2 次加/点 x N 点 x Xs
%     每一遍    逐点乘 2*b/3 次加/点；QRNS 逆向 2 次加/点；逆变换 2 次（2D 时含折叠 N 次加 x 2）
%     重组      第 2 遍起，每遍的有效输出（实部、虚部各 B=N/2 个）要乘一个常数比例再加到累加器：
%               按 16 位常数乘（16/3 次 16 位加）+ 1 次 16 位加计，即每个样本 16*(16/3+1) 个全加器。
%               这一项不在 b 位域里，单独以全加器计。
%   N <= 2^b 用 1D-HR；否则 2D-HR（L 取单遍成本最小的值，与 cx_fdc_fnt 一致）。

[~, base] = cx_fdc_fnt(N, b);
if strcmp(base.mode, '1D-HR')
    T = cx_fnt1d(N, b);
    fwdT = 2*T;  invT = 2*T;  pw = N;          % 每次（实+虚）正向 / 逆向变换成本；逐点乘点数
    fold = 0;
else
    L = base.L;  M = N / L;
    Tr = cx_fnt1d(M, b);  Tc = cx_fnt1d(2*L, b);
    fwdT = 2*(L*Tr + M*Tc);                    % 只有 L 个非零行
    invT = 2*(M*Tc + 2*L*Tr);
    pw = 2*N;                                  % 2L x M 个点
    fold = 2*N;
end
adds = Xs*(fwdT + 2*pw) + nPass*(invT + fold + pw*(2*b/3) + 2*pw);
FA = adds * b + (nPass - 1) * N * 16 * (16/3 + 1);
info = struct('mode', base.mode, 'adds_b', adds, 'Xs', Xs, 'nPass', nPass);
end
