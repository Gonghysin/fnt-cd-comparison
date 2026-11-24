function [y, info] = eq_fdc_fnt_mp(r, cfg, z, N, d)
%EQ_FDC_FNT_MP  多遍 FNT-TrDE：每一遍是一次真正的 FNT 圆卷积（fnt_cconv），结果在普通整数里重组。
%   [y, info] = eq_fdc_fnt_mp(r, cfg, z, N, d)
%   d  fnt_mp_design 的输出（离线设计：各遍整数滤波器 d.hN{p}、比例 d.a(p)、所用输入片 d.t(p)）
%   分块与 eq_fdc_fnt / eq_fdc_fft 完全相同（居中 overlap-save，重叠 50%）。
%   N > 2^b 时用 2D-FNT（L = 16），与 eq_fdc_fnt 一致。
%
%   info：passes（遍数）、maxAbs（各遍有效区卷积结果的最大绝对值）、limit、
%         nOverflow（有效区溢出个数，应为 0）、nMismatch（未溢出处 FNT 与精确整数卷积不一致的个数，应为 0）

o = d.o;
L = 0;  if N > 2^o.b, L = 16; end
X = mp_blocks(r, N);
[xs, wx] = mp_xslice(X, o);
P = round(N * cfg.blk.overlap);  B = N - P;  rows = P/2+1 : P/2+B;
lim = 2^o.b / 2;
acc = zeros(N, size(X, 2));
info.passes = numel(d.hN);  info.limit = lim;  info.maxAbs = 0;  info.nOverflow = 0;  info.nMismatch = 0;
for p = 1:numel(d.hN)
    xt = xs{d.t(p)};
    Yint = fnt_cconv(real(xt), imag(xt), d.hN{p}, o.b, L);
    acc = acc + wx(d.t(p)) * Yint / d.a(p);
    % 诊断（不属于算法）：与精确整数卷积对比
    ev = round(ifft(fft(xt) .* fft(d.hN{p})));  ev = ev(rows, :);
    over = abs(real(ev)) > lim | abs(imag(ev)) > lim;
    info.maxAbs = max(info.maxAbs, max(abs([real(ev(:)); imag(ev(:))])));
    info.nOverflow = info.nOverflow + nnz(over);
    info.nMismatch = info.nMismatch + nnz(~over & (Yint(rows, :) ~= ev));
end
y = reshape(acc(rows, :) / o.sx, [], 1);
end
