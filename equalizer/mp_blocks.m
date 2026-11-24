function X = mp_blocks(r, N)
%MP_BLOCKS  居中 overlap-save 分块（重叠 50%），与 eq_fdc_fft 相同：返回 N x nBlk。
r = r(:);  Ls = numel(r);  P = N/2;  B = N - P;  nb = Ls / B;
assert(mod(Ls, B) == 0, '块参数与信号长度不匹配');
start = (0:nb-1) * B - P/2;
X = r(mod(start + (0:N-1).', Ls) + 1);
end
