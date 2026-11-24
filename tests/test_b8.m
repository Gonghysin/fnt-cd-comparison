%% test_b8  b=8 改进探索的验收（fnt_cconv、hq_inband_opt、fnt_mp_design、eq_fdc_fnt_mp、cx_fdc_fnt_mp）
% 检查项：
%   1) fnt_cconv（1D 与 2D）在不溢出时与精确整数圆卷积逐点相等
%   2) hq_inband_opt 只会降低带内加权误差
%   3) 多遍方案在 80 km 跑真正的 FNT：无溢出、与精确卷积一致，且比单遍更准
%   4) cx_fdc_fnt_mp 在 1 遍时与 cx_fdc_fnt 相同，遍数越多越贵

addpath(fileparts(fileparts(mfilename('fullpath'))));  setup_paths;
fprintf('=== b=8 改进探索 验收 ===\n');

%% 1) fnt_cconv
rng(1);
for t = 1:2
    if t == 1, N = 256; L = 0; else, N = 2048; L = 16; end
    x = randi([-3 3], N, 4);  xi = randi([-3 3], N, 4);
    h = complex(zeros(N,1));  h(1:9) = complex(randi([-2 2], 9, 1), randi([-2 2], 9, 1));
    Y = fnt_cconv(x, xi, h, 8, L);
    E = round(ifft(fft(complex(x, xi)) .* fft(h)));
    assert(max(abs([real(E(:)); imag(E(:))])) <= 128 && all(Y(:) == E(:)), '检查 1 失败（N=%d）', N);
end
fprintf('[1] fnt_cconv 1D(N=256) 与 2D(N=2048, L=16) 均与精确卷积一致\n');

%% 2) 带内优化
N = 128;  tgt = 1.3 * exp(1j * 2*pi*rand(N,1)) .* (rand(N,1) < 0.5);
w = [ones(40,1); zeros(N-80,1); ones(40,1)] + 1e-3;
J = @(hh) sum(w .* abs(fft(hh - tgt)).^2);
h0 = round(tgt);  h1 = hq_inband_opt(h0, tgt, w);
assert(J(h1) <= J(h0) && all(h1 == round(h1)), '检查 2 失败');
fprintf('[2] 带内优化：加权误差 %.1f -> %.1f（降 %.1f dB）\n', J(h0), J(h1), 10*log10(J(h0)/J(h1)));

%% 3) 多遍方案（80 km）
cfg = cfg_fnt_cd();  rng(cfg.seed);  tx = gen_tx(cfg);
rng(cfg.seed + 999);  txc = gen_tx(cfg);
valid = (cfg.rx.edgeSym+1) : (cfg.sig.nSym - cfg.rx.edgeSym);
rs0 = rx_demod(tx.x, tx.h, cfg);
z = cfg.fiber.dist(1);  N = cfg.blk.N(1);
rng(cfg.seed + 901);  rcal = add_noise(apply_cd(txc.x, cfg, z), cfg, 8);
r0 = apply_cd(tx.x, cfg, z);
err = @(y) 20*log10(norm(sel(rx_demod(y, tx.h, cfg), valid) - rs0(valid)) / norm(rs0(valid)));
d1 = fnt_mp_design(rcal, cfg, z, N, struct('sx', 2.8));
d4 = fnt_mp_design(rcal, cfg, z, N, struct('sx', 22, 'Xs', 2, 'Hs', 2));
[y1, i1] = eq_fdc_fnt_mp(r0, cfg, z, N, d1);
[y4, i4] = eq_fdc_fnt_mp(r0, cfg, z, N, d4);
e1 = err(y1);  e4 = err(y4);
fprintf('[3] 80 km：1 遍 %.1f dB（%d 遍），输入 2 x 滤波器 2 共 %d 遍 %.1f dB；溢出 %d/%d，不一致 %d/%d\n', ...
    e1, i1.passes, i4.passes, e4, i1.nOverflow, i4.nOverflow, i1.nMismatch, i4.nMismatch);
assert(i1.nOverflow == 0 && i4.nOverflow == 0 && i1.nMismatch == 0 && i4.nMismatch == 0, '检查 3 失败：溢出或不一致');
assert(i4.passes == 4 && e4 < e1 - 6, '检查 3 失败：多遍没有明显更准');

%% 4) 复杂度
for N = [64 512 2048]
    assert(abs(cx_fdc_fnt_mp(N, 8, 1, 1) - cx_fdc_fnt(N, 8)) < 1e-6, '检查 4 失败：1 遍不等于原计数');
    assert(cx_fdc_fnt_mp(N, 8, 1, 2) > cx_fdc_fnt(N, 8) && cx_fdc_fnt_mp(N, 8, 2, 4) > cx_fdc_fnt_mp(N, 8, 1, 2), '检查 4 失败');
end
fprintf('[4] 多遍复杂度计数：1 遍与 cx_fdc_fnt 相同，遍数增加成本单调上升\n');
fprintf('test_b8 通过\n');

function v = sel(a, i), v = a(i); end
