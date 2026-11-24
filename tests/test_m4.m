%% test_m4  M4 验收：FNT 变换域色散补偿（方案 C）
% 检查项：
%   1) 正确性：FNT 结果与"精确整数卷积"逐点一致（nMismatch=0）；2D 与 1D 结果完全相同
%   2) b=16 精度：符号误差随位宽下降；(6,6) 位无溢出，误差 < -25 dB
%   3) b=16 的 BER（Eb/N0=6 dB，同一份噪声）与理想补偿、浮点 FFT-FDE 基本一致
%   4) b=8：在不溢出的前提下搜索最好的位宽组合，记录误差和 BER（不做通过条件，如实记录）
% 误差指标：补偿后符号率星座点与"无色散无噪声"参考的相对均方根误差（dB）。

addpath(fileparts(fileparts(mfilename('fullpath'))));  setup_paths;   % 自动加入项目路径

cfg = cfg_fnt_cd();
rng(cfg.seed);
tx = gen_tx(cfg);
nd = numel(cfg.fiber.dist);
valid = (cfg.rx.edgeSym+1) : (cfg.sig.nSym-cfg.rx.edgeSym);
rsym0 = rx_demod(tx.x, tx.h, cfg);
tic;

%% 1) 正确性
fprintf('=== M4 FNT 变换域色散补偿 ===\n\n[1] 正确性（无噪声）\n');
z = cfg.fiber.dist(4);  r = apply_cd(tx.x, cfg, z);          % 1120 km
[y1, i1] = eq_fdc_fnt(r, cfg, z, 512, struct('b', 16, 'xbit', 6, 'hbit', 6));          % 1D
[y2, i2] = eq_fdc_fnt(r, cfg, z, 512, struct('b', 16, 'xbit', 6, 'hbit', 6, 'L', 16)); % 2D, 32 x 32
fprintf('  b=16 N=512：1D 与 2D 的输出 %s；nMismatch 1D=%d 2D=%d\n', ...
    ternary(isequal(y1, y2), '完全相同', '不同'), i1.nMismatch, i2.nMismatch);
assert(isequal(y1, y2) && i1.nMismatch == 0 && i2.nMismatch == 0, 'M4 检查 1 失败');
[~, i3] = eq_fdc_fnt(r, cfg, z, 512, struct('b', 8, 'xbit', 3, 'hbit', 2));            % b=8 自动用 2D
fprintf('  b=8  N=512（自动 2D，L=16）：nMismatch=%d nOverflow=%d\n', i3.nMismatch, i3.nOverflow);
assert(i3.nMismatch == 0 && i3.nOverflow == 0, 'M4 检查 1 失败（b=8）');

%% 2) b=16 精度
cfgBits = [4 4; 5 5; 6 6];
errF = zeros(nd, 1);  errN = zeros(nd, size(cfgBits, 1));  maxAbs = zeros(nd, size(cfgBits, 1));
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    r = apply_cd(tx.x, cfg, z);
    errF(i) = symErr_dB(eq_fdc_fft(r, cfg, z, N), tx, cfg, rsym0, valid);
    for k = 1:size(cfgBits, 1)
        [y, inf_] = eq_fdc_fnt(r, cfg, z, N, struct('b', 16, 'xbit', cfgBits(k,1), 'hbit', cfgBits(k,2)));
        errN(i, k) = symErr_dB(y, tx, cfg, rsym0, valid);
        maxAbs(i, k) = inf_.maxAbs;
        assert(inf_.nMismatch == 0 && inf_.nOverflow == 0, 'M4 检查 2 失败：溢出或不一致 i=%d k=%d', i, k);
    end
end
fprintf('\n[2] b=16 符号误差(dB)，列为 (xbit,hbit)；最后一列为浮点 FFT-FDE；括号内为最大卷积值/限值\n');
hdr = arrayfun(@(a,b) sprintf('(%d,%d)', a, b), cfgBits(:,1), cfgBits(:,2), 'UniformOutput', false);
fprintf('%8s %6s', 'dist/km', 'N');  fprintf('%12s', hdr{:});  fprintf('%12s\n', 'float FFT');
for i = 1:nd
    fprintf('%8d %6d', cfg.fiber.dist_km(i), cfg.blk.N(i));  fprintf('%12.1f', errN(i,:));  fprintf('%12.1f\n', errF(i));
end
fprintf('  最大卷积值占限值(32768)的比例（6,6 位）：'); fprintf('%.2f ', maxAbs(:, 3)/32768); fprintf('\n');
assert(all(errN(:, 3) < -25), 'M4 检查 2 失败：b=16 (6,6) 误差过大');

%% 3) b=16 的 BER
EbN0 = 6;  rng(cfg.seed + 1);
fprintf('\n[3] Eb/N0 = %g dB，同一份噪声；比特错误数\n%8s %12s %12s %12s %12s\n', EbN0, 'dist/km', '理想逆色散', 'float FFT', 'FNT b16(6,6)', 'FNT b16(4,4)');
nErrI = zeros(nd, 1);  nErrF = nErrI;  nErr66 = nErrI;  nErr44 = nErrI;  rN = cell(nd, 1);
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    r = add_noise(apply_cd(tx.x, cfg, z), cfg, EbN0);
    rN{i} = r;
    nErrI(i)  = scoreOf(apply_cd(r, cfg, -z), tx, cfg).nErr;
    nErrF(i)  = scoreOf(eq_fdc_fft(r, cfg, z, N), tx, cfg).nErr;
    nErr66(i) = scoreOf(eq_fdc_fnt(r, cfg, z, N, struct('b', 16, 'xbit', 6, 'hbit', 6)), tx, cfg).nErr;
    nErr44(i) = scoreOf(eq_fdc_fnt(r, cfg, z, N, struct('b', 16, 'xbit', 4, 'hbit', 4)), tx, cfg).nErr;
    fprintf('%8d %12d %12d %12d %12d\n', cfg.fiber.dist_km(i), nErrI(i), nErrF(i), nErr66(i), nErr44(i));
    assert(abs(nErr66(i) - nErrI(i)) <= max(20, 0.12*nErrI(i)), 'M4 检查 3 失败：b=16 (6,6) 与理想补偿误码数差异过大');
end

%% 4) b=8：不溢出前提下搜索最好的位宽
xs = [3 4 5];  hs = [2 3 4];
fprintf('\n[4] b=8（q=257，限值 ±128）：每个距离在 xbit=%s、hbit=%s 里选"最大卷积值不超过限值的 90%%、误差最小"的组合\n', mat2str(xs), mat2str(hs));
fprintf('%8s %6s %10s %10s %12s %12s\n', 'dist/km', 'N', '最优(x,h)', '误差(dB)', '误码数', '理想误码数');
best = zeros(nd, 2);  err8 = zeros(nd, 1);  nErr8 = zeros(nd, 1);
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    r0 = apply_cd(tx.x, cfg, z);
    e = inf;  bx = 0;  bh = 0;
    for xb = xs
        for hb = hs
            [y, inf_] = eq_fdc_fnt(r0, cfg, z, N, struct('b', 8, 'xbit', xb, 'hbit', hb));
            if inf_.nOverflow == 0 && inf_.maxAbs <= 0.9*inf_.limit
                ee = symErr_dB(y, tx, cfg, rsym0, valid);
                if ee < e, e = ee;  bx = xb;  bh = hb; end
            end
        end
    end
    best(i, :) = [bx bh];  err8(i) = e;
    if bx > 0
        [y, inf_] = eq_fdc_fnt(rN{i}, cfg, z, N, struct('b', 8, 'xbit', bx, 'hbit', bh));
        nErr8(i) = scoreOf(y, tx, cfg).nErr;
        assert(inf_.nMismatch == 0, 'M4 检查 4 失败：b=8 与精确整数卷积不一致');
    else
        nErr8(i) = NaN;
    end
    fprintf('%8d %6d %10s %10.1f %12d %12d\n', cfg.fiber.dist_km(i), N, sprintf('(%d,%d)', bx, bh), e, nErr8(i), nErrI(i));
end

fprintf('\n用时 %.0f 秒\nM4 通过（检查 4 只记录不判定）\n', toc);

%% 图
fig = figure('Visible', 'off', 'Position', [100 100 900 420]);
semilogx(cfg.fiber.dist_km, errF, '-s', cfg.fiber.dist_km, errN(:,1), '-o', ...
    cfg.fiber.dist_km, errN(:,2), '-o', cfg.fiber.dist_km, errN(:,3), '-o', cfg.fiber.dist_km, err8, '-^');
grid on; xlabel('距离 (km)'); ylabel('符号相对误差 (dB)');
legend('浮点 FFT-FDE', 'FNT b=16 (4,4)', 'FNT b=16 (5,5)', 'FNT b=16 (6,6)', 'FNT b=8 最优', 'Location', 'best');
title('FNT-TrDE 的符号误差（无噪声）');
exportgraphics(fig, fullfile(cfg.out.dir, 'm4_check.png'), 'Resolution', 150);
close(fig);

%% 本地函数
function e = symErr_dB(y, tx, cfg, rsym0, valid)
rs = rx_demod(y, tx.h, cfg);
e = 20*log10(sqrt(mean(abs(rs(valid) - rsym0(valid)).^2)) / sqrt(mean(abs(rsym0(valid)).^2)) + eps);
end

function o = scoreOf(y, tx, cfg)
[~, hb] = rx_demod(y, tx.h, cfg);
o = score_ber(hb, tx.bits, cfg);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
