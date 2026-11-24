%% test_m2  M2 验收：时域 FIR（A）与 FFT-FDE（B）的浮点版
% 检查项：
%   0) cd_response 是纯相位，且 G(z1)*G(z2) = G(z1+z2)
%   1) FIR 抽头数扫描（无噪声）：符号误差随抽头数下降并饱和
%   2) FFT-FDE（无噪声）：各距离用对应的 N，符号误差
%   3) 加噪（Eb/N0=6 dB，同一份噪声）：理想补偿 / FIR / FFT-FDE 的 BER 应基本一致
% 误差指标：补偿后符号率星座点与"无色散无噪声"参考的相对均方根误差（dB，越小越好）。

addpath(fileparts(fileparts(mfilename('fullpath'))));  setup_paths;   % 自动加入项目路径

cfg = cfg_fnt_cd();
rng(cfg.seed);
tx = gen_tx(cfg);
nd = numel(cfg.fiber.dist);
valid = (cfg.rx.edgeSym+1) : (cfg.sig.nSym-cfg.rx.edgeSym);
rsym0 = rx_demod(tx.x, tx.h, cfg);                 % 参考符号

%% 0) cd_response 的基本性质：纯相位，且 G(z1)*G(z2) = G(z1+z2)
n0 = 4096;
G1 = cd_response(cfg, 560e3, n0);  G2 = cd_response(cfg, -240e3, n0);
e0 = max(abs(G1 .* G2 - cd_response(cfg, 320e3, n0)));
fprintf('[0] max||G|-1| = %.2e；G(z1)G(z2)-G(z1+z2) 最大差 = %.2e\n', max(abs(abs(G1)-1)), e0);
assert(max(abs(abs(G1)-1)) < 1e-12 && e0 < 1e-9, 'M2 检查 0 失败');

%% 1) FIR 抽头数扫描（无噪声）
mult = [0.5 0.75 1 1.25 1.5 2];
errFIR = zeros(nd, numel(mult));
errBound = zeros(nd, 1);
nTapUsed = zeros(nd, numel(mult));
for i = 1:nd
    z = cfg.fiber.dist(i);
    r = apply_cd(tx.x, cfg, z);
    for m = 1:numel(mult)
        nTap = max(3, round(mult(m) * cfg.td.nTap_Rs(i)));
        [y, h] = eq_fir(r, cfg, z, nTap);
        nTapUsed(i, m) = numel(h);
        errFIR(i, m) = symErr_dB(y, tx, cfg, rsym0, valid);
    end
    y = eq_fir(r, cfg, z, cfg.td.nTap_bound(i));
    errBound(i) = symErr_dB(y, tx, cfg, rsym0, valid);
end
fprintf('\n[1] FIR 符号误差(dB)，列为抽头数 = 倍数 x nTap_Rs；最后一列为式3-70上界\n');
fprintf('%8s', 'dist/km'); fprintf('%10.3gx', mult); fprintf('%12s\n', '上界');
for i = 1:nd
    fprintf('%8d', cfg.fiber.dist_km(i)); fprintf('%11.1f', errFIR(i,:)); fprintf('%12.1f\n', errBound(i));
end
fprintf('（抽头数）\n');
for i = 1:nd
    fprintf('%8d', cfg.fiber.dist_km(i)); fprintf('%11d', nTapUsed(i,:)); fprintf('%12d\n', cfg.td.nTap_bound(i));
end

%% 2) FFT-FDE（无噪声）
errFFT = zeros(nd, 1);
for i = 1:nd
    z = cfg.fiber.dist(i);
    r = apply_cd(tx.x, cfg, z);
    y = eq_fdc_fft(r, cfg, z, cfg.blk.N(i));
    errFFT(i) = symErr_dB(y, tx, cfg, rsym0, valid);
end
% 注意：不要拿 y 和理想逆色散在"样本级"比 max 误差。N 点网格上的 H 与全带理想逆色散
% 只在带外分量上差别大（带外功率约 -46 dB，但处理方式不同，样本差能到 -10 dB），
% 这些分量会被匹配滤波滤掉，所以以符号级误差为准。
fprintf('\n[2] FFT-FDE（N 按距离取，overlap 50%%）\n%8s %6s %14s\n', 'dist/km', 'N', '符号误差(dB)');
for i = 1:nd
    fprintf('%8d %6d %14.1f\n', cfg.fiber.dist_km(i), cfg.blk.N(i), errFFT(i));
end

%% 3) 加噪后 BER（同一份噪声）
EbN0 = 6;
rng(cfg.seed + 1);
fprintf('\n[3] Eb/N0 = %g dB，同一份噪声；比特错误数（BER）\n%8s %18s %18s %18s %18s\n', ...
    EbN0, 'dist/km', '理想逆色散', 'FIR(1.5x)', 'FIR(2x)', 'FFT-FDE');
for i = 1:nd
    z = cfg.fiber.dist(i);
    r = add_noise(apply_cd(tx.x, cfg, z), cfg, EbN0);
    oI = scoreOf(apply_cd(r, cfg, -z), tx, cfg);
    oA = scoreOf(eq_fir(r, cfg, z, round(1.5*cfg.td.nTap_Rs(i))), tx, cfg);
    oA2 = scoreOf(eq_fir(r, cfg, z, round(2*cfg.td.nTap_Rs(i))), tx, cfg);
    oB = scoreOf(eq_fdc_fft(r, cfg, z, cfg.blk.N(i)), tx, cfg);
    fprintf('%8d %10d(%.2e) %10d(%.2e) %10d(%.2e) %10d(%.2e)\n', cfg.fiber.dist_km(i), ...
        oI.nErr, oI.ber, oA.nErr, oA.ber, oA2.nErr, oA2.ber, oB.nErr, oB.ber);
    assert(abs(oB.nErr - oI.nErr) <= max(15, 0.1*oI.nErr), 'M2 检查 3 失败：FFT-FDE 与理想补偿的误码数差异过大');
    assert(abs(oA2.nErr - oI.nErr) <= max(15, 0.1*oI.nErr), 'M2 检查 3 失败：FIR(2x) 与理想补偿的误码数差异过大');
end
assert(all(errFFT < -30), 'M2 检查 2 失败：FFT-FDE 符号误差过大');
assert(all(errFIR(:, mult==2) < -28), 'M2 检查 1 失败：FIR(2x) 符号误差过大');

%% 作图
fig = figure('Visible', 'off', 'Position', [100 100 1100 450]);
subplot(1, 2, 1);
plot(mult, errFIR.', '-o'); grid on;
xlabel('抽头数 / nTap\_Rs'); ylabel('符号相对误差 (dB)');
legend(arrayfun(@(d) sprintf('%d km', d), cfg.fiber.dist_km, 'UniformOutput', false), 'Location', 'northeast');
title('A. 时域 FIR：抽头数扫描（无噪声）');
subplot(1, 2, 2);
plot(cfg.fiber.dist_km, errFFT, '-s', cfg.fiber.dist_km, errFIR(:, mult==1.5), '-o', ...
    cfg.fiber.dist_km, errBound, '--^');
set(gca, 'XScale', 'log'); grid on;
xlabel('距离 (km)'); ylabel('符号相对误差 (dB)');
legend('B. FFT-FDE（N 随距离，overlap 50%）', 'A. FIR 1.5x nTap\_Rs', 'A. FIR 式3-70 上界', 'Location', 'best');
title('各方案符号误差（无噪声）');
exportgraphics(fig, fullfile(cfg.out.dir, 'm2_check.png'), 'Resolution', 150);
close(fig);
fprintf('\n图已保存：%s\n', fullfile(cfg.out.dir, 'm2_check.png'));

%% 本地函数
function e = symErr_dB(y, tx, cfg, rsym0, valid)
rs = rx_demod(y, tx.h, cfg);
e = 20*log10(sqrt(mean(abs(rs(valid) - rsym0(valid)).^2)) / sqrt(mean(abs(rsym0(valid)).^2)) + eps);
end

function o = scoreOf(y, tx, cfg)
[~, hb] = rx_demod(y, tx.h, cfg);
o = score_ber(hb, tx.bits, cfg);
end
