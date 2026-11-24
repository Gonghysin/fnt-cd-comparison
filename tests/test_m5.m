%% test_m5  M5 验收：FFT-FDE 的位精确（定点）版
% 检查项：
%   1) 定点模型正确：字长很大（W=24）时，结果与浮点 eq_fdc_fft 几乎相同
%   2) 几乎无饱和：W=12..16 各距离 nSat<=10（输出峰值约占满量程 0.9）
%   3) 趋势：W 越大符号误差越小；截断与四舍五入的差别
%   4) BER（Eb/N0=6 dB，同一份噪声）：W=12..16 的定点 FFT，对照浮点 FFT 和 FNT(b=16, 6/6 位)
%   5) 缩放方案 'fwd' 与 'both'（Xing 字面说法）的差别
% 误差指标：补偿后符号率星座点与"无色散无噪声"参考的相对均方根误差（dB）。

addpath(fileparts(fileparts(mfilename('fullpath'))));  setup_paths;   % 自动加入项目路径

cfg = cfg_fnt_cd();
rng(cfg.seed);
tx = gen_tx(cfg);
nd = numel(cfg.fiber.dist);
valid = (cfg.rx.edgeSym+1) : (cfg.sig.nSym-cfg.rx.edgeSym);
rsym0 = rx_demod(tx.x, tx.h, cfg);
Ws = 12:16;
tic;

fprintf('=== M5 定点 FFT-FDE ===\n');

%% 1) 大字长时收敛到浮点
fprintf('\n[1] W=24 的定点结果与浮点 eq_fdc_fft 对比（符号误差 dB）\n%8s %6s %12s %12s %10s\n', 'dist/km', 'N', 'float', 'fixed W=24', 'nSat');
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    r = apply_cd(tx.x, cfg, z);
    eF = symErr_dB(eq_fdc_fft(r, cfg, z, N), tx, cfg, rsym0, valid);
    [y, inf_] = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', 24, 'rnd', true));
    e24 = symErr_dB(y, tx, cfg, rsym0, valid);
    fprintf('%8d %6d %12.1f %12.1f %10d\n', cfg.fiber.dist_km(i), N, eF, e24, inf_.nSat);
    assert(inf_.nSat == 0 && abs(e24 - eF) < 2, 'M5 检查 1 失败：W=24 未收敛到浮点');
end

%% 2)(3) W 扫描
errT = zeros(nd, numel(Ws));  errR = errT;  errF = zeros(nd, 1);
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    r = apply_cd(tx.x, cfg, z);
    errF(i) = symErr_dB(eq_fdc_fft(r, cfg, z, N), tx, cfg, rsym0, valid);
    for k = 1:numel(Ws)
        [y, inf_] = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', Ws(k), 'rnd', false));
        errT(i, k) = symErr_dB(y, tx, cfg, rsym0, valid);
        % 输出峰值约为满量程(xmax=3)的 0.86~0.9，W=12 截断时有极个别样本（约 260 万个里的 3 个）触顶，允许 <=10 个
        assert(inf_.nSat <= 10, 'M5 检查 2 失败：W=%d 距离 %d km 饱和 %d 次', Ws(k), cfg.fiber.dist_km(i), inf_.nSat);
        y = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', Ws(k), 'rnd', true));
        errR(i, k) = symErr_dB(y, tx, cfg, rsym0, valid);
    end
end
hdr = arrayfun(@(w) sprintf('W=%d', w), Ws, 'UniformOutput', false);
fprintf('\n[2] 符号误差(dB)：截断（trunc）\n%8s %6s', 'dist/km', 'N');  fprintf('%9s', hdr{:});  fprintf('%11s\n', 'float');
for i = 1:nd, fprintf('%8d %6d', cfg.fiber.dist_km(i), cfg.blk.N(i)); fprintf('%9.1f', errT(i,:)); fprintf('%11.1f\n', errF(i)); end
fprintf('\n[3] 符号误差(dB)：四舍五入（round）\n%8s %6s', 'dist/km', 'N');  fprintf('%9s', hdr{:});  fprintf('%11s\n', 'float');
for i = 1:nd, fprintf('%8d %6d', cfg.fiber.dist_km(i), cfg.blk.N(i)); fprintf('%9.1f', errR(i,:)); fprintf('%11.1f\n', errF(i)); end
assert(all(errT(:, end) < errT(:, 1) + 0.5), 'M5 检查 3 失败：字长增大后误差没有下降');

%% 4) BER
EbN0 = 6;  rng(cfg.seed + 1);
fprintf('\n[4] Eb/N0 = %g dB，同一份噪声；比特错误数（截断）\n%8s %10s %10s', EbN0, 'dist/km', '理想', 'float FFT');  fprintf('%10s', hdr{:});  fprintf('%14s\n', 'FNT b16(6,6)');
nErr = zeros(nd, numel(Ws) + 3);
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    r = add_noise(apply_cd(tx.x, cfg, z), cfg, EbN0);
    nErr(i, 1) = scoreOf(apply_cd(r, cfg, -z), tx, cfg).nErr;
    nErr(i, 2) = scoreOf(eq_fdc_fft(r, cfg, z, N), tx, cfg).nErr;
    for k = 1:numel(Ws)
        nErr(i, 2+k) = scoreOf(eq_fdc_fft_fixed(r, cfg, z, N, struct('W', Ws(k))), tx, cfg).nErr;
    end
    nErr(i, end) = scoreOf(eq_fdc_fnt(r, cfg, z, N, struct('b', 16, 'xbit', 6, 'hbit', 6)), tx, cfg).nErr;
    fprintf('%8d', cfg.fiber.dist_km(i));  fprintf('%10d', nErr(i, 1:end-1));  fprintf('%14d\n', nErr(i, end));
end

%% 5) 缩放方案
z = cfg.fiber.dist(6);  N = cfg.blk.N(6);
r = apply_cd(tx.x, cfg, z);
fprintf('\n[5] 缩放方案（4400 km，N=%d，符号误差 dB）\n%6s %10s %10s %8s\n', N, 'W', 'fwd', 'both', 'nSat(both)');
for W = [12 14 16 20]
    [y1, ~] = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', W, 'scheme', 'fwd'));
    [y2, i2] = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', W, 'scheme', 'both'));
    fprintf('%6d %10.1f %10.1f %8d\n', W, symErr_dB(y1, tx, cfg, rsym0, valid), symErr_dB(y2, tx, cfg, rsym0, valid), i2.nSat);
end

fprintf('\n用时 %.0f 秒\n', toc);

%% 图
fig = figure('Visible', 'off', 'Position', [100 100 900 420]);
plot(Ws, errT.', '-o'); hold on; set(gca, 'ColorOrderIndex', 1);
plot(Ws, errR.', '--x'); grid on;
xlabel('字长 W (bit)'); ylabel('符号相对误差 (dB)');
legend(arrayfun(@(d) sprintf('%d km', d), cfg.fiber.dist_km, 'UniformOutput', false), 'Location', 'northeast');
title('定点 FFT-FDE：实线=截断，虚线=四舍五入');
exportgraphics(fig, fullfile(cfg.out.dir, 'm5_check.png'), 'Resolution', 150);
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
