%% run_b8_improve  探索：怎样把 FNT b=8（q=257）的精度提上去，以及代价是多少
% 背景：M4/M6 里 b=8 直接量化（输入 3～4 位、系数 2～3 位）精度很差。瓶颈是动态范围：
%   有效区卷积结果必须落在 +-128 内，而色散滤波器很长，卷积结果 ~ 输入幅度 x 系数幅度 x sqrt(抽头数)。
% 这里试三类办法（实现见 fnt_mp_design / eq_fdc_fnt_mp）：
%   1) 单遍改进（不增加硬件）：带限滤波器、连续选比例、整数滤波器带内优化
%   2) 切片：输入或滤波器拆成高低两片，各做一遍 FNT，结果在普通整数里精确重组
%   3) 时间分段：滤波器按能量切成 K 段，每段单独一遍
% 流程：输入比例 sx 在另一组随机数据（标定数据）上选，不碰测试数据；
%       测试数据与 M6 完全相同（同一 tx、同一噪声），所以可以直接与 m6_results.mat 对比。
% 输出：outputs/b8_ber_vs_distance.png、b8_complexity.png、b8_results.mat，屏幕打印表格。
% 运行约 15～25 分钟（主要花在离线选比例和带内优化上）。

addpath(fileparts(mfilename('fullpath')));  setup_paths;
cfg0 = cfg_fnt_cd();                                   % 2^17 比特，用于标定
cfg = cfg0;  cfg.sig.nBits = 2^20;  cfg.sig.nSym = cfg.sig.nBits / cfg.sig.bitsPerSym;
nd = numel(cfg.fiber.dist);  Eb = cfg.noise.EbN0_work_dB;
m6 = load(fullfile(cfg.out.dir, 'm6_results.mat'));
tic;

% 方案：[Xs K Hs]（输入片数、时间段数、残差片数），遍数 = 乘积（空片会自动省掉）
cfgs = {'b8 改进·1 遍',              [1 1 1];
        'b8 2 遍（滤波器 2 片）',     [1 1 2];
        'b8 4 遍（2 段 x 2 片）',     [1 2 2];
        'b8 4 遍（输入 2 x 滤波器 2）', [2 1 2]};
nc = size(cfgs, 1);
sxGrid = 2.^(-1.5:0.25:5);

%% 标定数据（另一组随机比特和噪声）
rng(cfg.seed + 999);  txc = gen_tx(cfg0);
validc = (cfg0.rx.edgeSym+1) : (cfg0.sig.nSym - cfg0.rx.edgeSym);
rsc0 = rx_demod(txc.x, txc.h, cfg0);

%% 1) 每个方案、每个距离：在标定数据上选 sx，设计各遍滤波器
D = cell(nc, nd);  sxBest = nan(nc, nd);  eCal = nan(nc, nd);
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    rng(cfg.seed + 900 + i);  rcal = add_noise(apply_cd(txc.x, cfg0, z), cfg0, Eb);
    r0c = apply_cd(txc.x, cfg0, z);
    for c = 1:nc
        best = inf;
        for sx = sxGrid
            o = struct('sx', sx, 'Xs', cfgs{c,2}(1), 'K', cfgs{c,2}(2), 'Hs', cfgs{c,2}(3));
            d = fnt_mp_design(rcal, cfg0, z, N, o);
            if isempty(d.hN), continue; end
            e = symerr(quick_eval(r0c, cfg0, N, d), txc.h, cfg0, rsc0, validc);
            if e < best, best = e;  D{c,i} = d;  sxBest(c,i) = sx; end
        end
        eCal(c,i) = best;
    end
    fprintf('  标定 %4d km 完成（%.0f s）\n', cfg.fiber.dist_km(i), toc);
end

%% 2) 测试数据（与 M6 相同）：真正的 FNT 多遍补偿，BER + 诊断
rng(cfg.seed);  tx = gen_tx(cfg);
valid = (cfg.rx.edgeSym+1) : (cfg.sig.nSym - cfg.rx.edgeSym);
rs0 = rx_demod(tx.x, tx.h, cfg);
ber = nan(nd, nc);  eTest = ber;  nPass = ber;  nOver = ber;  nMis = ber;  FA = ber;
FA12R2 = nan(nd, 1);  FA12SR = FA12R2;  FAb8 = FA12R2;
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    rng(cfg.seed + i);  r = add_noise(apply_cd(tx.x, cfg, z), cfg, Eb);       % 与 run_ber_sweep 相同
    r0 = apply_cd(tx.x, cfg, z);
    for c = 1:nc
        d = D{c,i};
        [y, inf_] = eq_fdc_fnt_mp(r, cfg, z, N, d);
        [~, hb] = rx_demod(y, tx.h, cfg);  s = score_ber(hb, tx.bits, cfg);
        ber(i,c) = s.ber;  nOver(i,c) = inf_.nOverflow;  nMis(i,c) = inf_.nMismatch;  nPass(i,c) = inf_.passes;
        eTest(i,c) = symerr(eq_fdc_fnt_mp(r0, cfg, z, N, d), tx.h, cfg, rs0, valid);
        FA(i,c) = cx_fdc_fnt_mp(N, 8, max(d.t), inf_.passes);
    end
    FA12R2(i) = cx_fdc_fft(N, 12, 'R2');  FA12SR(i) = cx_fdc_fft(N, 12, 'SR');  FAb8(i) = cx_fdc_fnt(N, 8);
    fprintf('  测试 %4d km 完成（%.0f s）\n', cfg.fiber.dist_km(i), toc);
end
assert(all(nMis(:) == 0), 'FNT 结果与精确整数卷积不一致');

%% 打印
km = cfg.fiber.dist_km;
ref = {'理想逆色散', 1; 'B: float FFT', 3; 'B: FFT W=12', 4; 'B: FFT W=12 both', 11; 'C: FNT b16 6/6', 8; 'C: FNT b8（M6 原方案）', 9};
fprintf('\n[1] Eb/N0 = %g dB 的 BER（%d 比特，与 M6 同一份噪声）\n%-30s', Eb, cfg.sig.nBits, '方案');  fprintf('%10d', km);  fprintf('\n');
for k = 1:size(ref,1), fprintf('%-30s', ref{k,1});  fprintf('%10.2e', m6.ber1(:, ref{k,2}));  fprintf('\n'); end
for c = 1:nc, fprintf('%-30s', cfgs{c,1});  fprintf('%10.2e', ber(:,c));  fprintf('\n'); end
fprintf('\n[2] 无噪声符号误差（dB，测试数据）\n%-30s', '方案');  fprintf('%10d', km);  fprintf('\n');
for c = 1:nc, fprintf('%-30s', cfgs{c,1});  fprintf('%10.1f', eTest(:,c));  fprintf('\n'); end
fprintf('\n[3] 遍数 / 选中的 sx / 溢出数（应为 0）\n');
for c = 1:nc
    fprintf('%-30s 遍数 %s  sx %s  溢出 %s\n', cfgs{c,1}, mat2str(nPass(:,c).'), mat2str(round(sxBest(c,:)*100)/100), mat2str(nOver(:,c).'));
end
fprintf('\n[4] 每块 1 位全加器数，相对 12 位 SR-FFT 的比例（<1 表示更省）\n%-30s', '方案');  fprintf('%10d', km);  fprintf('\n');
fprintf('%-30s', 'C: FNT b8（M6 原方案）');  fprintf('%10.2f', FAb8 ./ FA12SR);  fprintf('\n');
for c = 1:nc, fprintf('%-30s', cfgs{c,1});  fprintf('%10.2f', FA(:,c) ./ FA12SR);  fprintf('\n'); end
fprintf('%-30s', '12 位 R2-FFT');  fprintf('%10.2f', FA12R2 ./ FA12SR);  fprintf('\n');

names = cfgs(:,1);
save(fullfile(cfg.out.dir, 'b8_results.mat'), 'names', 'cfgs', 'ber', 'eTest', 'eCal', 'nPass', 'sxBest', 'nOver', 'FA', 'FA12R2', 'FA12SR', 'FAb8', 'km', 'Eb');

%% 图
fl = @(b) max(b, 1e-7);
fig = figure('Visible', 'off', 'Position', [100 100 1000 600]);
st = {'k-', 'b-', 'm:', 'm-.', 'g-', 'r--'};
for k = 1:size(ref,1), semilogy(km, fl(m6.ber1(:, ref{k,2})), st{k}, 'Marker', 'o', 'MarkerSize', 4); hold on; end
cc = [1 .45 0; .85 .1 .1; .5 .2 .7; .1 .45 .85];
for c = 1:nc, semilogy(km, fl(ber(:,c)), '-', 'Color', cc(c,:), 'Marker', 's', 'MarkerSize', 5, 'LineWidth', 1.4); end
set(gca, 'XScale', 'log');  grid on;  xlabel('距离 (km)');  ylabel('BER');
lg = legend([ref(:,1); names], 'Location', 'eastoutside');
title(sprintf('FNT b=8 的改进（Eb/N0 = %g dB，%d 比特）', Eb, cfg.sig.nBits));
note = {'图例说明：', ...
        'B: FFT W=12：定点 FFT，正向每层右移 1 位；W=12 both：正逆向每层都右移 1 位（Xing 字面读法）', ...
        'C: FNT b8（M6 原方案）：直接量化，输入 3～4 位、系数 2～3 位', ...
        'b8 改进·1 遍：带限滤波器 + 连续选比例 + 整数滤波器带内优化，硬件不变', ...
        'b8 2 遍：滤波器拆成高低 2 片；2 段 x 2 片：滤波器再按时间切 2 段；输入 2 x 滤波器 2：输入也拆 2 片'};
drawnow;  ax = gca;  ax.Position(2) = 0.30;  ax.Position(4) = 0.62;  lg.Position(2) = ax.Position(2);
annotation(fig, 'textbox', [0.03 0.01 0.94 0.20], 'String', note, 'EdgeColor', 'none', 'FontSize', 9, 'VerticalAlignment', 'top');
exportgraphics(fig, fullfile(cfg.out.dir, 'b8_ber_vs_distance.png'), 'Resolution', 150);  close(fig);

fig = figure('Visible', 'off', 'Position', [100 100 900 460]);
semilogx(km, FAb8 ./ FA12SR, 'r--o'); hold on;
for c = 1:nc, semilogx(km, FA(:,c) ./ FA12SR, '-s', 'Color', cc(c,:), 'LineWidth', 1.4); end
semilogx(km, FA12R2 ./ FA12SR, 'b-o');  yline(1, 'k:', '12 位 SR-FFT');
grid on;  xlabel('距离 (km)');  ylabel('全加器数 / 12 位 SR-FFT');
legend(['C: FNT b8（M6 原方案）'; names; {'12 位 R2-FFT'}], 'Location', 'eastoutside');
title('每块复杂度（相对 12 位 split-radix FFT）');
exportgraphics(fig, fullfile(cfg.out.dir, 'b8_complexity.png'), 'Resolution', 150);  close(fig);
fprintf('\n用时 %.0f 秒\n', toc);

%% 本地函数
function e = symerr(y, h, cfg, rs0, valid)
rs = rx_demod(y, h, cfg);
e = 20*log10(norm(rs(valid) - rs0(valid)) / norm(rs0(valid)));
end

function y = quick_eval(r, cfg, N, d)
% 与 eq_fdc_fnt_mp 相同，但用 FFT 算精确整数卷积（不溢出时与 FNT 结果完全相同），只用于离线选比例
X = mp_blocks(r, N);  [xs, wx] = mp_xslice(X, d.o);
P = N/2;  B = N - P;  rows = P/2+1 : P/2+B;
acc = zeros(N, size(X, 2));
for p = 1:numel(d.hN)
    acc = acc + wx(d.t(p)) * round(ifft(fft(xs{d.t(p)}) .* fft(d.hN{p}))) / d.a(p);
end
y = reshape(acc(rows, :) / d.o.sx, [], 1);
end
