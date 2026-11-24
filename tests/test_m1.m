%% test_m1  M1 验收：发射、色散信道、加噪、接收判决
% 检查项：
%   1) 无色散无噪声：BER = 0
%   2) 不补偿：各距离 BER 明显变差，星座散开
%   3) 理想逆色散（apply_cd 取 -z）无噪声：BER = 0
%   4) 加噪：无色散 与 色散+理想补偿 两条 BER 曲线都落在 Q(sqrt(2Eb/N0)) 上
% 图不显示，保存到 outputs/。

addpath(fileparts(fileparts(mfilename('fullpath'))));  setup_paths;   % 自动加入项目路径

cfg = cfg_fnt_cd();
rng(cfg.seed);
tx = gen_tx(cfg);
fprintf('发射信号：%d 个样本，平均功率 %.4f\n', numel(tx.x), mean(abs(tx.x).^2));

%% 1) 无色散无噪声
o = run_rx(tx.x, tx, cfg);
fprintf('[1] 无色散无噪声：BER = %g（%d 个比特）\n', o.ber, o.nUsed);
assert(o.ber == 0, 'M1 检查 1 失败');

%% 2) 不补偿 / 3) 理想补偿（无噪声）
nd = numel(cfg.fiber.dist);
berRaw = zeros(1, nd); berIdeal = zeros(1, nd); maxErr = zeros(1, nd);
rsymRaw = cell(1, nd);
for i = 1:nd
    z = cfg.fiber.dist(i);
    r = apply_cd(tx.x, cfg, z);
    oRaw = run_rx(r, tx, cfg);
    berRaw(i) = oRaw.ber;  rsymRaw{i} = oRaw.rsym;
    rc = apply_cd(r, cfg, -z);
    oId = run_rx(rc, tx, cfg);
    berIdeal(i) = oId.ber;
    maxErr(i) = max(abs(rc - tx.x));
    fprintf('[2/3] %5d km：不补偿 BER = %.4f；理想补偿 BER = %g，max|rc-x| = %.2e\n', ...
        cfg.fiber.dist_km(i), berRaw(i), berIdeal(i), maxErr(i));
end
assert(all(berIdeal == 0) && all(maxErr < 1e-9), 'M1 检查 3 失败');
assert(berRaw(end) > 0.3, 'M1 检查 2 失败：4400 km 不补偿应接近随机');

%% 4) 加噪后与理论对比
EbN0 = 0:1:8;
theory = 0.5 * erfc(sqrt(10.^(EbN0/10)));      % QPSK：Q(sqrt(2Eb/N0))
berA = zeros(size(EbN0));  berB = zeros(size(EbN0));
zTest = cfg.fiber.dist(end);                    % 4400 km
for i = 1:numel(EbN0)
    nA = add_noise(tx.x, cfg, EbN0(i));
    berA(i) = run_rx(nA, tx, cfg).ber;
    nB = add_noise(apply_cd(tx.x, cfg, zTest), cfg, EbN0(i));
    berB(i) = run_rx(apply_cd(nB, cfg, -zTest), tx, cfg).ber;
end
fprintf('\n%8s %12s %12s %12s\n', 'Eb/N0', '理论', '无色散', '色散+理想补偿');
fprintf('%8.1f %12.3e %12.3e %12.3e\n', [EbN0; theory; berA; berB]);

%% 作图
fig = figure('Visible', 'off', 'Position', [100 100 1100 800]);
sel = [1 3 6];
for j = 1:3
    subplot(2, 3, j);
    s = rsymRaw{sel(j)}(cfg.rx.edgeSym+1 : cfg.rx.edgeSym+3000);
    plot(real(s), imag(s), '.', 'MarkerSize', 3); axis equal; grid on;
    title(sprintf('不补偿 %d km', cfg.fiber.dist_km(sel(j))));
end
subplot(2, 3, 4:6);
semilogy(EbN0, theory, 'k-', EbN0, max(berA, 1e-6), 'bo', EbN0, max(berB, 1e-6), 'rx');
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('理论 Q(\surd(2E_b/N_0))', '无色散', sprintf('%d km + 理想补偿', cfg.fiber.dist_km(end)), ...
    'Location', 'southwest');
title(sprintf('M1 检查（%d 比特，最低可见 BER \\approx %.0e）', cfg.sig.nBits, 1/cfg.sig.nBits));
exportgraphics(fig, fullfile(cfg.out.dir, 'm1_check.png'), 'Resolution', 150);
close(fig);
fprintf('\n图已保存：%s\n', fullfile(cfg.out.dir, 'm1_check.png'));

%% 本地函数：接收 + 评分（接收机只拿到 r 和 h；tx.bits 只用于评分）
function o = run_rx(r, tx, cfg)
[rsym, hatBits] = rx_demod(r, tx.h, cfg);
o = score_ber(hatBits, tx.bits, cfg);
o.rsym = rsym;
end
