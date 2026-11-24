%% run_ber_sweep  M6：扫距离和 Eb/N0，比较各色散补偿方案的 BER（对照 Xing 图 19 的趋势）
% 方案：理想逆色散、A 时域 FIR、B 浮点 FFT-FDE、B 定点 FFT-FDE（W=12/14/16）、
%       C FNT-TrDE（b=16 输入/系数 4/4 位与 6/6 位；b=8 取不溢出的最好组合）。
% 输出：outputs/m6_ber_vs_distance.png、m6_ber_vs_ebn0.png、m6_results.mat，并在屏幕打印表格。
% 同一个距离、同一个 Eb/N0 下所有方案吃同一份接收序列（同一份噪声），差异只来自补偿方法。
% 100 万比特 + 6 个距离 + 多个方案，运行几分钟。

addpath(fileparts(mfilename('fullpath')));  setup_paths;

cfg0 = cfg_fnt_cd();                         % 2^17 比特：只用来搜 b=8 的位宽
cfg = cfg0;
cfg.sig.nBits = 2^20;  cfg.sig.nSym = cfg.sig.nBits / cfg.sig.bitsPerSym;   % 100 万比特
nd = numel(cfg.fiber.dist);
EbWork = cfg.noise.EbN0_work_dB;
EbList = cfg.noise.EbN0_sweep_dB;
tic;

%% 0) b=8 的位宽：在小信号上无噪声搜索（不溢出且最大值 <= 90% 限值，符号误差最小）
rng(cfg0.seed);  tx0 = gen_tx(cfg0);
valid0 = (cfg0.rx.edgeSym+1) : (cfg0.sig.nSym-cfg0.rx.edgeSym);
rsym00 = rx_demod(tx0.x, tx0.h, cfg0);
bits8 = zeros(nd, 2);
for i = 1:nd
    z = cfg0.fiber.dist(i);  N = cfg0.blk.N(i);
    r0 = apply_cd(tx0.x, cfg0, z);
    e = inf;
    for xb = 3:5
        for hb = 2:4
            [y, inf_] = eq_fdc_fnt(r0, cfg0, z, N, struct('b', 8, 'xbit', xb, 'hbit', hb));
            if inf_.nOverflow == 0 && inf_.maxAbs <= 0.9*inf_.limit
                rs = rx_demod(y, tx0.h, cfg0);
                ee = 20*log10(sqrt(mean(abs(rs(valid0) - rsym00(valid0)).^2)) / sqrt(mean(abs(rsym00(valid0)).^2)));
                if ee < e, e = ee;  bits8(i, :) = [xb hb]; end
            end
        end
    end
    if bits8(i, 1) == 0, bits8(i, :) = [3 2]; end       % 找不到不溢出的组合时取最粗的，结果里标注
end
fprintf('b=8 位宽（xbit,hbit）：'); fprintf('(%d,%d) ', bits8.'); fprintf('\n');

%% 方案表
names = {'理想逆色散', 'A: FIR 2x', 'B: float FFT', 'B: FFT W=12', 'B: FFT W=14', 'B: FFT W=16', ...
         'C: FNT b16 4/4', 'C: FNT b16 6/6', 'C: FNT b8', 'B: FFT W=16 both', 'B: FFT W=12 both'};
run_scheme = @(k, r, z, N, i) run_one(k, r, cfg, z, N, i, bits8);
nS = numel(names);

%% 1) 距离扫描（固定 Eb/N0）
rng(cfg.seed);  tx = gen_tx(cfg);
ber1 = nan(nd, nS);  nErr1 = nan(nd, nS);
for i = 1:nd
    z = cfg.fiber.dist(i);  N = cfg.blk.N(i);
    rng(cfg.seed + i);
    r = add_noise(apply_cd(tx.x, cfg, z), cfg, EbWork);
    for k = 1:nS
        y = run_scheme(k, r, z, N, i);
        o = scoreOf(y, tx, cfg);
        ber1(i, k) = o.ber;  nErr1(i, k) = o.nErr;
    end
    fprintf('  距离 %4d km 完成（%.0f s）\n', cfg.fiber.dist_km(i), toc);
end
fprintf('\n[1] Eb/N0 = %g dB 时的 BER（%d 比特）\n%8s', EbWork, cfg.sig.nBits, 'dist/km');
fprintf('%15s', names{:});  fprintf('\n');
for i = 1:nd, fprintf('%8d', cfg.fiber.dist_km(i)); fprintf('%15.2e', ber1(i, :)); fprintf('\n'); end
fprintf('误码数\n');
for i = 1:nd, fprintf('%8d', cfg.fiber.dist_km(i)); fprintf('%15d', nErr1(i, :)); fprintf('\n'); end

%% 2) Eb/N0 扫描（最难的 4400 km）
i2 = nd;  z = cfg.fiber.dist(i2);  N = cfg.blk.N(i2);
ks = [1 3 4 5 6 7 8 9 10 11];                 % 不含 FIR（4400 km 的 FIR 太慢）
ber2 = nan(numel(EbList), nS);
for j = 1:numel(EbList)
    rng(cfg.seed + 100 + j);
    r = add_noise(apply_cd(tx.x, cfg, z), cfg, EbList(j));
    for k = ks
        y = run_scheme(k, r, z, N, i2);
        ber2(j, k) = scoreOf(y, tx, cfg).ber;
    end
    fprintf('  Eb/N0 %g dB 完成（%.0f s）\n', EbList(j), toc);
end
theory = 0.5 * erfc(sqrt(10.^(EbList/10)));
fprintf('\n[2] %d km 的 BER vs Eb/N0\n%8s %12s', cfg.fiber.dist_km(i2), 'Eb/N0', '理论');
fprintf('%15s', names{ks});  fprintf('\n');
for j = 1:numel(EbList), fprintf('%8.1f %12.2e', EbList(j), theory(j)); fprintf('%15.2e', ber2(j, ks)); fprintf('\n'); end

% 代价：达到目标 BER 所需的 Eb/N0，相对 B: float FFT 的差
tg = [1e-3 1e-4];
need = nan(numel(tg), nS);
for t = 1:numel(tg)
    for k = ks, need(t, k) = ebn0_at(EbList, ber2(:, k), tg(t)); end
end
fprintf('\n达到目标 BER 所需的 Eb/N0（dB）与相对 B: float FFT 的代价（dB）\n%8s', '目标BER');
fprintf('%15s', names{ks});  fprintf('\n');
for t = 1:numel(tg)
    fprintf('%8.0e', tg(t));  fprintf('%15.2f', need(t, ks));  fprintf('\n');
    fprintf('%8s', '代价');     fprintf('%15.2f', need(t, ks) - need(t, 3));  fprintf('\n');
end

save(fullfile(cfg.out.dir, 'm6_results.mat'), 'names', 'ber1', 'nErr1', 'ber2', 'EbList', 'EbWork', 'bits8', 'need', 'tg');

%% 图（不画 'both' 两条线：它们是论文字面缩放的试验，只在表里保留）
fl = @(b) max(b, 1e-7);
mk = {'-', '-', '-', ':', '--', '-', '--', '-', '-'};
col = [0 0 0; 0 .75 .75; 0 0 1; .75 0 .75; .75 0 .75; .75 0 .75; 0 .6 0; 0 .6 0; 1 0 0];
kp1 = 1:9;                                    % 距离图：理想、FIR、B、C
kp2 = ks(ks <= 9);                            % Eb/N0 图：同上，去掉 FIR
desc = {'理想逆色散：直接乘 G^{-1}，无量化，作上限', ...
        'A: FIR 2x：时域 FIR，抽头数取估计值的 2 倍', ...
        'B: float FFT：浮点 FFT 频域补偿（overlap-save）', ...
        'B: FFT W=12/14/16：定点 FFT，数据通路字长 12/14/16 位', ...
        'C: FNT b16 4/4、6/6：模数 65537 的 FNT，输入/滤波器系数分别量化为 4/4 位、6/6 位整数（xbit/hbit）', ...
        'C: FNT b8：模数 257 的 FNT，位宽取该距离下不溢出的最佳 (xbit,hbit)'};

fig = figure('Visible', 'off', 'Position', [100 100 1000 560]);
for k = kp1, semilogy(cfg.fiber.dist_km, fl(ber1(:, k)), mk{k}, 'Color', col(k,:), 'Marker', 'o', 'MarkerSize', 4); hold on; end
set(gca, 'XScale', 'log'); grid on;
xlabel('距离 (km)'); ylabel('BER'); lg = legend(names(kp1), 'Location', 'eastoutside');
title(sprintf('BER vs 距离（Eb/N0 = %g dB，%d 比特）', EbWork, cfg.sig.nBits));
add_note(fig, gca, lg, desc);
exportgraphics(fig, fullfile(cfg.out.dir, 'm6_ber_vs_distance.png'), 'Resolution', 150);
close(fig);

fig = figure('Visible', 'off', 'Position', [100 100 1000 560]);
semilogy(EbList, theory, 'k--'); hold on;
for k = kp2, semilogy(EbList, fl(ber2(:, k)), mk{k}, 'Color', col(k,:), 'Marker', 'o', 'MarkerSize', 4); end
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
lg = legend(['理论 QPSK', names(kp2)], 'Location', 'eastoutside');
title(sprintf('BER vs Eb/N0（%d km，N=%d）', cfg.fiber.dist_km(i2), N));
add_note(fig, gca, lg, desc([1 3:6]));
exportgraphics(fig, fullfile(cfg.out.dir, 'm6_ber_vs_ebn0.png'), 'Resolution', 150);
close(fig);
fprintf('\n用时 %.0f 秒；图和结果已保存到 %s\n', toc, cfg.out.dir);

%% 本地函数
function y = run_one(k, r, cfg, z, N, i, bits8)
switch k
    case 1, y = apply_cd(r, cfg, -z);
    case 2, y = eq_fir(r, cfg, z, 2*cfg.td.nTap_Rs(i));
    case 3, y = eq_fdc_fft(r, cfg, z, N);
    case 4, y = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', 12));
    case 5, y = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', 14));
    case 6, y = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', 16));
    case 7, y = eq_fdc_fnt(r, cfg, z, N, struct('b', 16, 'xbit', 4, 'hbit', 4));
    case 8, y = eq_fdc_fnt(r, cfg, z, N, struct('b', 16, 'xbit', 6, 'hbit', 6));
    case 9, y = eq_fdc_fnt(r, cfg, z, N, struct('b', 8, 'xbit', bits8(i,1), 'hbit', bits8(i,2)));
    case 10, y = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', 16, 'scheme', 'both'));   % 每层正逆向都截 1 位的字面理解
    case 11, y = eq_fdc_fft_fixed(r, cfg, z, N, struct('W', 12, 'scheme', 'both'));
end
end

function o = scoreOf(y, tx, cfg)
[~, hb] = rx_demod(y, tx.h, cfg);
o = score_ber(hb, tx.bits, cfg);
end

function e = ebn0_at(eb, ber, target)
% 在 log10(BER)-Eb/N0 曲线上线性插值，找首次降到 target 以下的 Eb/N0；够不到返回 NaN
lb = log10(max(ber(:), 1e-12));
k = find(lb <= log10(target), 1);
if isempty(k), e = NaN; return; end
if k == 1, e = eb(1); return; end
e = eb(k-1) + (log10(target) - lb(k-1)) / (lb(k) - lb(k-1)) * (eb(k) - eb(k-1));
end

function add_note(fig, ax, lg, desc)
% 在图下方加一块图例说明文字
drawnow;
ax.Position(2) = 0.30;  ax.Position(4) = 0.62;
lg.Position(2) = ax.Position(2);
annotation(fig, 'textbox', [0.03 0.01 0.94 0.20], 'String', ['图例说明：' desc], ...
    'EdgeColor', 'none', 'FontSize', 9, 'VerticalAlignment', 'top', 'Interpreter', 'tex');
end
