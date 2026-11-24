%% run_complexity  M7：按 Xing 第 VI 节规则数一个 overlap-save 块的 1 位全加器数，对照表 VIII
addpath(fileparts(mfilename('fullpath')));  setup_paths;
cfg = cfg_fnt_cd();
Ns = cfg.blk.N;
paper.R2 = [7.64e4 1.82e5 4.21e5 9.59e5 2.15e6 4.77e6;   % 12 bit
            8.69e4 2.07e5 4.80e5 1.09e6 2.45e6 5.43e6;
            9.81e4 2.33e5 5.42e5 1.23e6 2.77e6 6.14e6;
            1.10e5 2.62e5 6.07e5 1.38e6 3.10e6 6.88e6;
            1.22e5 2.91e5 6.76e5 1.54e6 3.46e6 7.67e6];
paper.SR = [5.59e4 1.33e5 3.09e5 7.03e5 1.58e6 3.50e6;
            6.32e4 1.51e5 3.49e5 7.96e5 1.79e6 3.96e6;
            7.10e4 1.69e5 3.93e5 8.95e5 2.01e6 4.46e6;
            7.92e4 1.89e5 4.38e5 9.99e5 2.24e6 4.98e6;
            8.78e4 2.09e5 4.86e5 1.11e6 2.49e6 5.52e6];
paper.F8  = [2.00e4 4.65e4 1.05e5 2.89e5 6.43e5 1.42e6];
paper.F16 = [4.27e4 1.04e5 2.46e5 5.56e5 1.22e6 2.65e6];
Ws = 12:16;
mine.R2 = zeros(5, 6);  mine.SR = mine.R2;
for i = 1:5, for k = 1:6
    mine.R2(i, k) = cx_fdc_fft(Ns(k), Ws(i), 'R2');
    mine.SR(i, k) = cx_fdc_fft(Ns(k), Ws(i), 'SR');
end, end
mine.F8 = zeros(1, 6);  mine.F16 = mine.F8;  inf8 = cell(1, 6);  inf16 = inf8;
for k = 1:6
    [mine.F8(k),  inf8{k}]  = cx_fdc_fnt(Ns(k), 8);
    [mine.F16(k), inf16{k}] = cx_fdc_fnt(Ns(k), 16);
end
fprintf('=== 每个 overlap-save 块的 1 位全加器数（本仓库 / 论文表 VIII / 相对差）===\n%-14s', 'N');  fprintf('%12d', Ns);  fprintf('\n');
pr = @(nm, a, p) fprintf('%-14s%s\n%-14s%s\n%-14s%s\n', nm, sprintf('%12.3g', a), '  论文', sprintf('%12.3g', p), '  相对差', sprintf('%11.1f%%', 100*(a-p)./p));
for i = 1:5
    pr(sprintf('%dbit R2', Ws(i)), mine.R2(i, :), paper.R2(i, :));
    pr(sprintf('%dbit SR', Ws(i)), mine.SR(i, :), paper.SR(i, :));
end
pr('FNT b=8', mine.F8, paper.F8);
pr('FNT b=16', mine.F16, paper.F16);
fprintf('\n每个 N 选中的方案（R 为基，L 为 2D 的参数）：\n');
for k = 1:6, fprintf('  N=%4d  b=8 %s R=%s L=%d | b=16 %s R=%s L=%d\n', Ns(k), inf8{k}.mode, mat2str(inf8{k}.R), inf8{k}.L, inf16{k}.mode, mat2str(inf16{k}.R), inf16{k}.L); end
fprintf('\n节省比例（相对 12bit R2 / 12bit SR；论文：b=8 约 72%% / 62%%，b=16 对 16bit SR 约 51%%）\n');
fprintf('%-26s', 'N'); fprintf('%8d', Ns); fprintf('\n');
fprintf('%-26s', 'FNT b=8  vs 12b R2');  fprintf('%7.0f%%', 100*(1 - mine.F8./mine.R2(1,:)));  fprintf('\n');
fprintf('%-26s', 'FNT b=8  vs 12b SR');  fprintf('%7.0f%%', 100*(1 - mine.F8./mine.SR(1,:)));  fprintf('\n');
fprintf('%-26s', 'FNT b=16 vs 16b SR');  fprintf('%7.0f%%', 100*(1 - mine.F16./mine.SR(5,:)));  fprintf('\n');
save(fullfile(cfg.out.dir, 'm7_results.mat'), 'Ns', 'Ws', 'mine', 'paper', 'inf8', 'inf16');

fig = figure('Visible', 'off', 'Position', [100 100 900 460]);
loglog(Ns, mine.R2(1,:), 'b-o', Ns, mine.SR(1,:), 'b--o', Ns, mine.SR(5,:), 'm--o', Ns, mine.F8, 'r-s', Ns, mine.F16, 'g-s', 'LineWidth', 1.2); hold on;
loglog(Ns, paper.F8, 'r:x', Ns, paper.F16, 'g:x', 'LineWidth', 1.2);
grid on; xlabel('变换长度 N'); ylabel('1 位全加器数 / 块');
legend({'12bit R2 FFT', '12bit SR FFT', '16bit SR FFT', 'FNT b=8（本仓库）', 'FNT b=16（本仓库）', 'FNT b=8（论文）', 'FNT b=16（论文）'}, 'Location', 'northwest');
title('复杂度对比（实线/虚线：本仓库计数；点线：论文表 VIII）');
exportgraphics(fig, fullfile(cfg.out.dir, 'm7_complexity.png'), 'Resolution', 150);  close(fig);
