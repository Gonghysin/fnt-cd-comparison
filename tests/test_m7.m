%% test_m7  M7 验收：复杂度计数
% 检查项：
%   1) FFT-FDE（R2、SR；12～16 位；N=64～2048）与 Xing 表 VIII 相差 < 0.5%（表里只有 3 位有效数字）
%   2) FNT 1D-HR 的情形（b=16 全部 N，b=8 的 N=64/128/256）与表 VIII 相差 < 6%
%   3) FNT b=8 的 2D 情形（N>=512）：只报告差距，不做断言（本仓库未实现 Xing 第 V-C 节的 2D 优化，是偏保守的上界）
%   4) 单项手算：N=64、W=12 的 R2 块 = 76392 个全加器（手算见 cx_fdc_fft 注释）

addpath(fileparts(fileparts(mfilename('fullpath'))));  setup_paths;
Ns = [64 128 256 512 1024 2048];
R2 = [7.64e4 1.82e5 4.21e5 9.59e5 2.15e6 4.77e6; 8.69e4 2.07e5 4.80e5 1.09e6 2.45e6 5.43e6; 9.81e4 2.33e5 5.42e5 1.23e6 2.77e6 6.14e6; ...
      1.10e5 2.62e5 6.07e5 1.38e6 3.10e6 6.88e6; 1.22e5 2.91e5 6.76e5 1.54e6 3.46e6 7.67e6];
SR = [5.59e4 1.33e5 3.09e5 7.03e5 1.58e6 3.50e6; 6.32e4 1.51e5 3.49e5 7.96e5 1.79e6 3.96e6; 7.10e4 1.69e5 3.93e5 8.95e5 2.01e6 4.46e6; ...
      7.92e4 1.89e5 4.38e5 9.99e5 2.24e6 4.98e6; 8.78e4 2.09e5 4.86e5 1.11e6 2.49e6 5.52e6];
F8  = [2.00e4 4.65e4 1.05e5 2.89e5 6.43e5 1.42e6];
F16 = [4.27e4 1.04e5 2.46e5 5.56e5 1.22e6 2.65e6];

fprintf('=== M7 复杂度计数 ===\n');
assert(cx_fdc_fft(64, 12, 'R2') == 76392, 'M7 检查 4 失败');
for i = 1:5, for k = 1:6
    a = cx_fdc_fft(Ns(k), 11+i, 'R2');  b = cx_fdc_fft(Ns(k), 11+i, 'SR');
    assert(abs(a/R2(i,k) - 1) < 0.005 && abs(b/SR(i,k) - 1) < 0.005, 'M7 检查 1 失败：W=%d N=%d', 11+i, Ns(k));
end, end
fprintf('[1] FFT-FDE R2/SR 与表 VIII 全部 60 项相差 < 0.5%%\n');
e8 = zeros(1, 6);  e16 = e8;
for k = 1:6
    e8(k)  = cx_fdc_fnt(Ns(k), 8)  / F8(k)  - 1;
    e16(k) = cx_fdc_fnt(Ns(k), 16) / F16(k) - 1;
end
assert(all(abs(e16) < 0.06) && all(abs(e8(1:3)) < 0.06), 'M7 检查 2 失败');
fprintf('[2] FNT 1D-HR 相对表 VIII：b=16 %s%%；b=8(N=64..256) %s%%\n', mat2str(100*e16, 3), mat2str(100*e8(1:3), 3));
fprintf('[3] FNT b=8 的 2D 情形（N=512..2048）相对表 VIII：%s%%（未做 2D 优化，偏保守）\n', mat2str(100*e8(4:6), 3));
fprintf('M7 通过\n');
