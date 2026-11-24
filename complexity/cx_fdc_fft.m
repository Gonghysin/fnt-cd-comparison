function [FA, parts] = cx_fdc_fft(N, W, kind)
%CX_FDC_FFT  方案 B 一个 overlap-save 块的复杂度，单位 1 位全加器（M7，Xing 第 VI 节规则）。
%   [FA, parts] = cx_fdc_fft(N, W, kind)
%   N 变换点数；W 数据位宽；kind 'R2'（radix-2）或 'SR'（split-radix）。
%
%   一个块 = 正向 FFT + N 个复数逐点乘 + 逆向 IFFT。规则：
%     复数乘 = 3 次实数乘 + 3 次实数加；复数加 = 2 次实数加；
%     1 次 W 位实数乘 = W/3 次 W 位实数加；W 位实数加 = W 个全加器。
%     所以 1 次复数乘 = (W+3) 次实数加。
%   变换里的复数乘只数"非平凡"的（旋转因子不是 1、-1、±j）：
%     R2：  (N/2)log2N - N + 1
%     SR：  M(N) = M(N/2) + 2M(N/4) + 2(N/4 - 1)，M(1)=M(2)=M(4)=0（W8 按一次完整复数乘算）
%   复数加两种都是 N*log2N 次/次变换。逐点乘 N 次复数乘（滤波器是常数，同样按 W 位乘）。
%   这两条计数式与 Xing 表 VIII 的 FFT 行（12～16 位，N=64～2048）逐位吻合，见 test_m7。

logN = log2(N);
switch kind
    case 'R2', cm = N/2*logN - N + 1;
    case 'SR', cm = sr_cmul(N);
    otherwise, error('kind 只能是 R2 或 SR');
end
cadd = N * logN;
radd_per_cmul = W + 3;                           % W/3*3 + 3
perFFT = cm * radd_per_cmul + cadd * 2;          % 以 W 位实数加为单位
pw = N * radd_per_cmul;
realAdds = 2*perFFT + pw;
FA = realAdds * W;
parts = struct('cmulPerFFT', cm, 'caddPerFFT', cadd, 'realAdds', realAdds);
end

function m = sr_cmul(N)
if N <= 4, m = 0; return; end
m = sr_cmul(N/2) + 2*sr_cmul(N/4) + 2*(N/4 - 1);
end
