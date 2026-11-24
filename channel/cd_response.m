function G = cd_response(cfg, z, n)
%CD_RESPONSE  长度 n 的 FFT 频率网格上的色散传递函数（M2）。
%   G = cd_response(cfg, z, n)
%   G(z,w) = exp(-j*beta2*z/2 * w^2)，w = 2*pi*k*fs/n，k 按 FFT 顺序折叠为
%   0..n/2-1, -n/2..-1。z 取负值得到理想逆色散，用作补偿器的频率响应。
%
%   色散公式只在这一个文件里：apply_cd（信道）、eq_fir、eq_fdc_fft（补偿器）都调用它。
%   n 不必等于信号长度：信道用信号长度，FIR 用 2^16，FFT-FDE 用块长 N。

k = (0:n-1).';                                      % 末尾的 .' 是转置，把行向量变成列向量
k(k >= ceil(n/2)) = k(k >= ceil(n/2)) - n;          % FFT 频率序号 0..n/2-1, -n/2..-1
w = 2*pi * k * cfg.sig.fs / n;                      % 角频率 rad/s

G = exp(-1j * cfg.fiber.beta2 * z / 2 * w.^2);
end
