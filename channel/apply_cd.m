function y = apply_cd(x, cfg, z)
%APPLY_CD  在频域对信号施加光纤色散（纯相位项），M1。
%   y = apply_cd(x, cfg, z)
%   z 为传输距离（m）。z 取负值就是理想逆色散（用来做"完美补偿"的参照）。
%
%   传递函数（Xing 式 13）：
%       G(z,w) = exp(j*D*lambda^2*z/(4*pi*c) * w^2) = exp(-j*beta2*z/2 * w^2)
%   频率轴和 G 的计算在 cd_response.m 里。
%
%   这里用整段信号做一次 FFT，等价于**圆卷积**信道：信号首尾会互相"卷"到对方，
%   只影响两端各约 nTap/2 个样本。BER 统计时丢弃两端（cfg.rx.edgeSym）。

x = x(:);
G = cd_response(cfg, z, numel(x));
y = ifft(fft(x) .* G);
end
