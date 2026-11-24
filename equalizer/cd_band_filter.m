function h = cd_band_filter(cfg, z, N)
%CD_BAND_FILTER  带限逆色散滤波器（N 点时域，零时刻在下标 1）。
%   信号带内 |f| <= Rs(1+beta)/2 等于 G^{-1}，外侧 3 GHz 内余弦滚降到 0，再往外为 0。
%   带外本来就会被接收匹配滤波器滤掉，所以对理想结果没有影响，但抽头能量更集中、卷积结果更小。
H = cd_response(cfg, -z, N);
k = (0:N-1).'; k(k >= ceil(N/2)) = k(k >= ceil(N/2)) - N; f = abs(k) * cfg.sig.fs / N;
f1 = cfg.sig.Rs*(1+cfg.pulse.rolloff)/2;  f2 = f1 + 3e9;
w = ones(N, 1); w(f >= f2) = 0; m = f > f1 & f < f2; w(m) = 0.5*(1+cos(pi*(f(m)-f1)/(f2-f1)));
h = ifft(H .* w);
end
