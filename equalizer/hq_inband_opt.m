function hN = hq_inband_opt(hN, t, w)
%HQ_INBAND_OPT  整数滤波器的带内优化：逐抽头试 +-1（实部、虚部），使
%   J = sum_f w(f) |FFT(hN - t)(f)|^2   减小，直到一整轮没有改动（坐标下降，离线）。
%   hN 初值（一般为 round(t)）；t 目标（实数比例下的理想系数）；w 频域权重（N x 1，>= 0）。
%   J 是二次型 e' Q e（Q 为循环矩阵，第一列 ifft(w)）；把 e_n 改 delta（delta = +-1 或 +-j）时
%   J 的变化 = 2*delta*Re(conj(part)*(Qe)_n) + Q_nn，只在变化为负时接受。
N = numel(hN);  qcol = ifft(w);  qcol = qcol(:);  q0 = real(qcol(1));
Qe = ifft(w .* fft(hN - t));
for sweep = 1:20
    changed = 0;
    for n = 1:N
        for part = [1 1i]
            g = real(conj(part) * Qe(n));
            if abs(g) > q0/2 + 1e-12
                dlt = -sign(g) * part;
                hN(n) = hN(n) + dlt;
                Qe = Qe + dlt * circshift(qcol, n-1);
                changed = changed + 1;
            end
        end
    end
    if changed == 0, break; end
end
end
