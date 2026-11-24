function d = fnt_mp_design(rcal, cfg, z, N, o)
%FNT_MP_DESIGN  多遍 FNT-TrDE（eq_fdc_fnt_mp）的离线设计：生成各遍的整数滤波器和还原权重。
%   d = fnt_mp_design(rcal, cfg, z, N, o)
%   rcal  标定用的接收信号（另一组随机数据，离线），只用来定每一遍的比例，不参与补偿
%   o     结构体（缺省值见 mp_defaults）：
%           b     8      q = 2^b+1
%           sx    1      输入比例：xq = round(x*sx)
%           Xs    1      输入切片数（1 或 2）。2 片：xq = S*xh + xl，xl 在 [-S/2, S/2]
%           S     8      输入切片的基（2 的幂，重组只需移位）
%           K     1      滤波器按能量切成的时间段数
%           Hs    1      每段滤波器的"残差切片"数：第 1 片量化后的误差再放大量化成第 2 片……
%           hopt  true   带内优化整数滤波器（见 hq_inband_opt）
%           fill  0.9    每一遍在标定信号上的输出峰值占 +-(q-1)/2 的比例
%
%   三个思路（第 1 条单遍就有，第 2、3 条要多遍）：
%   1) 带限滤波器 + 带内优化：滤波器只在信号带宽内（|f| <= Rs(1+beta)/2，外侧 3 GHz 余弦滚降）
%      等于逆色散，带外置 0，抽头能量更集中；量化后再按信号带内加权误差逐抽头 +-1 微调。
%   2) 切片（精确扩展精度）：输入或滤波器拆成高、低两部分，各自做一遍 FNT 卷积，结果在
%      普通整数里按权重加回去。每一遍都不溢出，所以结果精确。
%   3) 时间分段：滤波器按能量切成 K 段，每段抽头少，卷积结果小，可以用更大的比例。
%   每一遍（输入片 t, 滤波器段 k, 残差片 s）各自定比例 a，使该遍在标定信号上的输出峰值
%   = fill*(q-1)/2；最终 y = sum wx(t) * (x_t (*) hN) / a / sx。

o = mp_defaults(o);
lim = (2^o.b) / 2;
Xc = mp_blocks(rcal, N);
xs = mp_xslice(Xc, o);
h = cd_band_filter(cfg, z, N);
n = (0:N-1).'; n(n >= N/2) = n(n >= N/2) - N;           % 以 0 为中心的时间下标
[~, ord] = sort(n); e = cumsum(abs(h(ord)).^2); e = e / e(end);
seg = zeros(N, 1); seg(ord) = min(o.K, 1 + floor(e * o.K - 1e-12));
w = inband_weight(cfg, N);
d.o = o;  d.hN = {};  d.a = [];  d.t = [];
for t = 1:numel(xs)
    FX = fft(xs{t});
    for k = 1:o.K
        res = h .* (seg == k);
        for s = 1:o.Hs
            c0 = ifft(FX .* fft(res));
            pk = max(abs([real(c0(:)); imag(c0(:))]));
            if pk < 1e-9, break; end                         % 该输入片全为 0（sx 太小），这一遍省掉
            a = o.fill * lim / pk;
            for tries = 1:40
                hN = round(res * a);
                if o.hopt, hN = hq_inband_opt(hN, res * a, w); end
                cc = round(ifft(FX .* fft(hN)));
                if max(abs([real(cc(:)); imag(cc(:))])) <= o.fill * lim, break; end
                a = a * 0.95;
            end
            d.hN{end+1} = hN;  d.a(end+1) = a;  d.t(end+1) = t;
            res = res - hN / a;                              % 剩下的误差交给下一片
        end
    end
end
end

function o = mp_defaults(opts)
o.b = 8;  o.sx = 1;  o.Xs = 1;  o.S = 8;  o.K = 1;  o.Hs = 1;  o.hopt = true;  o.fill = 0.9;
f = fieldnames(opts);
for k = 1:numel(f), o.(f{k}) = opts.(f{k}); end
end

function w = inband_weight(cfg, N)
% 带内误差权重：RRC 功率谱的平方（误差经接收匹配滤波后才影响判决），加一个很小的底
k = (0:N-1).'; k(k >= ceil(N/2)) = k(k >= ceil(N/2)) - N; f = abs(k) * cfg.sig.fs / N;
a = cfg.sig.Rs*(1-cfg.pulse.rolloff)/2;  b = cfg.sig.Rs*(1+cfg.pulse.rolloff)/2;
r2 = ones(N, 1); r2(f >= b) = 0; m = f > a & f < b; r2(m) = 0.5*(1+cos(pi*(f(m)-a)/(b-a)));
w = r2.^2 + 1e-3;
end
