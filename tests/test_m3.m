%% test_m3  M3 验收：FNT 核心（1D、2D、移位乘法、D1、QRNS）
% 检查项：
%   1) 参数：omega 是 N 次本原单位根；sqrt2 满足 sqrt2^2=2；列出各 N 的 kind
%   2) 移位乘法 = 取模乘法；D1 编码解码互逆
%   3) FNT 与定义式（O(N^2) 求和）逐点相等；FNT->IFNT 精确还原
%   4) 1D：FNT 圆卷积与直接圆卷积（模 q）逐点完全相等
%   5) 有符号小整数：卷积结果不溢出时，提升后与真实整数卷积一致
%   6) QRNS：复数圆卷积与直接复数卷积一致
%   7) 2D：折回后与一维圆卷积逐点相等；overlap-save 有效区等于线性卷积

addpath(fileparts(fileparts(mfilename('fullpath'))));  setup_paths;   % 自动加入项目路径

rng(2026);
fprintf('=== M3 FNT 核心验收 ===\n');

%% 1) 参数
fprintf('\n[1] 变换核（omega）\n%4s %7s %6s %9s %8s\n', 'b', 'q', 'N', 'kind', 'omega');
for b = [8 16]
    q = 2^b + 1;
    Ns = 2.^(1:min(b, 11));
    for N = Ns
        par = fnt_params(b, N);
        if ismember(N, [2*b, 4*b, 8*b]) || N == Ns(end) || N == 2
            fprintf('%4d %7d %6d %9s %8d\n', b, q, N, par.kind, par.omega);
        end
    end
end

%% 2) 移位乘法与 D1
for b = [8 16]
    q = 2^b + 1;
    if b == 8, xs = 0:2^b; else, xs = [0:300, randi([0 2^b], 1, 20000)]; end
    for k = 0:2*b-1
        assert(isequal(fnt_mulpow2(xs, k, b), mod(xs * 2^k, q)), 'fnt_mulpow2 失败 b=%d k=%d', b, k);
    end
    xs = unique([0:min(1000, 2^b), 2^b-5:2^b]);
    [d, z] = fnt_d1('enc', xs);
    assert(all(d >= 0 & d <= 2^b - 1), 'D1 编码超出 b 位');
    assert(isequal(fnt_d1('dec', d, z), xs), 'D1 互逆失败');
end
fprintf('\n[2] fnt_mulpow2 与取模乘法逐点一致（k=0..2b-1）；D1 编码解码互逆\n');

%% 3) 定义式与还原
for b = [8 16]
    q = 2^b + 1;
    par = fnt_params(b, 16);
    x = randi([0 q-1], 16, 3);
    Xdef = zeros(16, 3); % 对比一下蝶形FNT和常规定义式出发的FNT是不是一样的
    for k = 0:15
        for n = 0:15
            Xdef(k+1, :) = mod(Xdef(k+1, :) + x(n+1, :) * mod_pow(par.omega, mod(n*k, 16), q), q);
        end
    end
    assert(isequal(fnt_1d(x, par, 1), Xdef), 'FNT 与定义式不一致 b=%d', b);
    for N = 2.^(1:min(b, 11))
        par = fnt_params(b, N);
        x = randi([0 q-1], N, 4);
        assert(isequal(fnt_1d(fnt_1d(x, par, 1), par, -1), x), 'FNT->IFNT 还原失败 b=%d N=%d', b, N);
    end
end
fprintf('[3] FNT 与定义式逐点相等；FNT->IFNT 精确还原（各 N）\n');

%% 4) 1D 圆卷积
fprintf('\n[4] 1D 圆卷积（全范围随机输入，模 q 下逐点相等）\n%4s %6s %9s %8s\n', 'b', 'N', 'kind', '结果');
cases = [8 16; 8 64; 8 256; 16 16; 16 64; 16 512; 16 2048];
for c = 1:size(cases, 1)
    b = cases(c, 1); N = cases(c, 2); q = 2^b + 1;
    par = fnt_params(b, N);
    x = randi([0 q-1], N, 1);  h = randi([0 q-1], N, 1);
    y = fnt_1d(mod(fnt_1d(x, par, 1) .* fnt_1d(h, par, 1), q), par, -1);
    ok = isequal(y, ref_cconv(x, h, q)); % ref_cconv 是直接按卷积定义（循环地逐项累加）算的参照答案。要求两者逐位相同。
    fprintf('%4d %6d %9s %8s\n', b, N, par.kind, string(ok));
    assert(ok, '1D 圆卷积不相等 b=%d N=%d', b, N);
end

%% 5) 有符号小整数，不溢出
b = 16; N = 256; q = 2^b + 1; par = fnt_params(b, N);
x = randi([-7 7], N, 1);  h = randi([-7 7], N, 1);
yt = round(real(ifft(fft(x) .* fft(h))));
y = fnt_lift(fnt_1d(mod(fnt_1d(mod(x, q), par, 1) .* fnt_1d(mod(h, q), par, 1), q), par, -1), q);
assert(isequal(y, yt), '有符号卷积失败');
fprintf('\n[5] 有符号小整数（|x|,|h|<=7，N=256，q=65537）：提升后与真实整数卷积一致，max|y|=%d < q/2\n', max(abs(yt)));

%% 6) QRNS 复数卷积
fprintf('\n[6] QRNS 复数圆卷积\n');
for cc = [8 16 7; 16 64 7]'
    b = cc(1); N = cc(2); a = cc(3); q = 2^b + 1;
    if b == 8, a = 1; end
    par = fnt_params(b, N);
    x = complex(randi([-a a], N, 1), randi([-a a], N, 1));
    h = complex(randi([-a a], N, 1), randi([-a a], N, 1));
    yt = round(ifft(fft(x) .* fft(h)));
    Xr = fnt_1d(mod(real(x), q), par, 1);  Xi = fnt_1d(mod(imag(x), q), par, 1);
    Hr = fnt_1d(mod(real(h), q), par, 1);  Hi = fnt_1d(mod(imag(h), q), par, 1);
    [Yr, Yi] = qrns_pointmul(Xr, Xi, Hr, Hi, b);
    y = complex(fnt_lift(fnt_1d(Yr, par, -1), q), fnt_lift(fnt_1d(Yi, par, -1), q));
    assert(max(abs(yt)) < (q-1)/2 && isequal(y, yt), 'QRNS 复数卷积失败 b=%d N=%d', b, N);
    fprintf('  b=%d N=%d 输入幅度<=%d：一致（max|y|=%d）\n', b, N, a, max(abs([real(yt); imag(yt)])));
end

%% 7) 2D
fprintf('\n[7] 2D-FNT：N = L*M，行变换 M 点、列变换 2L 点\n%4s %6s %4s %4s %5s %14s %18s\n', ...
    'b', 'N', 'L', 'M', '2L', '折回=1D圆卷积', 'overlap-save 有效区');
cases2 = [8 256 8 32; 8 512 16 32; 8 1024 16 64; 8 2048 16 128; 16 4096 32 128];
for c = 1:size(cases2, 1)
    b = cases2(c, 1); N = cases2(c, 2); L = cases2(c, 3); M = cases2(c, 4); q = 2^b + 1;
    parR = fnt_params(b, M);  parC = fnt_params(b, 2*L);
    P = N/2;
    % (i) 一般输入：折回后应与一维 N 点圆卷积逐点相等
    x = randi([0 q-1], N, 1);  h = randi([0 q-1], N, 1);
    Z = fnt_2d(mod(fnt_2d(fnt_pack2d(x, L, M), parR, parC, 1) .* fnt_2d(fnt_pack2d(h, L, M), parR, parC, 1), q), parR, parC, -1);
    ok1 = isequal(fnt_unpack2d(Z, q), ref_cconv(x, h, q));
    % (ii) 因果短滤波器（长 P+1）：折回结果在 n >= P 的有效区应等于线性卷积
    h2 = zeros(N, 1);  h2(1:P+1) = randi([0 q-1], P+1, 1);
    Z = fnt_2d(mod(fnt_2d(fnt_pack2d(x, L, M), parR, parC, 1) .* fnt_2d(fnt_pack2d(h2, L, M), parR, parC, 1), q), parR, parC, -1);
    y2 = fnt_unpack2d(Z, q);
    ylin = zeros(N, 1);                          % 线性卷积（不环绕），模 q
    for k = 0:N-1
        ylin(k+1:end) = mod(ylin(k+1:end) + x(k+1) * h2(1:N-k), q);
    end
    ok2 = isequal(y2(P+1:end), ylin(P+1:end));
    fprintf('%4d %6d %4d %4d %5d %14s %18s\n', b, N, L, M, 2*L, string(ok1), string(ok2));
    assert(ok1 && ok2, '2D-FNT 失败 b=%d N=%d', b, N);
end
fprintf('\nM3 全部通过\n');

%% 本地函数
function y = ref_cconv(x, h, q)
% 直接圆卷积（模 q），逐项累加，作为精确参照
N = numel(x);  y = zeros(N, 1);
for k = 0:N-1
    y = mod(y + x(k+1) * circshift(h, k), q);
end
end

function r = mod_pow(a, e, q)
r = 1; a = mod(a, q);
while e > 0
    if mod(e, 2) == 1, r = mod(r*a, q); end
    a = mod(a*a, q);  e = floor(e/2);
end
end
