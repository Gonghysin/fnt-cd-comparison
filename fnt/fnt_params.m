function par = fnt_params(b, N)
%FNT_PARAMS  费马数变换（FNT）的参数：模数 q = 2^b+1、长度 N、变换核 omega（M3）。
%   par = fnt_params(b, N)
%   b   位宽（8 -> q=257，16 -> q=65537）；N 变换长度，2 的幂，且 N <= 2^b。
%
%   变换核 omega 必须是模 q 下的 N 次本原单位根（阶恰为 N），且乘 omega 的幂要尽量是移位：
%     N <= 2b       : omega = 2^(2b/N)        （2 的阶是 2b）           kind = 'pow2'
%     2b < N <= 4b  : omega = (sqrt2)^(4b/N)  （sqrt2 的阶是 4b）       kind = 'sqrt2'
%     N > 4b        : omega = 3^((q-1)/N)     （3 是 q 的本原根）       kind = 'general'
%   其中 sqrt2 = 2^(b/4) * (2^(b/2) - 1) mod q，满足 sqrt2^2 = 2 mod q（Xing 第 II 节）。
%   kind = 'general' 的变换里，乘变换核需要真乘法（对应论文里的"乘法数"）。
%
%   输出字段：b, q, N, logN, omega, omegaInv, Ninv, kind,
%            W (omega^0..omega^(N/2-1))、Winv（逆变换用）、rev（位反转下标，0 基）。
%   所有运算用 double 存整数，q < 2^17，乘积 < 2^34，远小于 2^53，结果是精确的。

assert(b >= 4 && mod(b, 4) == 0, 'b 需为 4 的倍数（sqrt2 的公式要求 b/4 为整数）');
q = 2^b + 1;
assert(N >= 2 && N == 2^round(log2(N)) && N <= 2^b, 'N 需为 2 的幂且 N <= 2^b');

% 选出合适的\omega 本质有点像查表
if N <= 2*b
    omega = modpow(2, 2*b/N, q);   kind = 'pow2';
elseif N <= 4*b
    s2 = mod(2^(b/4) * (2^(b/2) - 1), q);
    assert(mod(s2*s2, q) == 2, 'sqrt2 公式验证失败');
    omega = modpow(s2, 4*b/N, q);  kind = 'sqrt2';
else
    omega = modpow(3, (q-1)/N, q); kind = 'general';
end
assert(modpow(omega, N, q) == 1 && modpow(omega, N/2, q) == q-1, 'omega 不是 N 次本原单位根');

par.b = b;  par.q = q;  par.N = N;  par.logN = log2(N);
par.omega    = omega;
par.omegaInv = modpow(omega, q-2, q);       % 费马小定理求逆 费马小定理告诉我们：$a^{q-1}\equiv1$，所以 $a\cdot a^{q-2}\equiv1$，即 $a^{q-2}$ 就是 $a$ 的逆。
par.Ninv     = modpow(N, q-2, q);
par.kind     = kind;

% 旋转因子表 W 。 Wi 是逆变换用的
W = ones(N/2, 1);  Wi = ones(N/2, 1); 
for k = 2:N/2
    W(k)  = mod(W(k-1)  * par.omega,    q);
    Wi(k) = mod(Wi(k-1) * par.omegaInv, q);
end
par.W = W;  par.Winv = Wi;

% 位反转表rev
% 降低FNT这个变换本身的运算复杂度用的。不理解也没事
nb = par.logN;  rev = zeros(N, 1);
for i = 0:N-1
    v = i; r = 0;
    for t = 1:nb
        r = bitor(bitshift(r, 1), bitand(v, 1));
        v = bitshift(v, -1);
    end
    rev(i+1) = r;
end
par.rev = rev;
end

function r = modpow(a, e, q)
% r = $a^e \bmod q$
r = 1; a = mod(a, q);
while e > 0
    if mod(e, 2) == 1, r = mod(r*a, q); end
    a = mod(a*a, q);
    e = floor(e/2);
end
end
