function y = fnt_mulpow2(x, k, b)
%FNT_MULPOW2  模 q=2^b+1 下乘 2^k，只用"移位 + 一次减法"实现（M3）。
%   y = fnt_mulpow2(x, k, b)    x 为 [0, 2^b] 内的整数（可以是数组）
%
%   原理：2^b = -1 (mod q)，所以 2^k 的周期是 2b：
%     k >= b 时   2^k = -2^(k-b)，先算 2^(k-b) 再取负；
%     k' < b 时   把 x 拆成高、低两段：x = hi*2^(b-k') + lo，
%                 x*2^k' = hi*2^b + lo*2^k' = lo*2^k' - hi   （循环左移，高位取反回填）
%   即：低 b-k' 位左移 k' 位，减去溢出到高位的部分。这就是"反相循环移位"。

q = 2^b + 1;
k = mod(k, 2*b);
neg = k >= b;
if neg, k = k - b; end
lo = mod(x, 2^(b-k));
hi = floor(x / 2^(b-k));
y = lo * 2^k - hi;
if neg, y = -y; end
y = mod(y, q);
end
