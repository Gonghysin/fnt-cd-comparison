function c = cx_twcost(v, b)
%CX_TWCOST  模 q=2^b+1 下乘常数 v 的成本，单位 b 位加法（Xing 第 IV 节）。
%   v 为列向量（0..q-1）。分三类：
%     v = ±2^k           -> 0（循环移位 + 取反）
%     v = ±2^k * sqrt2   -> 1（两次移位 + 一次减法，k 为奇数次幂的情形）
%     其他常数           -> b/3（CSD 平均）
q = 2^b + 1;
p = 1;  P = zeros(2*b, 1);
for t = 1:2*b, P(t) = p;  p = mod(p*2, q); end         % 2^t，t=0..2b-1（含负号：2^b = -1）
s2 = mod(2^(b/4) * (2^(b/2) - 1), q);
S = mod(P * s2, q);
c = (b/3) * ones(size(v));
c(ismember(v, S)) = 1;
c(ismember(v, P)) = 0;
end
