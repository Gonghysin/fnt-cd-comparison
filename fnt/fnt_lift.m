function y = fnt_lift(x, q)
%FNT_LIFT  把模 q 的结果 [0, q-1] 还原成有符号整数 [-(q-1)/2, (q-1)/2]（M3）。
%   y = fnt_lift(x, q)
%   FNT 卷积只在"真实卷积结果的绝对值 < q/2"时才能还原出正确的有符号值，
%   否则发生溢出（这就是 Xing 里位宽约束 x_bit+h_bit+log2(M_d/2)+2 <= pb 的来源）。

y = x;
y(x > (q-1)/2) = x(x > (q-1)/2) - q;
end
