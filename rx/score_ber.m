function out = score_ber(hatBits, bits, cfg)
%SCORE_BER  评分：把判决比特与发送比特（标准答案）比对，统计 BER（M1）。
%   out = score_ber(hatBits, bits, cfg)
%
%   这一步只在仿真里存在：真实接收机不知道发了什么。
%   统计时丢弃首尾各 cfg.rx.edgeSym 个符号，避开圆卷积和滤波的边缘效应。
%
%   输出 out：ber, nErr, nUsed（参与统计的比特数）

nSym = cfg.sig.nSym;
edge = cfg.rx.edgeSym;

idxSym  = (edge+1) : (nSym-edge);
idxBits = reshape([2*idxSym-1; 2*idxSym], [], 1);

out.nErr  = sum(hatBits(idxBits) ~= bits(idxBits));
out.nUsed = numel(idxBits);
out.ber   = out.nErr / out.nUsed;
end
