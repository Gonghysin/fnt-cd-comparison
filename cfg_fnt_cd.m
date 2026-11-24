function cfg = cfg_fnt_cd()
%CFG_FNT_CD  FNT / FFT / 时域 三种静态色散补偿对比仿真的参数集中处（M0）。
%   cfg = cfg_fnt_cd()  返回参数结构体；不带输出参数运行时打印一张汇总表。
%
%   参数来源：Xing 等, JLT 42(15) 5190 (2024) 第 VI 节；时域抽头数公式来自
%   余建军《基于数字信号处理的相干光通信技术》式 3-70。
%   标注"待定"的项目是论文里没读到明确值、暂取的临时值，见 README.md 第 6 节。
%
%   所有量用国际单位（Hz, s, m），仅在名字里带后缀的才例外（_km, _dB）。

%% 随机种子
cfg.seed = 2026;                       % 脚本入口处 rng(cfg.seed)

%% 信号
cfg.sig.Rb         = 56e9;             % 比特率 bit/s（Xing 第 VI 节：56 Gbit/s）
cfg.sig.bitsPerSym = 2;                % QPSK
cfg.sig.Rs         = cfg.sig.Rb / cfg.sig.bitsPerSym;   % 符号率 28 GBaud（推算值，论文未直接写）
cfg.sig.sps        = 2;                % 每符号采样数 2 SPS ( sampers per symbol )
cfg.sig.fs         = cfg.sig.sps * cfg.sig.Rs;          % 采样率 56 GSa/s
cfg.sig.Ts         = 1 / cfg.sig.fs;                    % 采样间隔
cfg.sig.nBits      = 2^17;             % 发送的比特总数 : 对应 Xing 里 PRBS（伪随机二进制序列）的长度（Xing：2^17）
cfg.sig.nSym       = cfg.sig.nBits / cfg.sig.bitsPerSym; % 符号总数

%% 脉冲成形（待定）
cfg.pulse.rolloff = 0.2;               % 滚降系数 $\beta$ 待定：论文未见，暂与 pol_demux_cma 一致
cfg.pulse.span    = 20;                % 滤波器长度 RRC 长度（符号数）2 SPS 时对应 $20\times2+1=41$ 个抽头

%% 接收端统计
cfg.rx.edgeSym = 1000;                 % BER 统计丢弃首尾各 1000 个符号（避开圆卷积/滤波边缘；> nTap_edge/sps）

%% 光纤色散
cfg.fiber.c     = 299792458;          % 光速 m/s
cfg.fiber.lambda = 1550e-9;            % 波长 m
cfg.fiber.D_ps   = 17;                 % 色散参数 ps/(nm*km) 1 nm 的光谱宽度在 1 km 光纤里会产生 17 ps 的时延差。
cfg.fiber.D      = cfg.fiber.D_ps * 1e-6;   % 转成 s/m^2：1 ps/(nm*km) = 1e-6 s/m^2
cfg.fiber.beta2  = -cfg.fiber.D * cfg.fiber.lambda^2 / (2*pi*cfg.fiber.c);  % s^2/m 群速度色散系数 $\beta_2$，s²/m，与 $D$ 的关系：$\beta_2=-\dfrac{D\lambda^2}{2\pi c}$
cfg.fiber.dist_km = [80 240 560 1120 2240 4400];     % Xing 第 VI 节
cfg.fiber.dist    = cfg.fiber.dist_km * 1e3;         % m

% 色散对信号的作用就是在频域乘一个纯相位项（Xing 式 13）：
% $$G(z,\omega)=\exp\!\Big(j\frac{D\lambda^2 z}{4\pi c}\omega^2\Big)=\exp\!\Big(-j\frac{\beta_2 z}{2}\omega^2\Big)$$

%% 噪声（待定）
% Xing 第 VI 节没有读到明确的 OSNR/SNR。先扫 Eb/N0，选浮点 BER 约 1e-4 的点作工作点。
cfg.noise.EbN0_sweep_dB = 5:1:10;      % M6：BER-Eb/N0 曲线的扫描范围（高于 10 dB 时 BER<4e-6，100 万比特里误码太少）
cfg.noise.EbN0_work_dB  = 8;           % M6 定：理论 BER 约 2e-4（Eb/N0 8 dB），100 万比特里约 200 个误码

%% 分块（overlap-save）
% 与距离一一对应：80/240/560/1120/2240/4400 km -> N = 64/128/.../2048
cfg.blk.N       = [64 128 256 512 1024 2048];        % 变换长度
cfg.blk.overlap = 0.5;                               % 重叠比例（Xing：50%）
cfg.blk.P       = cfg.blk.N * cfg.blk.overlap;       % 重叠样本数
cfg.blk.B       = cfg.blk.N - cfg.blk.P;             % 每块有效（新增）样本数

%% 时域 FIR 的抽头数
% (1) 余建军式 3-70：N = 2*floor(|D|*lambda^2*z / (2*c*T^2)) + 1，T 取采样间隔 Ts。
%     书里说这是抽头数的"上界"，对应信号占满整个奈奎斯特带宽 fs 的情形。
cfg.td.nTap_bound = 2 * floor(abs(cfg.fiber.D) * cfg.fiber.lambda^2 .* cfg.fiber.dist ...
                              ./ (2 * cfg.fiber.c * cfg.sig.Ts^2)) + 1;
% (2) 按信号实际带宽 B 估计：色散造成的时延展宽 = |D|*lambda^2*z*B/c，换算成样本数再 +1
%     （B = fs 时与式 3-70 一致，这是我自己按式 3-69 推的，不是书里的公式）。
%     B 取符号率 Rs（等效奈奎斯特带宽）和带边 Rs*(1+rolloff) 两种。
dl = abs(cfg.fiber.D) * cfg.fiber.lambda^2 .* cfg.fiber.dist / cfg.fiber.c;   % s/Hz
cfg.td.nTap_Rs   = ceil(dl * cfg.sig.Rs                       * cfg.sig.fs) + 1; % 默认取值，M2 里先用它
cfg.td.nTap_edge = ceil(dl * cfg.sig.Rs*(1+cfg.pulse.rolloff) * cfg.sig.fs) + 1;

%% FFT-FDE 位宽
cfg.fft.bits     = 12:16;              % 定点位宽扫描
cfg.fft.floatRef = true;               % 同时跑浮点 FFT 作参考

%% FNT 参数：F_n = 2^b + 1，b = 2^n（n=3 -> F=257，n=4 -> F=65537）
cfg.fnt.n   = [3 4];
cfg.fnt.b   = 2 .^ cfg.fnt.n;          % [8 16]
cfg.fnt.Fn  = 2 .^ cfg.fnt.b + 1;      % [257 65537]
% 无乘法变换核的最大长度（Xing 表 II 与第 II 节）：核为 2 时 N<=2^(n+1)，核为 sqrt(2) 时 N<=2^(n+2)
cfg.fnt.NmaxPow2    = 2 .^ (cfg.fnt.n + 1);
cfg.fnt.NmaxSqrt2   = 2 .^ (cfg.fnt.n + 2);

%% 复杂度计数规则（Xing 第 VI 节）
cfg.cx.RM_per_complexMul = 3;          % 复数乘 = 3 次实数乘 + 3 次实数加
cfg.cx.RA_per_complexMul = 3;
cfg.cx.RA_per_complexAdd = 2;          % 复数加 = 2 次实数加
cfg.cx.RA_per_RM_factor  = 1/3;        % 1 次 b 位实数乘 折合 b/3 次 b 位实数加
% b 位实数加 = b 个 1 位全加器（所以最终单位是 1 位全加器）

%% 输出
here = fileparts(mfilename('fullpath'));
cfg.out.dir = fullfile(here, 'outputs');
if ~exist(cfg.out.dir, 'dir'); mkdir(cfg.out.dir); end

%% 一致性检查
% 1) 重叠样本数 P 需不小于时域等效滤波器长度 - 1，否则 overlap-save 去重叠后有残余混叠。
%    以 nTap_Rs 为准（信号带宽按符号率）；nTap_edge（带边）只作提示。
cfg.chk.overlapOK     = cfg.blk.P >= (cfg.td.nTap_Rs   - 1);
cfg.chk.overlapEdgeOK = cfg.blk.P >= (cfg.td.nTap_edge - 1);
% 2) FNT 的 1D 变换长度上限是 2^b（模数 F_n 的限制）。行：n=3,4；列：各个 N。
%    为 false 的位置不能直接用 1D-FNT，需要 2D-FNT（Xing 第 V-C 节）。
cfg.chk.fnt1D = cfg.blk.N <= (2 .^ cfg.fnt.b(:));

%% 无输出参数时打印汇总
if nargout == 0
    fprintf('=== cfg_fnt_cd 汇总 ===\n');
    fprintf('符号率 %.1f GBaud，采样率 %.1f GSa/s，PRBS %d 比特 (%d 符号)\n', ...
        cfg.sig.Rs/1e9, cfg.sig.fs/1e9, cfg.sig.nBits, cfg.sig.nSym);
    fprintf('beta2 = %.4e s^2/m（D = %g ps/nm/km，lambda = %g nm）\n', ...
        cfg.fiber.beta2, cfg.fiber.D_ps, cfg.fiber.lambda*1e9);
    fprintf('FNT: F_3=%d(b=%d), F_4=%d(b=%d)；无乘法核最大长度 2^(n+1)=%s, 2^(n+2)=%s\n', ...
        cfg.fnt.Fn(1), cfg.fnt.b(1), cfg.fnt.Fn(2), cfg.fnt.b(2), ...
        mat2str(cfg.fnt.NmaxPow2), mat2str(cfg.fnt.NmaxSqrt2));
    fprintf('\n%8s %6s %6s %6s | %8s %8s %8s | %8s %8s\n', 'dist/km', 'N', 'P', 'B', ...
        '式3-70', '带宽Rs', '带边', 'P够(Rs)', 'P够(带边)');
    for k = 1:numel(cfg.fiber.dist_km)
        fprintf('%8d %6d %6d %6d | %8d %8d %8d | %8s %8s\n', cfg.fiber.dist_km(k), ...
            cfg.blk.N(k), cfg.blk.P(k), cfg.blk.B(k), cfg.td.nTap_bound(k), ...
            cfg.td.nTap_Rs(k), cfg.td.nTap_edge(k), ...
            string(cfg.chk.overlapOK(k)), string(cfg.chk.overlapEdgeOK(k)));
    end
    fprintf('\n1D-FNT 能否直接用（N <= 2^b）：\n%8s', 'N');
    fprintf('%6d', cfg.blk.N); fprintf('\n');
    for i = 1:numel(cfg.fnt.b)
        fprintf('%8s', sprintf('b=%d', cfg.fnt.b(i)));
        fprintf('%6s', string(cfg.chk.fnt1D(i, :))); fprintf('\n');
    end
    clear cfg
end
end
