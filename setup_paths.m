function setup_paths()
%SETUP_PATHS  把本项目的子文件夹加入 MATLAB 搜索路径。
%   在项目根目录运行一次 setup_paths 即可；tests/ 下的脚本会自己调用它。
%
%   MATLAB 按"文件名"找函数，所以整个项目里函数名不能重名；FNT 相关的函数
%   统一用 fnt_ 前缀，就是为了避免和别处的函数冲突。

root = fileparts(mfilename('fullpath'));
addpath(root, ...
        fullfile(root, 'channel'), ...     % 发射、色散信道、加噪
        fullfile(root, 'rx'), ...          % 接收判决、评分
        fullfile(root, 'equalizer'), ...   % 色散补偿方案 A、B、C
        fullfile(root, 'fnt'), ...         % FNT 核心
        fullfile(root, 'complexity'), ...  % 复杂度计数（M7）
        fullfile(root, 'tests'));          % 验收脚本
end
