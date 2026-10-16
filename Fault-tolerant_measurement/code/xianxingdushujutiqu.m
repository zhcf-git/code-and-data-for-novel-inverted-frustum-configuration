%% ================================================================
%  六分量动态加载数据：830 Hz 数据峰值提取与六分量测量幅值计算
%
%  功能：
%  1）读取 830 Hz 六分量动态加载数据；
%  2）在 830 ± 1 Hz 范围内寻找 signal_Flc 的频谱峰值；
%  3）以输入力峰值频点作为实际激励频点；
%  4）提取该频点处 24 路响应复数频谱；
%  5）选取最接近实际峰值频率的 H_complex(:,:,idxH)；
%  6）计算校正后的六分量测量幅值；
%  7）按加载等级生成工作区变量。
%
%  工作区输出：
%       Result
%       LevelData
%       Result_1N、Result_2N、...
%       ActualAppliedLoad_1N、ActualAppliedLoad6_1N、FM_1N
%
%  说明：
%       ActualAppliedLoad：
%           6 × 1 实际施加载荷主分量
%           顺序：[Fx; Fy; Fz; Mx; My; Mz]
%           前三项单位 N，后三项单位 N.m
%
%       ActualAppliedLoad6：
%           6 × 6 实际施加载荷矩阵
%           行：激励分量 Fx、Fy、Fz、Mx、My、Mz
%           列：理论输入分量 Fx、Fy、Fz、Mx、My、Mz
%
%       FM：
%           6 × 6 六分量测量幅值矩阵
%           行：激励分量 Fx、Fy、Fz、Mx、My、Mz
%           列：测量分量 Fx、Fy、Fz、Mx、My、Mz
% ================================================================

clear;
clc;
close all;

%% ===================== 1. 基本参数 =====================

% ====== 数据路径：放置 87hz_fx_1n.mat 等文件的文件夹 ======
dataPath = 'D:\file\code\matlab\app\bd_freqall_260516\data\liner\830';
addpath(dataPath);

% ====== 标定矩阵文件 ======
% 如果 H_complex_jizhenqi_new.mat 不在 dataPath 或当前工作路径下，请改为绝对路径
HFile = 'H_q24_ch_01_02_03_04_05_06_07_08_09_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24.mat';

% ====== 只处理指定目标频率 / Hz ======
targetFreqOnly = 830;

% ====== 峰值搜索半带宽 / Hz ======
% 实际提取 signal_Flc 在 targetFreqOnly ± peakSearchHalfWidth 范围内的峰值
peakSearchHalfWidth = 1;

% ====== 截取时间 / s ======
T_start = 3;
T_end   = 13;

% ====== 24 路测力台响应通道 ======
num_channel = 24;

% ====== 输入力信号名称 ======
forceSignalName = 'signal_Flc';

% ====== 力矩加载力臂 / m ======
arm_Mx = 0.318;
arm_My = 0.318;
arm_Mz = 0.25;

% ====== 是否使用单边峰值幅值 ======
% 对正弦信号，FFT/N 后正频率幅值约为实际峰值的一半，因此这里乘 2
% 若需要 RMS 幅值，可在峰值基础上除以 sqrt(2)
useSingleSidedPeak = true;

% ====== 分量名称 ======
componentNames = {'Fx','Fy','Fz','Mx','My','Mz'};

%% ===================== 2. 加载 H_complex =====================

if exist(fullfile(dataPath, HFile), 'file')
    load(fullfile(dataPath, HFile));     % H_complex, freq_list 或只有 H_complex
else
    load(HFile);                         % H_complex, freq_list 或只有 H_complex
end

if ~exist('H_complex', 'var')
    error('H_complex_jizhenqi_new.mat 中缺少变量 H_complex。');
end

if size(H_complex, 1) ~= 6 || size(H_complex, 2) ~= num_channel
    error('H_complex 尺寸应为 6 × %d × N，当前为 %d × %d × %d。', ...
        num_channel, size(H_complex,1), size(H_complex,2), size(H_complex,3));
end

% 如果文件中没有 freq_list，则默认第三维对应 1:1000 Hz
if exist('freq_list', 'var')
    freq_list = freq_list(:);
else
    freq_list = (1:size(H_complex,3)).';
end

if length(freq_list) ~= size(H_complex,3)
    error('freq_list 长度与 H_complex 第三维频点数不一致。');
end

if targetFreqOnly < min(freq_list) || targetFreqOnly > max(freq_list)
    error('目标频率 %.3f Hz 超出 H_complex 标定频率范围 %.3f ~ %.3f Hz。', ...
        targetFreqOnly, min(freq_list), max(freq_list));
end

%% ===================== 3. 搜索数据文件 =====================

fileList = dir(fullfile(dataPath, '*.mat'));

validFiles = [];

for i = 1:length(fileList)

    name_i = fileList(i).name;

    % 排除标定矩阵文件
    if strcmpi(name_i, 'H_complex_jizhenqi_new.mat')
        continue;
    end

    % 匹配文件名，例如：830hz_fx_1n.mat
    token = regexp(name_i, ...
        '^(?<freq>\d+)\s*hz_(?<comp>fx|fy|fz|mx|my|mz)_(?<level>\d+)\s*n\.mat$', ...
        'names', 'ignorecase');

    if ~isempty(token)

        thisFreq = str2double(token.freq);

        % 只处理 830 Hz 数据
        if thisFreq ~= targetFreqOnly
            continue;
        end

        fileList(i).targetFreq = thisFreq;
        fileList(i).targetComp = upperFirst(token.comp);
        fileList(i).levelIndex = str2double(token.level);

        validFiles = [validFiles; fileList(i)];
    end
end

if isempty(validFiles)
    error('未找到符合命名规则的 %.0f Hz 文件，例如 %.0fhz_fx_1n.mat。', ...
        targetFreqOnly, targetFreqOnly);
end

fprintf('共找到 %d 个 %.0f Hz 数据文件。\n', length(validFiles), targetFreqOnly);

%% ===================== 4. 循环处理每个文件 =====================

Result = table;

wait_h = waitbar(0, '正在提取实际施加力峰值、理论施加载荷和六分量测量幅值...');

for fileIdx = 1:length(validFiles)

    fileName = validFiles(fileIdx).name;
    filePath = fullfile(dataPath, fileName);

    targetFreq = validFiles(fileIdx).targetFreq;
    targetComp = validFiles(fileIdx).targetComp;
    levelIndex = validFiles(fileIdx).levelIndex;

    data = load(filePath);

    %% ---------- 读取参考时间轴 ----------
    if ~isfield(data, 'signal_in1')
        error('文件 %s 中缺少 signal_in1。', fileName);
    end

    t_ref = getTimeVector(data.signal_in1);

    dt = mean(diff(t_ref));
    fs = 1 / dt;

    if fs <= 2 * targetFreq
        warning('文件 %s 的采样率 %.3f Hz 低于 2 × %.3f Hz，可能存在混叠风险。', ...
            fileName, fs, targetFreq);
    end

    idx_start = round(T_start * fs) + 1;
    idx_end   = round(T_end * fs);

    if idx_end > length(t_ref)
        error('文件 %s 数据长度不足 %.2f s。', fileName, T_end);
    end

    N = idx_end - idx_start + 1;

    % FFT 频率轴
    f_full = (0:N-1).' * fs / N;

    %% ---------- 提取输入力，并在 targetFreq ± 1 Hz 范围内寻找峰值 ----------
    if isfield(data, forceSignalName)
        ts_force = data.(forceSignalName);
    else
        error('文件 %s 中缺少输入力变量 %s。', fileName, forceSignalName);
    end

    forceData = getDataVector(ts_force);

    if idx_end > length(forceData)
        error('文件 %s 中输入力信号长度不足 %.2f s。', fileName, T_end);
    end

    f_cut = forceData(idx_start:idx_end);
    f_cut = f_cut - mean(f_cut);

    F_fft = fft(f_cut) / N;

    % 只在目标频率 ± peakSearchHalfWidth Hz 范围内找输入力频谱峰值
    freqBandIdx = find(f_full >= targetFreq - peakSearchHalfWidth & ...
                       f_full <= targetFreq + peakSearchHalfWidth);

    if isempty(freqBandIdx)
        error('文件 %s 中未找到 %.3f ± %.3f Hz 范围内的 FFT 频点。', ...
            fileName, targetFreq, peakSearchHalfWidth);
    end

    [~, localPeakIdx] = max(abs(F_fft(freqBandIdx)));
    idxFreqData = freqBandIdx(localPeakIdx);

    actualFreq_data = f_full(idxFreqData);
    freqErr_data = actualFreq_data - targetFreq;

    F_complex = F_fft(idxFreqData);

    if useSingleSidedPeak
        ActualAppliedForce_N = 2 * abs(F_complex);
    else
        ActualAppliedForce_N = abs(F_complex);
    end

    %% ---------- 根据输入力实际峰值频率选取对应 H 矩阵 ----------
    [freqErr_H, idxH] = min(abs(freq_list - actualFreq_data));
    actualFreq_H = freq_list(idxH);

    H_k = H_complex(:, :, idxH);      % 6 × 24

    fprintf('文件 %-20s | 目标 %.3f Hz | 输入力峰值 %.3f Hz | H矩阵 %.3f Hz | 力 %.6f N\n', ...
        fileName, targetFreq, actualFreq_data, actualFreq_H, ActualAppliedForce_N);

    %% ---------- 提取 24 路响应信号在输入力峰值频点处的复数频谱 ----------
    U_complex = zeros(num_channel, 1, 'like', 1i);

    for ch = 1:num_channel

        signalName = sprintf('signal_in%d', ch);

        if ~isfield(data, signalName)
            error('文件 %s 中缺少变量 %s。', fileName, signalName);
        end

        x = getDataVector(data.(signalName));

        if idx_end > length(x)
            error('文件 %s 中 %s 数据长度不足 %.2f s。', ...
                fileName, signalName, T_end);
        end

        x_cut = x(idx_start:idx_end);
        x_cut = x_cut - mean(x_cut);

        X_fft = fft(x_cut) / N;

        % 所有响应通道均取输入力峰值对应的同一频点
        U_complex(ch) = X_fft(idxFreqData);
    end

    %% ---------- 根据激励分量计算实际施加载荷 ----------
    % ActualAppliedLoad_one 顺序：
    % [Fx; Fy; Fz; Mx; My; Mz]
    %
    % Fx、Fy、Fz 单位：N
    % Mx、My、Mz 单位：N.m

    ActualAppliedLoad_one = zeros(6, 1);

    switch targetComp

        case 'Fx'
            ActualAppliedLoad_one(1) = ActualAppliedForce_N;

        case 'Fy'
            ActualAppliedLoad_one(2) = ActualAppliedForce_N;

        case 'Fz'
            ActualAppliedLoad_one(3) = ActualAppliedForce_N;

        case 'Mx'
            ActualAppliedLoad_one(4) = ActualAppliedForce_N * arm_Mx;

        case 'My'
            ActualAppliedLoad_one(5) = ActualAppliedForce_N * arm_My;

        case 'Mz'
            ActualAppliedLoad_one(6) = ActualAppliedForce_N * arm_Mz;

        otherwise
            error('未知激励分量：%s。', targetComp);
    end

    % 当前试验的实际主施加载荷值
    % 力激励时单位为 N，力矩激励时单位为 N.m
    ActualAppliedInputValue = max(abs(ActualAppliedLoad_one));

    %% ---------- 使用 H_complex 计算六分量测量幅值 ----------
    FM_complex = H_k * U_complex;     % 6 × 1

    if useSingleSidedPeak
        FM = 2 * abs(FM_complex);
    else
        FM = abs(FM_complex);
    end

    % FM 顺序：
    % FM(1) = Fx / N
    % FM(2) = Fy / N
    % FM(3) = Fz / N
    % FM(4) = Mx / N.m
    % FM(5) = My / N.m
    % FM(6) = Mz / N.m

    %% ---------- 写入总表 ----------
    oneRow = table;

    oneRow.FileName = string(fileName);

    oneRow.TargetFreq_Hz = targetFreq;
    oneRow.TargetComponent = string(targetComp);
    oneRow.LevelIndex = levelIndex;

    oneRow.ActualPeakFreq_Hz = actualFreq_data;
    oneRow.HMatrix_Freq_Hz = actualFreq_H;
    oneRow.FreqError_Data_Hz = freqErr_data;
    oneRow.FreqError_HMatrix_Hz = freqErr_H;

    oneRow.ActualAppliedForce_N = ActualAppliedForce_N;
    oneRow.ActualAppliedInputValue = ActualAppliedInputValue;

    oneRow.Input_Fx_N  = ActualAppliedLoad_one(1);
    oneRow.Input_Fy_N  = ActualAppliedLoad_one(2);
    oneRow.Input_Fz_N  = ActualAppliedLoad_one(3);
    oneRow.Input_Mx_Nm = ActualAppliedLoad_one(4);
    oneRow.Input_My_Nm = ActualAppliedLoad_one(5);
    oneRow.Input_Mz_Nm = ActualAppliedLoad_one(6);

    oneRow.Fx_Amp_N  = FM(1);
    oneRow.Fy_Amp_N  = FM(2);
    oneRow.Fz_Amp_N  = FM(3);
    oneRow.Mx_Amp_Nm = FM(4);
    oneRow.My_Amp_Nm = FM(5);
    oneRow.Mz_Amp_Nm = FM(6);

    Result = [Result; oneRow];

    waitbar(fileIdx / length(validFiles), wait_h, ...
        sprintf('正在处理 %s，进度 %d / %d', ...
        fileName, fileIdx, length(validFiles)));
end

close(wait_h);

%% ===================== 5. 总表排序 =====================

compOrder = categorical(Result.TargetComponent, ...
    {'Fx','Fy','Fz','Mx','My','Mz'}, ...
    'Ordinal', true);

Result.CompOrder = compOrder;

Result = sortrows(Result, ...
    {'TargetFreq_Hz','LevelIndex','CompOrder'});

Result.CompOrder = [];

%% ===================== 6. 按加载等级整理到工作区 =====================

levelList = unique(Result.LevelIndex).';

LevelData = struct();

for li = 1:numel(levelList)

    levelID = levelList(li);

    % 当前等级的数据
    T_level = Result(Result.LevelIndex == levelID, :);

    % 当前等级内按激励分量排序
    compOrder = categorical(T_level.TargetComponent, ...
        {'Fx','Fy','Fz','Mx','My','Mz'}, ...
        'Ordinal', true);

    T_level.CompOrder = compOrder;
    T_level = sortrows(T_level, 'CompOrder');
    T_level.CompOrder = [];

    % ============================================================
    % 实际施加载荷主分量，6 × 1
    %
    % 顺序：[Fx; Fy; Fz; Mx; My; Mz]
    %
    % 前三项单位为 N，后三项单位为 N.m
    % ============================================================
    ActualAppliedLoad = T_level.ActualAppliedInputValue;

    % ============================================================
    % 实际施加载荷完整矩阵，6 × 6
    %
    % 行：激励分量 Fx、Fy、Fz、Mx、My、Mz
    % 列：理论输入分量 Fx、Fy、Fz、Mx、My、Mz
    %
    % 理想情况下只有主对角线非零
    % ============================================================
    ActualAppliedLoad6 = [ ...
        T_level.Input_Fx_N, ...
        T_level.Input_Fy_N, ...
        T_level.Input_Fz_N, ...
        T_level.Input_Mx_Nm, ...
        T_level.Input_My_Nm, ...
        T_level.Input_Mz_Nm];

    % ============================================================
    % 测量六分量幅值矩阵，6 × 6
    %
    % 行：激励分量 Fx、Fy、Fz、Mx、My、Mz
    % 列：测量分量 Fx、Fy、Fz、Mx、My、Mz
    % ============================================================
    FM = [ ...
        T_level.Fx_Amp_N, ...
        T_level.Fy_Amp_N, ...
        T_level.Fz_Amp_N, ...
        T_level.Mx_Amp_Nm, ...
        T_level.My_Amp_Nm, ...
        T_level.Mz_Amp_Nm];

    rowComponent = cellstr(T_level.TargetComponent);
    colComponent = {'Fx','Fy','Fz','Mx','My','Mz'};

    % 保存到结构体，方便工作区展开查看
    fieldName = sprintf('Level_%dN', levelID);

    LevelData.(fieldName).ResultTable = T_level;

    LevelData.(fieldName).ActualAppliedLoad = ActualAppliedLoad;
    LevelData.(fieldName).ActualAppliedLoad6 = ActualAppliedLoad6;

    LevelData.(fieldName).FM = FM;

    LevelData.(fieldName).RowComponent = rowComponent;
    LevelData.(fieldName).ColComponent = colComponent;

    % 同时在工作区生成单独变量，方便直接双击查看
    assignin('base', sprintf('Result_%dN', levelID), T_level);
    assignin('base', sprintf('ActualAppliedLoad_%dN', levelID), ActualAppliedLoad);
    assignin('base', sprintf('ActualAppliedLoad6_%dN', levelID), ActualAppliedLoad6);
    assignin('base', sprintf('FM_%dN', levelID), FM);
end

%% ===================== 7. 显示结果 =====================

fprintf('\n处理完成！\n');
fprintf('已在工作区生成以下主要变量：\n');
fprintf('1）Result：全部试验汇总表\n');
fprintf('2）LevelData：按加载等级整理后的结构体\n');

for li = 1:numel(levelList)
    levelID = levelList(li);
    fprintf('3）Result_%dN、ActualAppliedLoad_%dN、ActualAppliedLoad6_%dN、FM_%dN\n', ...
        levelID, levelID, levelID, levelID);
end

fprintf('\nActualAppliedLoad_%dN 的含义：\n', levelList(1));
fprintf('顺序为 [Fx; Fy; Fz; Mx; My; Mz]\n');
fprintf('前三项单位为 N，后三项单位为 N.m\n\n');

fprintf('ActualAppliedLoad6_%dN 的矩阵含义：\n', levelList(1));
fprintf('行：激励分量 Fx、Fy、Fz、Mx、My、Mz\n');
fprintf('列：理论输入分量 Fx、Fy、Fz、Mx、My、Mz\n\n');

fprintf('FM_%dN 的矩阵含义：\n', levelList(1));
fprintf('行：激励分量 Fx、Fy、Fz、Mx、My、Mz\n');
fprintf('列：测量分量 Fx、Fy、Fz、Mx、My、Mz\n\n');

disp(Result);

%% ================================================================
%  局部函数
% ================================================================

function t = getTimeVector(ts)
% 从 timeseries 或结构体中读取时间向量

    if isa(ts, 'timeseries')
        t = ts.Time(:);
    elseif isstruct(ts) && isfield(ts, 'Time')
        t = ts.Time(:);
    else
        error('无法识别信号时间格式。信号应为 timeseries 或包含 Time 字段的结构体。');
    end
end

function x = getDataVector(ts)
% 从 timeseries 或结构体中读取数据向量

    if isa(ts, 'timeseries')
        x = ts.Data;
    elseif isstruct(ts) && isfield(ts, 'Data')
        x = ts.Data;
    else
        error('无法识别信号数据格式。信号应为 timeseries 或包含 Data 字段的结构体。');
    end

    x = squeeze(x);

    if size(x, 2) > 1 && size(x, 1) == 1
        x = x.';
    end

    if size(x, 2) > 1
        x = x(:, 1);
    end

    x = x(:);
end

function compName = upperFirst(compName)
% 将 fx 转为 Fx，mz 转为 Mz

    compName = lower(compName);

    if length(compName) ~= 2
        return;
    end

    compName = [upper(compName(1)), compName(2)];
end