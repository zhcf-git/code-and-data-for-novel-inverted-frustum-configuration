%% ================================================================
%  手动选择单个验证文件
%  使用 H_complex(f) 验证实际测量结果
%
%  逻辑：
%       预测值 yc = H_k * U_k
%       参考值 yj = C(:, reference_position) * F_k
%       Normalized output ratio (-)
%  U_k 来自所选文件的 24 路电压信号 FFT
%  F_k 来自所选文件的力锤信号 FFT
%  C 来自 c6_50_2605.mat，尺寸为 6 × 50
%相比于7  这个是手动选文件和在c矩阵中对应的
%  ================================================================

clear;
clc;

%% -------------------- 基本参数 --------------------
dynamic_data_path = 'D:\file\code\matlab\app\bd_freqall_260516\data\bd_jizhenqi';
addpath(dynamic_data_path);

load('H_complex_jizhenqi_new.mat');     % H_complex, freq_list

%% -------------------- 加载 C 矩阵 --------------------
C_data = load('c6_50_2605.mat');
C_names = fieldnames(C_data);
C = C_data.(C_names{1});                % 默认读取 mat 文件中的第一个变量


%% ---------- 手动设置该文件对应 C 矩阵的第几列 ----------
reference_position = 18;     % 这里改成当前文件对应的点号，范围 1~50


if size(C, 1) ~= 6 || size(C, 2) ~= 50
    error('C 矩阵尺寸应为 6 × 50，当前尺寸为 %d × %d。', ...
        size(C, 1), size(C, 2));
end

freq_range = [1, 1000];

T_start = 3;
T_end   = 13;

num_position = 1;
num_channel = 24;

%% ---------- 手动选择验证文件 ----------
[file, path] = uigetfile('*.mat', '请选择一个验证文件');

if isequal(file, 0)
    error('未选择验证文件。');
end

selected_file = fullfile(path, file);





if reference_position < 1 || reference_position > size(C, 2)
    error('reference_position 必须在 1 到 %d 之间。', size(C, 2));
end

freq_idx_H_list = find(freq_list >= freq_range(1) & freq_list <= freq_range(2));
freq_used = freq_list(freq_idx_H_list);
num_freq = length(freq_used);

if num_freq == 0
    error('H_complex 中没有 %.3f Hz 到 %.3f Hz 范围内的频点。', ...
        freq_range(1), freq_range(2));
end

fprintf('使用 H_complex 中 %.3f Hz 到 %.3f Hz，共 %d 个频点。\n', ...
    freq_used(1), freq_used(end), num_freq);

fprintf('当前验证文件：%s\n', selected_file);
fprintf('对应 C 矩阵第 %d 列。\n', reference_position);

%% -------------------- 预分配 --------------------
U_complex = zeros(num_channel, num_position, num_freq, 'like', 1i);
F_complex = zeros(1, num_position, num_freq, 'like', 1i);

fs_record = nan(num_position, 1);
N_record  = nan(num_position, 1);

actual_freq_data = nan(num_position, num_freq);
freq_error_data  = nan(num_position, num_freq);

%% ================================================================
%  提取所选验证文件频率区间内的复数频谱
%  ================================================================

data = load(selected_file);

if ~isfield(data, 'signal_in1')
    error('所选文件中缺少变量：signal_in1');
end

ts_ref = data.signal_in1;
t_ref = ts_ref.Time(:);

dt = mean(diff(t_ref));
fs = 1 / dt;

idx_start = round(T_start * fs) + 1;
idx_end   = round(T_end * fs);

if idx_end > length(t_ref)
    error('所选文件长度不足 %.2f s。', T_end);
end

N = idx_end - idx_start + 1;

fs_record(1) = fs;
N_record(1)  = N;

f_full = (0:N-1) * fs / N;

freq_idx_data_list = zeros(num_freq, 1);

for k = 1:num_freq
    [~, freq_idx_data_list(k)] = min(abs(f_full - freq_used(k)));
end

actual_freq_data(1, :) = f_full(freq_idx_data_list);
freq_error_data(1, :) = f_full(freq_idx_data_list) - freq_used;

%% ---------- 24 路电压信号 ----------
wait_yz = waitbar(0, '正在提取24路电压信号频谱...');

for ch = 1:num_channel

    signal_name = sprintf('signal_in%d', ch);

    if ~isfield(data, signal_name)
        error('所选文件中缺少变量：%s', signal_name);
    end

    ts = data.(signal_name);
    x = ts.Data(:);

    if length(x) < idx_end
        error('所选文件中 %s 数据长度不足。', signal_name);
    end

    x_cut = x(idx_start:idx_end);
    x_cut = x_cut - mean(x_cut);

    X = fft(x_cut) / N;

    U_complex(ch, 1, :) = reshape(X(freq_idx_data_list), ...
        1, 1, num_freq);

    waitbar(ch / num_channel, wait_yz, ...
        sprintf('正在处理第 %d / %d 路电压...', ...
        ch, num_channel));
end

close(wait_yz);

%% ---------- 力锤信号 ----------
if isfield(data, 'signal_Flc')
    ts_lc = data.signal_Flc;
else
    force_name = sprintf('signal_Flc%d', reference_position);

    if ~isfield(data, force_name)
        error('所选文件中缺少力锤变量：signal_Flc 或 %s', force_name);
    end

    ts_lc = data.(force_name);
end

value_f = ts_lc.Data(:);

if length(value_f) < idx_end
    error('所选文件中力锤通道数据长度不足。');
end

value_f_cut = value_f(idx_start:idx_end);
value_f_cut = value_f_cut - mean(value_f_cut);

F_fft = fft(value_f_cut) / N;

F_complex(1, 1, :) = reshape(F_fft(freq_idx_data_list), ...
    1, 1, num_freq);

%% ================================================================
%  用 H_complex 计算预测六分量，并与 C(:,reference_position)*F 对比
%  ================================================================

FM_sj_amp = zeros(6, num_position, num_freq);
FM_yc_amp = zeros(6, num_position, num_freq);
FM_diff_amp = zeros(6, num_position, num_freq);
FM_rel_error = zeros(6, num_position, num_freq);

FM_sj_complex_all = zeros(6, num_position, num_freq, 'like', 1i);
FM_yc_complex_all = zeros(6, num_position, num_freq, 'like', 1i);

F_abs_all = abs(squeeze(F_complex(1, 1, :)));
F_eps_ratio = 1e-6;
F_eps = max(F_abs_all) * F_eps_ratio;

wait_calc = waitbar(0, '正在计算各频点六分量误差...');

for k = 1:num_freq

    H_k = H_complex(:, :, freq_idx_H_list(k));   % 6 × 24
    U_k = U_complex(:, 1, k);                    % 24 × 1
    F_k = F_complex(1, 1, k);                    % 1 × 1

    if abs(F_k) > F_eps

        %% ---------- 预测值 yc：实际测量形式 H * U ----------
        FM_yc_complex = H_k * U_k;               % 6 × 1

        %% ---------- 参考值 yj：C 矩阵对应列 * 输入力 ----------
        FM_sj_complex = C(:, reference_position) .* F_k;   % 6 × 1

    else

        FM_yc_complex = nan(6, 1);
        FM_sj_complex = nan(6, 1);

    end

    %% ---------- 保存复数结果 ----------
    FM_yc_complex_all(:, 1, k) = FM_yc_complex;
    FM_sj_complex_all(:, 1, k) = FM_sj_complex;

    %% ---------- 幅值误差 ----------
    FM_sj_amp(:, 1, k) = abs(FM_sj_complex);
    FM_yc_amp(:, 1, k) = abs(FM_yc_complex);

    FM_diff_amp(:, 1, k) = FM_sj_amp(:, 1, k) - FM_yc_amp(:, 1, k);

    temp_rel_error = abs(FM_diff_amp(:, 1, k)) ./ abs(FM_sj_amp(:, 1, k));
    temp_rel_error(abs(FM_sj_amp(:, 1, k)) < eps) = NaN;

    FM_rel_error(:, 1, k) = temp_rel_error;

    waitbar(k / num_freq, wait_calc, ...
        sprintf('正在计算第 %d / %d 个频点：%.3f Hz', ...
        k, num_freq, freq_used(k)));

end

close(wait_calc);





%% ================================================================
%  额外模块：生成该文件的六分量耦合干扰验证曲线
%  ================================================================

FM_yc_one = squeeze(FM_yc_amp(:, 1, :));   % 6 × num_freq
FM_sj_one = squeeze(FM_sj_amp(:, 1, :));   % 6 × num_freq

tol = 1e-12;

ratio_one = nan(size(FM_yc_one));          % 6 × num_freq

for k = 1:num_freq

    %% ===== 力部分：Fx Fy Fz =====
    force_ref_now = abs(FM_sj_one(1:3, k));
    force_est_now = abs(FM_yc_one(1:3, k));

    force_valid = find(force_ref_now > tol);

    if ~isempty(force_valid)

        main_force_idx = force_valid(1);

        main_force_ref = force_ref_now(main_force_idx);
        main_force_est = force_est_now(main_force_idx);

        for f = 1:3
            if f == main_force_idx
                ratio_one(f, k) = force_est_now(f) / main_force_ref;
            else
                if main_force_est > tol
                    ratio_one(f, k) = force_est_now(f) / main_force_est;
                else
                    ratio_one(f, k) = NaN;
                end
            end
        end

    end

    %% ===== 力矩部分：Mx My Mz =====
    moment_ref_now = abs(FM_sj_one(4:6, k));
    moment_est_now = abs(FM_yc_one(4:6, k));

    max_moment_est = max(moment_est_now);

    for m = 1:3

        row = 3 + m;

        if moment_ref_now(m) > tol
            ratio_one(row, k) = moment_est_now(m) / moment_ref_now(m);
        else
            if max_moment_est > tol
                ratio_one(row, k) = moment_est_now(m) / max_moment_est;
            else
                ratio_one(row, k) = NaN;
            end
        end

    end

end

%% ---------- 绘制六分量耦合干扰曲线 ----------
figure('Name', sprintf('%s 耦合干扰曲线', file), ...
    'NumberTitle', 'off', ...
    'Position', [150, 100, 850, 600]);

tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

comp_names = {'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};

for comp = 1:6

    nexttile;

    plot(freq_used, ratio_one(comp, :), 'LineWidth', 1.3);
    hold on;
    yline(1, '--', 'LineWidth', 1.2);
    grid on;

    xlabel('Frequency / Hz');
    ylabel('耦合干扰比');
    title(comp_names{comp});

    ylim([0 2]);
    xlim([min(freq_used), max(freq_used)]);

end

sgtitle(sprintf('%s 对应 C 矩阵第 %d 列的六分量耦合干扰曲线', ...
    file, reference_position), ...
    'Interpreter', 'none');


% % % 
% % % %% ================================================================
% % % %  额外模块：分别绘制六分量耦合干扰曲线
% % % %  每个分量单独一个 figure
% % % %  ================================================================
% % % 
% % % comp_names = {'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};
% % % comp_latex = {'$F_x$', '$F_y$', '$F_z$', '$M_x$', '$M_y$', '$M_z$'};
% % % 
% % % for comp = 1:6
% % % 
% % %     figure('Name', sprintf('%s Coupling Interference - %s', file, comp_names{comp}), ...
% % %         'NumberTitle', 'off', ...
% % %         'Position', [200 + 40 * comp, 120 + 30 * comp, 750, 480]);
% % % 
% % %     plot(freq_used, ratio_one(comp, :), ...
% % %         'Color', [0.8500 0.3250 0.0980], ...
% % %         'LineWidth', 2);
% % %     hold on;
% % % 
% % %     yline(1, '--', ...
% % %         'Color', [0 0 0], ...
% % %         'LineWidth', 1.5);
% % %     hold on;
% % % 
% % %     grid on;
% % %     box on;
% % % 
% % %     xlabel('Frequency (Hz)', 'Interpreter', 'latex');
% % % 
% % %     ylabel({['Response - Component ', comp_latex{comp}]; ...
% % %         'Coupling interference ratio'}, ...
% % %         'Interpreter', 'latex');
% % % 
% % %     LegendString = [comp_latex(comp); {'Reference = 1'}];
% % % 
% % %     L1 = legend(LegendString, 'Location', 'best');
% % %     set(L1, 'Interpreter', 'latex');
% % % 
% % % 
% % %     set(gca, 'TickLabelInterpreter', 'latex');
% % % 
% % %     xlim([min(freq_used), max(freq_used)]);
% % %     ylim([0 2]);
% % % 
% % % end



%% ================================================================
%  额外模块：力分量 Fx Fy Fz 同时显示
%  ================================================================

figure('Name', sprintf('%s Coupling Interference - Force Components', file), ...
       'NumberTitle', 'off', ...
       'Position', [260, 180, 780, 500]);

plot(freq_used, ratio_one(1, :), ...
     'Color', [0.8500 0.3250 0.0980], ...
     'LineWidth', 2);
hold on;

plot(freq_used, ratio_one(2, :), ...
     'Color', [0 0.4470 0.7410], ...
     'LineWidth', 2);
hold on;

plot(freq_used, ratio_one(3, :), ...
     'Color', [0.4660 0.6740 0.1880], ...
     'LineWidth', 2);
hold on;

yline(1, '--', ...
      'Color', [0 0 0], ...
      'LineWidth', 1.5);
hold on;

grid on;
box on;

xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel({'Response - Force Components'; ...
        'Coupling interference ratio'}, ...
        'Interpreter', 'latex');

LegendString = {'$F_x$', '$F_y$', '$F_z$', 'Reference = 1'};

L1 = legend(LegendString, 'Location', 'best');
set(L1, 'Interpreter', 'latex');

set(gca, 'TickLabelInterpreter', 'latex');

xlim([min(freq_used), max(freq_used)]);
ylim([0 2]);


%% ================================================================
%  额外模块：力矩分量 Mx My Mz 同时显示
%  ================================================================

figure('Name', sprintf('%s Coupling Interference - Moment Components', file), ...
       'NumberTitle', 'off', ...
       'Position', [320, 220, 780, 500]);

plot(freq_used, ratio_one(4, :), ...
     'Color', [0.8500 0.3250 0.0980], ...
     'LineWidth', 2);
hold on;

plot(freq_used, ratio_one(5, :), ...
     'Color', [0 0.4470 0.7410], ...
     'LineWidth', 2);
hold on;

plot(freq_used, ratio_one(6, :), ...
     'Color', [0.4660 0.6740 0.1880], ...
     'LineWidth', 2);
hold on;

yline(1, '--', ...
      'Color', [0 0 0], ...
      'LineWidth', 1.5);
hold on;

grid on;
box on;

xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel({'Response - Moment Components'; ...
        'Coupling interference ratio'}, ...
        'Interpreter', 'latex');

LegendString = {'$M_x$', '$M_y$', '$M_z$', 'Reference = 1'};

L1 = legend(LegendString, 'Location', 'best');
set(L1, 'Interpreter', 'latex');

set(gca, 'TickLabelInterpreter', 'latex');

xlim([min(freq_used), max(freq_used)]);
ylim([0 2]);