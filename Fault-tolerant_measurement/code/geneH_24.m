%% ================================================================
%  复数频率相关标定矩阵 H_complex(f)
%  使用 Ui/Fi 方法，单次数据：1.mat ~ 50.mat，  相较于力锤0520不做平均
%  截取时间段：3s ~ 13s
%  ================================================================

clear;
clc;

%% -------------------- 路径设置 --------------------
calib_data_path = 'D:\file\code\matlab\app\bd_freqall_260516\data\bd_jizhenqi';
addpath(calib_data_path);

%% -------------------- 基本参数 --------------------
number_position = 50;

freq_start = 1;
freq_end   = 1000;
freq_step  = 1;
freq_list  = freq_start:freq_step:freq_end;
num_freq   = length(freq_list);

%% -------------------- 截取时间设置 --------------------
T_start = 3;
T_end   = 13;

%% -------------------- 加载 C 矩阵 --------------------
load('c6_50_2605.mat');

cols_to_extract = 1:50;
C = c6_50_2605(:, cols_to_extract);

%% -------------------- 有效通道选择 --------------------
rows_to_extract = 1:24;

num_channel = length(rows_to_extract);

%% -------------------- 预分配 --------------------
Ui_all = nan(24, number_position, num_freq, 'like', 1i);

fs_record = nan(number_position, 1);
N_record  = nan(number_position, 1);

actual_freq_calib = nan(number_position, num_freq);
freq_error_calib  = nan(number_position, num_freq);

F_eps_ratio = 1e-6;

%% -------------------- 提取每个作用点的 Ui = U/F --------------------
h = waitbar(0, '正在提取 1.mat ~ 50.mat 的 Ui = U/F 频响...');

for positionnumber = 1:number_position

    calib_filename = sprintf('%d.mat', positionnumber);
    calib_file_fullpath = fullfile(calib_data_path, calib_filename);

    if ~exist(calib_file_fullpath, 'file')
        warning('文件不存在，跳过：%s', calib_file_fullpath);
        continue;
    end

    data = load(calib_file_fullpath);

    %% ---------- 用第 1 通道确定时间轴和采样率 ----------
    if ~isfield(data, 'signal_in1')
        warning('文件 %s 缺少 signal_in1，跳过。', calib_filename);
        continue;
    end

    ts_ref = data.signal_in1;
    t_ref = ts_ref.Time(:);

    dt = mean(diff(t_ref));
    fs = 1 / dt;

    idx_start = round(T_start * fs) + 1;
    idx_end   = round(T_end * fs);

    if idx_end > length(t_ref)
        warning('文件 %s 长度不足 %.2f s，跳过。', calib_filename, T_end);
        continue;
    end

    N = idx_end - idx_start + 1;

    fs_record(positionnumber) = fs;
    N_record(positionnumber)  = N;

    f_full = (0:N-1) * fs / N;

    if freq_end > fs/2
        warning('文件 %s：freq_end = %.2f Hz 超过 Nyquist = %.2f Hz。', ...
                calib_filename, freq_end, fs/2);
    end

    %% ---------- freq_list 对应 FFT bin ----------
    freq_idx = zeros(1, num_freq);

    for k = 1:num_freq
        [~, freq_idx(k)] = min(abs(f_full - freq_list(k)));
    end

    actual_freq_calib(positionnumber, :) = f_full(freq_idx);
    freq_error_calib(positionnumber, :) = f_full(freq_idx) - freq_list;

    %% ---------- 读取 24 个电压通道 FFT ----------
    U_fft_all_ch = zeros(24, num_freq, 'like', 1i);

    for ch = 1:24

        signal_name = sprintf('signal_in%d', ch);

        if ~isfield(data, signal_name)
            warning('文件 %s 缺少变量 %s，该通道记为 NaN。', ...
                    calib_filename, signal_name);
            U_fft_all_ch(ch, :) = nan;
            continue;
        end

        ts = data.(signal_name);
        x = ts.Data(:);

        if length(x) < idx_end
            warning('文件 %s 中 %s 数据长度不足，该通道记为 NaN。', ...
                    calib_filename, signal_name);
            U_fft_all_ch(ch, :) = nan;
            continue;
        end

        x_cut = x(idx_start:idx_end);
        x_cut = x_cut - mean(x_cut);

        X = fft(x_cut);

        U_fft_all_ch(ch, :) = X(freq_idx);

    end

    %% ---------- 读取力锤通道 FFT ----------
    if isfield(data, 'signal_Flc')
        ts_lc = data.signal_Flc;
    else
        possible_name = sprintf('signal_Flc%d', positionnumber);

        if isfield(data, possible_name)
            ts_lc = data.(possible_name);
        else
            warning('文件 %s 未找到 signal_Flc 或 %s，跳过。', ...
                    calib_filename, possible_name);
            continue;
        end
    end

    value_f = ts_lc.Data(:);

    if length(value_f) < idx_end
        warning('文件 %s 力锤通道长度不足，跳过。', calib_filename);
        continue;
    end

    value_f_cut = value_f(idx_start:idx_end);
    value_f_cut = value_f_cut - mean(value_f_cut);

    F_fft = fft(value_f_cut);
    F_k_all = F_fft(freq_idx);

    %% ---------- 计算 Ui = U/F ----------
    F_eps = max(abs(F_k_all)) * F_eps_ratio;

    if F_eps == 0 || isnan(F_eps)
        warning('文件 %s 力锤频谱异常，跳过。', calib_filename);
        continue;
    end

    for k = 1:num_freq

        F_k = F_k_all(k);

        if abs(F_k) <= F_eps
            Ui_all(:, positionnumber, k) = nan;
        else
            Ui_all(:, positionnumber, k) = U_fft_all_ch(:, k) ./ F_k;
        end

    end

    waitbar(positionnumber / number_position, h, ...
        sprintf('正在处理：%s  (%d / %d)', ...
        calib_filename, positionnumber, number_position));

end

close(h);

%% -------------------- 不做平均，直接使用 Ui_all --------------------
Ui_calib = Ui_all(rows_to_extract, cols_to_extract, :);

%% -------------------- 逐频点计算 H_complex(f) --------------------
H_complex = nan(6, num_channel, num_freq, 'like', 1i);

cond_U = nan(1, num_freq);
rank_U = nan(1, num_freq);
valid_position_count = zeros(1, num_freq);
ui_energy = nan(1, num_freq);

h = waitbar(0, '正在用 Ui/Fi 计算 H_complex(f)...');

for k = 1:num_freq

    Ui_k = Ui_calib(:, :, k);
    Fi_k = C;

    valid_pos = all(isfinite(Ui_k), 1);

    valid_position_count(k) = sum(valid_pos);

    if sum(valid_pos) < 6
        warning('%.1f Hz 有效作用点少于 6 个，跳过该频点。', freq_list(k));
        continue;
    end

    Ui_valid = Ui_k(:, valid_pos);
    Fi_valid = Fi_k(:, valid_pos);

    H_complex(:, :, k) = Fi_valid * pinv(Ui_valid);

    rank_U(k) = rank(Ui_valid);

    s = svd(Ui_valid);

    if isempty(s) || min(s) == 0
        cond_U(k) = inf;
    else
        cond_U(k) = max(s) / min(s);
    end

    ui_energy(k) = norm(Ui_valid, 'fro');

    waitbar(k / num_freq, h, ...
        sprintf('正在计算 %.1f Hz 的 H_complex(f)...', freq_list(k)));

end

close(h);

%% -------------------- 保存 H_complex(f) --------------------
save('H_complex.mat', ...
    'H_complex', ...
    'freq_list', ...
    'rows_to_extract', ...
    'cols_to_extract', ...
    'C', ...
    'Ui_all', ...
    'Ui_calib', ...
    'fs_record', ...
    'N_record', ...
    'actual_freq_calib', ...
    'freq_error_calib', ...
    'cond_U', ...
    'rank_U', ...
    'valid_position_count', ...
    'ui_energy', ...
    'F_eps_ratio', ...
    'T_start', ...
    'T_end', ...
    '-v7.3');

fprintf('\nUi/Fi 单次数据 1.mat ~ 50.mat 复数频率相关标定矩阵计算完成！\n');
fprintf('H_complex 维度：%d × %d × %d\n', ...
        size(H_complex,1), size(H_complex,2), size(H_complex,3));
fprintf('对应含义：6分量 × %d个电压通道 × %d个频率点\n', ...
        num_channel, num_freq);
fprintf('频率范围：%.1f Hz 到 %.1f Hz\n', freq_list(1), freq_list(end));
fprintf('截取时间段：%.2f s 到 %.2f s\n', T_start, T_end);

fprintf('\n标定阶段频率匹配误差：\n');
fprintf('最大绝对误差 = %.6f Hz\n', ...
        max(abs(freq_error_calib(:)), [], 'omitnan'));

%% -------------------- 标定质量诊断窗口 --------------------
figure('Name', 'Ui/Fi 标定矩阵质量诊断', ...
       'NumberTitle', 'off', ...
       'Position', [50, 50, 800, 600]);

tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(freq_list, ui_energy, 'LineWidth', 1.2);
grid on;
xlabel('目标频率 / Hz');
ylabel('||Ui||_F');
title('Ui 频响矩阵能量');

nexttile;
semilogy(freq_list, cond_U, 'LineWidth', 1.2);
grid on;
xlabel('目标频率 / Hz');
ylabel('cond(Ui)');
title('Ui_k 条件数');

nexttile;
plot(freq_list, rank_U, 'LineWidth', 1.2);
grid on;
xlabel('目标频率 / Hz');
ylabel('rank(Ui)');
title('Ui_k 秩');

nexttile;
plot(freq_list, valid_position_count, 'LineWidth', 1.2);
grid on;
xlabel('目标频率 / Hz');
ylabel('有效作用点数量');
title('每个频点参与 Ui/Fi 标定的有效作用点数量');

%% -------------------- 绘制 H_complex 某一条频率曲线 --------------------
comp_idx = 1;
ch_idx   = 3;

H_curve = squeeze(H_complex(comp_idx, ch_idx, :));
H_amp = abs(H_curve);

figure;
plot(freq_list, H_amp, 'LineWidth', 1.5);
grid on;
xlabel('Frequency / Hz');
ylabel('|H(f)|');
title(sprintf('H\\_{complex} 分量 %d 电压通道 %d 幅值曲线', comp_idx, ch_idx));
xlim([min(freq_list), max(freq_list)]);

H_phase_deg = rad2deg(angle(H_curve));

figure;
plot(freq_list, H_phase_deg, 'LineWidth', 1.5);
grid on;
xlabel('Frequency / Hz');
ylabel('Phase / deg');
title(sprintf('H\\_{complex} 分量 %d 电压通道 %d 相位曲线', comp_idx, ch_idx));
xlim([min(freq_list), max(freq_list)]);

%% -------------------- 绘制 H_complex Frobenius 范数 --------------------
W_fro_norm = nan(1, num_freq);

for k = 1:num_freq
    if all(isfinite(H_complex(:,:,k)), 'all')
        W_fro_norm(k) = norm(H_complex(:,:,k), 'fro');
    end
end

figure;
plot(freq_list, W_fro_norm, 'LineWidth', 1.5);
grid on;
xlabel('Frequency / Hz');
ylabel('||H(f)||_F');
title('标定矩阵 H(f) 的 Frobenius 范数随频率变化曲线');
xlim([min(freq_list), max(freq_list)]);

%% ================================================================
%  PART 2：用 H_complex(f) 处理真实动态 24 路电压信号
%  ================================================================

dynamic_data_path = 'D:\file\code\matlab\app\bd_freqall_260516\data\bd_jizhenqi';
dynamic_filename = '34.mat';
val_position = 34;


dynamic_file_fullpath = fullfile(dynamic_data_path, dynamic_filename);

if exist(dynamic_file_fullpath, 'file')

    addpath(dynamic_data_path);

    data_dyn = load(dynamic_file_fullpath);

    %% ---------- 用第一个有效通道确定时间轴 ----------
    first_ch = rows_to_extract(1);
    first_signal_name = sprintf('signal_in%d', first_ch);

    if ~isfield(data_dyn, first_signal_name)
        error('真实动态数据中缺少变量：%s', first_signal_name);
    end

    ts_ref = data_dyn.(first_signal_name);
    t_full = ts_ref.Time(:);

    dt = mean(diff(t_full));
    fs = 1 / dt;

    idx_start = round(T_start * fs) + 1;
    idx_end   = round(T_end * fs);

    if idx_end > length(t_full)
        error('真实动态数据长度不足 %.2f s。', T_end);
    end

    t = t_full(idx_start:idx_end);
    N = length(t);

    f_full = (0:N-1) * fs / N;
    nyquist_freq = fs / 2;

    if max(freq_list) > nyquist_freq
        warning('H_complex 的最高频率 %.2f Hz 超过当前动态信号 Nyquist 频率 %.2f Hz，将自动忽略超出部分。', ...
                max(freq_list), nyquist_freq);
    end

    %% ---------- 读取有效通道电压并 FFT ----------
    U_full = zeros(num_channel, N, 'like', 1i);

    for ch_idx = 1:num_channel

        ch = rows_to_extract(ch_idx);
        signal_name = sprintf('signal_in%d', ch);

        if ~isfield(data_dyn, signal_name)
            error('真实动态数据中缺少变量：%s', signal_name);
        end

        ts = data_dyn.(signal_name);
        x_full = ts.Data(:);

        if length(x_full) < idx_end
            error('真实动态数据中 %s 的长度不足 %.2f s。', signal_name, T_end);
        end

        x = x_full(idx_start:idx_end);
        x = x - mean(x);

        U_full(ch_idx, :) = fft(x).';

    end

    %% ---------- 根据 H_complex(f) 预测六分量完整频谱 ----------
    FM_full = zeros(6, N, 'like', 1i);

    FM_freq_list = zeros(6, num_freq, 'like', 1i);
    actual_freq_list_dynamic = nan(1, num_freq);
    freq_error_dynamic = nan(1, num_freq);

    for k = 1:num_freq

        target_f = freq_list(k);

        if target_f <= 0 || target_f > nyquist_freq
            continue;
        end

        if any(~isfinite(H_complex(:,:,k)), 'all')
            continue;
        end

        [~, idx_pos] = min(abs(f_full - target_f));

        U_vec = U_full(:, idx_pos);
        FM_vec = H_complex(:, :, k) * U_vec;

        FM_full(:, idx_pos) = FM_vec;

        FM_freq_list(:, k) = FM_vec;
        actual_freq_list_dynamic(k) = f_full(idx_pos);
        freq_error_dynamic(k) = actual_freq_list_dynamic(k) - target_f;

        if idx_pos > 1
            idx_neg = N - idx_pos + 2;

            if idx_neg <= N && idx_neg ~= idx_pos
                FM_full(:, idx_neg) = conj(FM_vec);
            end
        end

    end

    fprintf('\n动态阶段频率匹配误差：\n');
    fprintf('最大绝对误差 = %.6f Hz\n', ...
            max(abs(freq_error_dynamic), [], 'omitnan'));

    %% ---------- IFFT 得到六分量时域信号 ----------
    fm_time = real(ifft(FM_full, [], 2));

    Fx = fm_time(1, :).';
    Fy = fm_time(2, :).';
    Fz = fm_time(3, :).';
    Mx = fm_time(4, :).';
    My = fm_time(5, :).';
    Mz = fm_time(6, :).';

    %% ---------- 频域幅值和相位 ----------
    FM_amp_freq_list = abs(FM_freq_list) / N * 2;
    FM_phase_freq_list_deg = rad2deg(angle(FM_freq_list));

    component_names = {'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};

    plot_colors = {
        [1, 0.5, 0], ...
        'b', ...
        'g', ...
        [0.6, 0.3, 0.1], ...
        'm', ...
        [0, 0.75, 0.75]
    };

    %% ---------- 窗口 1：六分量频域幅值结果 ----------
    fig_freq = figure('Name', '六分量频域校正结果', ...
                      'NumberTitle', 'off', ...
                      'Position', [100, 100, 800, 600]);

    main_panel = uipanel('Parent', fig_freq, ...
                         'Position', [0, 0.1, 1, 0.9]);

    control_panel = uipanel('Parent', fig_freq, ...
                            'Position', [0, 0, 1, 0.1], ...
                            'BackgroundColor', [0.94 0.94 0.94]);

    f_plot = freq_list;
    amp_plot = FM_amp_freq_list;

    ax_freq = gobjects(1, 6);

    for comp = 1:6

        ax_freq(comp) = subplot(2, 3, comp, 'Parent', main_panel);

        plot(f_plot, amp_plot(comp, :), '.-', ...
             'Color', plot_colors{comp}, ...
             'LineWidth', 1.2);

        grid on;
        title(['频域幅值: ', component_names{comp}]);
        xlabel('频率 / Hz');
        ylabel('|X(f)|');
        xlim([1 800]);

    end

    uicontrol(control_panel, ...
              'Style', 'text', ...
              'String', 'X轴范围:', ...
              'Position', [210, 20, 70, 25], ...
              'BackgroundColor', [0.94 0.94 0.94]);

    start_freq_edit = uicontrol(control_panel, ...
                                'Style', 'edit', ...
                                'String', '1', ...
                                'Position', [285, 20, 60, 25]);

    end_freq_edit = uicontrol(control_panel, ...
                              'Style', 'edit', ...
                              'String', '600', ...
                              'Position', [355, 20, 60, 25]);

    uicontrol(control_panel, ...
              'Style', 'pushbutton', ...
              'String', '更新范围', ...
              'Position', [440, 20, 90, 25], ...
              'Callback', @(src, event) update_xlim_freq(ax_freq, ...
                                                         start_freq_edit, ...
                                                         end_freq_edit));

    %% ---------- 窗口 2：六分量时域结果 ----------
    fig_time = figure('Name', '六分量时域重构结果', ...
                      'NumberTitle', 'off', ...
                      'Position', [50, 60, 800, 600]);

    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    fm_cell = {Fx, Fy, Fz, Mx, My, Mz};

    for comp = 1:6

        nexttile;

        plot(t, fm_cell{comp}, ...
             'Color', plot_colors{comp}, ...
             'LineWidth', 1.2);

        grid on;
        title(['时域: ', component_names{comp}]);
        xlabel('Time / s');
        ylabel(component_names{comp});

    end

    sgtitle(sprintf('重构后的六分量信号, %.2f s 到 %.2f s', T_start, T_end));

    %% ---------- 保存真实动态六分量结果 ----------
    dynamic_result.t = t;
    dynamic_result.T_start = T_start;
    dynamic_result.T_end = T_end;

    dynamic_result.fm_time = fm_time;

    dynamic_result.Fx = Fx;
    dynamic_result.Fy = Fy;
    dynamic_result.Fz = Fz;
    dynamic_result.Mx = Mx;
    dynamic_result.My = My;
    dynamic_result.Mz = Mz;

    dynamic_result.FM_full = FM_full;
    dynamic_result.FM_freq_list = FM_freq_list;
    dynamic_result.FM_amp_freq_list = FM_amp_freq_list;
    dynamic_result.FM_phase_freq_list_deg = FM_phase_freq_list_deg;

    dynamic_result.freq_list = freq_list;
    dynamic_result.actual_freq_list_dynamic = actual_freq_list_dynamic;
    dynamic_result.freq_error_dynamic = freq_error_dynamic;

    dynamic_result.rows_to_extract = rows_to_extract;
    dynamic_result.source_file = dynamic_filename;

    % save('dynamic_6axis_result_UiFi_no_average_1to50.mat', ...
    %      'dynamic_result', ...
    %      '-v7.3');

    fprintf('\n真实动态信号六分量重构完成！\n');
    fprintf('输出文件：dynamic_6axis_result_UiFi_no_average_1to50.mat\n');

else

    fprintf('\n未找到真实动态数据文件：\n%s\n', dynamic_file_fullpath);
    fprintf('已完成 H_complex(f) Ui/Fi 不平均标定矩阵生成。\n');

end

%% ================================================================
%  PART 3：提取 87Hz 复数标定矩阵
%  ================================================================

% % [~, idx_87] = min(abs(freq_list - 87));
% % H_87_complex = H_complex(:, :, idx_87);
% % 
% % save('H_87_complex_from_Hfreq_UiFi_no_average_1to50.mat', ...
% %      'H_87_complex');
% % 
% % fprintf('\n已额外保存 H_87_complex_from_Hfreq_UiFi_no_average_1to50.mat。\n');

%% ================================================================
%  PART 2.5：验证曲线 |F_est| / |F_ref|
%  ================================================================

if exist('data_dyn', 'var')

    col_idx = find(cols_to_extract == val_position, 1);

    if isempty(col_idx)
        error('val_position = %d 不在 cols_to_extract 中，无法用 C 做参考。', val_position);
    end

    C_val = C(:, col_idx);

    %% ---------- 读取真实数据中的力锤通道 ----------
    if isfield(data_dyn, 'signal_Flc')
        ts_lc = data_dyn.signal_Flc;
    else
        possible_name = sprintf('signal_Flc%d', val_position);

        if isfield(data_dyn, possible_name)
            ts_lc = data_dyn.(possible_name);
        else
            error('真实数据中未找到 signal_Flc 或 %s。', possible_name);
        end
    end

    f_hammer_full = ts_lc.Data(:);

    if length(f_hammer_full) < idx_end
        error('力锤通道长度不足 %.2f s。', T_end);
    end

    f_hammer = f_hammer_full(idx_start:idx_end);
    f_hammer = f_hammer - mean(f_hammer);

    F_hammer_fft = fft(f_hammer);

    %% ---------- 提取 freq_list 对应的力锤频谱 ----------
    F_ref_freq_list = zeros(6, num_freq, 'like', 1i);

    for k = 1:num_freq

        target_f = freq_list(k);

        if target_f <= 0 || target_f > nyquist_freq
            F_ref_freq_list(:, k) = nan;
            continue;
        end

        [~, idx_pos] = min(abs(f_full - target_f));

        F_ref_freq_list(:, k) = C_val * F_hammer_fft(idx_pos);

    end

    %% ---------- 计算验证比值曲线 ----------
    ratio_amp = abs(FM_freq_list) ./ abs(F_ref_freq_list);

    F_ref_abs = abs(F_ref_freq_list);
    F_ref_eps = max(F_ref_abs(:), [], 'omitnan') * 1e-6;
    ratio_amp(F_ref_abs < F_ref_eps) = nan;

    %% ---------- 绘制六分量验证曲线 ----------
    figure('Name', 'H_complex 验证曲线 |F_重构| / |F_实际|', ...
           'NumberTitle', 'off', ...
           'Position', [150, 100, 800, 600]);

    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for comp = 1:6

        nexttile;

        plot(freq_list, ratio_amp(comp, :), 'LineWidth', 1.3);
        hold on;
        yline(1, '--', 'LineWidth', 1.2);
        grid on;

        xlabel('Frequency / Hz');
        ylabel('|F_{est}| / |F_{ref}|');
        title(['验证曲线: ', component_names{comp}]);

        ylim([0 2]);
        xlim([min(freq_list), max(freq_list)]);

    end

    sgtitle(sprintf('验证数据 %s, 作用点 %d, %.2f s 到 %.2f s, H\\_complex 重构幅值比', ...
            dynamic_filename, val_position, T_start, T_end));

end

%% ================================================================
%  本地函数：更新六个子图的 X 轴范围
%  ================================================================

function update_xlim_freq(ax_array, start_freq_edit, end_freq_edit)

    f_start = str2double(start_freq_edit.String);
    f_end   = str2double(end_freq_edit.String);

    if isnan(f_start) || isnan(f_end)
        errordlg('请输入数字形式的频率范围。', '频率范围错误');
        return;
    end

    if f_start <= 0
        errordlg('起始频率必须大于 0。', '频率范围错误');
        return;
    end

    if f_end <= f_start
        errordlg('终止频率必须大于起始频率。', '频率范围错误');
        return;
    end

    for i = 1:length(ax_array)
        if isgraphics(ax_array(i))
            xlim(ax_array(i), [f_start, f_end]);
        end
    end

end