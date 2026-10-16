%% ========================================================================
%  任意通道频率相关复数标定矩阵 H_q(f) + 动态六分量重构与验证
%  1.mat ~ 50.mat, Ui = FFT(U)/FFT(F), 不做重复试验平均
%  先缓存24通道原始频响，改变 selected_channels 后无需重复FFT
%  MATLAB R2019b+ (脚本末尾含局部函数)
%% ========================================================================
clear; clc;

%% ===== 0. 用户参数：通常只需要修改此区域 =====
calib_data_path = 'D:\file\code\matlab\app\bd_freqall_260516\data\bd_jizhenqi';
dynamic_data_path = calib_data_path;
c_file = 'c6_50_2605.mat';
c_variable = 'c6_50_2605';
output_dir = fullfile(pwd, 'H_selected_output');
if ~exist(output_dir,'dir'), mkdir(output_dir); end

number_position = 50;
cols_to_extract = 1:50;            % C 的列、作用点/文件编号，按相同顺序排列


% % % % 保留全部24通道
% % % selected_channels = 1:24;
% % % % 保留前18个通道
% % % selected_channels = 1:18;
% % % % 任意12个通道
% % % selected_channels = [1 3 5 7 9 11 13 15 17 19 21 23];

selected_channels = [10 11 12 22 23 24];  % <<< 任意保留通道编号
freq_list = 1:1:1000;
T_start = 3;
T_end = 13;
F_eps_ratio = 1e-6;
force_reextract = false;            % true: 强制重新提取所有原始频响
run_dynamic = true;
dynamic_filename = '34.mat';
val_position = 34;                  % 若用于验证，应是与 dynamic_filename 对应的加载点
make_plots = true;

%% ===== 1. 输入检查 =====
selected_channels = selected_channels(:).';
cols_to_extract = cols_to_extract(:).';
freq_list = freq_list(:).';
assert(~isempty(selected_channels) && all(selected_channels==fix(selected_channels)) && ...
    all(selected_channels>=1 & selected_channels<=24) && ...
    numel(unique(selected_channels))==numel(selected_channels), ...
    'selected_channels 必须为1~24内不重复的整数。');
assert(numel(selected_channels)>=6, '独立重构六分量通常至少需要6个通道。');
assert(~isempty(cols_to_extract) && all(cols_to_extract==fix(cols_to_extract)) && ...
    all(cols_to_extract>=1 & cols_to_extract<=number_position) && ...
    numel(unique(cols_to_extract))==numel(cols_to_extract), ...
    'cols_to_extract 必须为1~number_position内不重复的整数。');
assert(all(isfinite(freq_list)) && all(freq_list>0) && ...
    all(diff(freq_list)>0), 'freq_list 必须为严格递增的正频率。');
assert(T_start>=0 && T_end>T_start, '检查时间截取范围。');
q = numel(selected_channels);
nf = numel(freq_list);

s = load(c_file, c_variable);
assert(isfield(s,c_variable), '未找到 C 变量 %s',c_variable);
C_all = s.(c_variable);
assert(size(C_all,1)==6 && size(C_all,2)>=max(cols_to_extract), ...
    'C 必须是6行且列数覆盖所选加载点。');
C = C_all(:,cols_to_extract);

channel_tag = sprintf('%02d_',selected_channels);
channel_tag(end) = [];
cache_file = fullfile(output_dir,'Ui_all24_cache.mat');
H_file = fullfile(output_dir,sprintf('H_q%d_ch_%s.mat',q,channel_tag));

%% ===== 2. 获取原始24通道频响（缓存或重新提取） =====
% 缓存与所选通道无关：24 × number_position × nf
cache_ok = false;
if ~force_reextract && exist(cache_file,'file')
    meta = load(cache_file,'cache_meta');
    if isfield(meta,'cache_meta')
        old = meta.cache_meta;
        cache_ok = isequal(old.path,calib_data_path) && ...
            isequal(old.number_position,number_position) && ...
            isequal(old.freq_list,freq_list) && ...
            isequal(old.T_start,T_start) && isequal(old.T_end,T_end) && ...
            isequal(old.F_eps_ratio,F_eps_ratio);
    end
end
if cache_ok
    fprintf('使用已有24通道频响缓存：%s\n',cache_file);
    z = load(cache_file,'Ui_all','fs_record','N_record', ...
        'actual_freq_calib','freq_error_calib');
    Ui_all = z.Ui_all;
    fs_record = z.fs_record;
    N_record = z.N_record;
    actual_freq_calib = z.actual_freq_calib;
    freq_error_calib = z.freq_error_calib;
    assert(isequal(size(Ui_all),[24 number_position nf]), '缓存维度不匹配。');
else
    Ui_all = complex(nan(24,number_position,nf));
    fs_record = nan(number_position,1);
    N_record = nan(number_position,1);
    actual_freq_calib = nan(number_position,nf);
    freq_error_calib = nan(number_position,nf);

    for p = 1:number_position
        fp = fullfile(calib_data_path,sprintf('%d.mat',p));
        if ~isfile(fp), warning('缺少标定文件：%s',fp); continue; end
        d = load(fp);
        if ~isfield(d,'signal_in1')
            warning('%s 缺少 signal_in1 时间参考，跳过。',fp); continue;
        end
        t0 = d.signal_in1.Time(:);
        if numel(t0)<2, warning('%s 时间轴过短。',fp); continue; end
        fs = 1/mean(diff(t0));
        if ~isfinite(fs) || fs<=0, warning('%s 采样率非法。',fp); continue; end
        i1 = round(T_start*fs)+1;
        i2 = round(T_end*fs);
        if i1<1 || i2>numel(t0) || i2<i1
            warning('%s 时间截取范围不在数据内。',fp); continue;
        end
        N = i2-i1+1;
        fs_record(p) = fs;
        N_record(p) = N;
        valid_freq = freq_list<=fs/2;
        fgrid = (0:floor(N/2))*fs/N; % 只选非负频率 FFT bin
        bin = zeros(1,nf);
        for k=find(valid_freq)
            [~,bin(k)] = min(abs(fgrid-freq_list(k)));
        end
        actual_freq_calib(p,valid_freq) = fgrid(bin(valid_freq));
        freq_error_calib(p,valid_freq) = fgrid(bin(valid_freq))-freq_list(valid_freq);

        if isfield(d,'signal_Flc')
            tsF = d.signal_Flc;
        elseif isfield(d,sprintf('signal_Flc%d',p))
            tsF = d.(sprintf('signal_Flc%d',p));
        else
            warning('%s 缺少力参考通道，跳过。',fp); continue;
        end
        fdata = tsF.Data(:);
        if numel(fdata)<i2, warning('%s 力参考长度不足。',fp); continue; end
        fcut = fdata(i1:i2);
        fcut = fcut-mean(fcut);
        Ffft = fft(fcut);
        Fspec = complex(nan(1,nf));
        Fspec(valid_freq) = Ffft(bin(valid_freq));
        finiteF = abs(Fspec(isfinite(Fspec)));
        if isempty(finiteF) || max(finiteF)==0
            warning('%s 力频谱无效。',fp); continue;
        end
        epsF = max(finiteF)*F_eps_ratio;
        good = valid_freq & isfinite(Fspec) & abs(Fspec)>epsF;
        if ~any(good), warning('%s 无有效频点。',fp); continue; end

        for ch = 1:24
            fn = sprintf('signal_in%d',ch);
            if ~isfield(d,fn), warning('%s 缺少 %s。',fp,fn); continue; end
            x = d.(fn).Data(:);
            if numel(x)<i2, warning('%s 中 %s 数据长度不足。',fp,fn); continue; end
            x = x(i1:i2);
            x = x-mean(x);
            X = fft(x);
          U_selected = X(bin(good));
F_selected = Fspec(good);

Ui_all(ch,p,good) = reshape(U_selected(:).' ./ F_selected(:).', 1, 1, []);
        end
        fprintf('原始频响提取 %d/%d：%s\n',p,number_position,fp);
    end
    cache_meta = struct('path',calib_data_path,'number_position',number_position, ...
        'freq_list',freq_list,'T_start',T_start,'T_end',T_end, ...
        'F_eps_ratio',F_eps_ratio);
    save(cache_file,'Ui_all','fs_record','N_record','actual_freq_calib', ...
        'freq_error_calib','cache_meta','-v7.3');
    fprintf('24通道原始频响缓存已保存：%s\n',cache_file);
end
% 提示：若 1.mat~50.mat 内容被替换，手动设 force_reextract=true。

%% ===== 3. 对任意选定通道重新计算 H_q(f) =====
Ui_calib = Ui_all(selected_channels,cols_to_extract,:);
H_complex = complex(nan(6,q,nf));
cond_U = nan(1,nf);
rank_U = nan(1,nf);
valid_position_count = zeros(1,nf);
ui_energy = nan(1,nf);
for k = 1:nf
    U = Ui_calib(:,:,k);
    good_pos = all(isfinite(U),1) & all(isfinite(C),1);
    valid_position_count(k) = nnz(good_pos);
    if nnz(good_pos)<6
        continue;
    end
    U = U(:,good_pos);
    F = C(:,good_pos);
    r = rank(U);
    rank_U(k) = r;
    sigma = svd(U);
    if isempty(sigma) || min(sigma)==0
        cond_U(k) = inf;
    else
        cond_U(k) = max(sigma)/min(sigma);
    end
    ui_energy(k) = norm(U,'fro');
    if r<6
        % 不用低秩响应生成看似正常但不具备六维可辨识性的 H
        continue;
    end
    H_complex(:,:,k) = F*pinv(U);
end
valid_H = squeeze(all(all(isfinite(H_complex),1),2)).';
save(H_file,'H_complex','freq_list','selected_channels','cols_to_extract', ...
    'C','Ui_calib','cond_U','rank_U','valid_position_count','ui_energy', ...
    'fs_record','N_record','actual_freq_calib','freq_error_calib', ...
    'F_eps_ratio','T_start','T_end','-v7.3');
fprintf('\n标定完成：H=%d × %d × %d，有效频点 %d/%d\n',6,q,nf,nnz(valid_H),nf);
fprintf('输出：%s\n',H_file);

if make_plots
    figure('Name','标定质量诊断');
    tiledlayout(4,1,'TileSpacing','compact');
    nexttile; plot(freq_list,ui_energy); grid on; ylabel('||U||_F'); title('响应能量');
    nexttile; semilogy(freq_list,cond_U); grid on; ylabel('cond(U)'); title('条件数');
    nexttile; plot(freq_list,rank_U); grid on; ylabel('rank(U)'); title('矩阵秩');
    nexttile; plot(freq_list,valid_position_count); grid on;
    ylabel('加载点数'); xlabel('Frequency / Hz'); title('有效标定加载点');

    comp_idx = 1;
    plot_channel = selected_channels(1);
    local_ch = find(selected_channels==plot_channel,1);
    hc = squeeze(H_complex(comp_idx,local_ch,:));
    figure('Name','选定H元素的幅值与相位');
    tiledlayout(2,1);
    nexttile; plot(freq_list,abs(hc)); grid on;
    ylabel('|H|'); title(sprintf('分量%d - 原始通道%d',comp_idx,plot_channel));
    nexttile; plot(freq_list,rad2deg(angle(hc))); grid on;
    ylabel('Phase / deg'); xlabel('Frequency / Hz');

    hnorm = nan(1,nf);
    for k=1:nf
        if all(isfinite(H_complex(:,:,k)),'all')
            hnorm(k)=norm(H_complex(:,:,k),'fro');
        end
    end
    figure('Name','H的Frobenius范数'); plot(freq_list,hnorm);
    grid on; xlabel('Frequency / Hz'); ylabel('||H||_F');
end

%% ===== 4. 动态数据重构：只读取 selected_channels =====
if run_dynamic
    dynamic_file_fullpath = fullfile(dynamic_data_path,dynamic_filename);
    if ~isfile(dynamic_file_fullpath)
        warning('动态文件不存在：%s。仅输出标定矩阵。',dynamic_file_fullpath);
    else
        d = load(dynamic_file_fullpath);
        ref_name = sprintf('signal_in%d',selected_channels(1));
        assert(isfield(d,ref_name),'动态文件缺少 %s',ref_name);
        tt = d.(ref_name).Time(:);
        fs = 1/mean(diff(tt));
        assert(isfinite(fs) && fs>0,'动态采样率无效。');
        i1 = round(T_start*fs)+1;
        i2 = round(T_end*fs);
        assert(i1>=1 && i2<=numel(tt) && i2>=i1,'动态数据时间长度不足。');
        t = tt(i1:i2);
        N = numel(t);
        fgrid = (0:floor(N/2))*fs/N;
        nyquist_freq = fs/2;
        U_full = complex(zeros(q,N));
        for j=1:q
            ch = selected_channels(j);
            fn = sprintf('signal_in%d',ch);
            assert(isfield(d,fn),'动态文件缺少 %s',fn);
            x = d.(fn).Data(:);
            assert(numel(x)>=i2,'%s 动态数据长度不足',fn);
            x = x(i1:i2); x = x-mean(x);
            U_full(j,:) = fft(x).';
        end

        FM_full = complex(zeros(6,N));
        FM_freq_list = complex(nan(6,nf));
        actual_freq_list_dynamic = nan(1,nf);
        freq_error_dynamic = nan(1,nf);
        for k=1:nf
            f = freq_list(k);
            if f>nyquist_freq || ~valid_H(k), continue; end
            [~,idx_pos] = min(abs(fgrid-f));
            U_vec = U_full(:,idx_pos);
            FM_vec = H_complex(:,:,k)*U_vec;
            FM_full(:,idx_pos) = FM_vec;
            FM_freq_list(:,k) = FM_vec;
            actual_freq_list_dynamic(k) = fgrid(idx_pos);
            freq_error_dynamic(k) = fgrid(idx_pos)-f;
            idx_neg = N-idx_pos+2;
            if idx_pos>1 && idx_neg<=N && idx_neg~=idx_pos
                FM_full(:,idx_neg) = conj(FM_vec);
            end
        end
        fm_time = real(ifft(FM_full,[],2));
        Fx = fm_time(1,:).'; Fy = fm_time(2,:).'; Fz = fm_time(3,:).';
        Mx = fm_time(4,:).'; My = fm_time(5,:).'; Mz = fm_time(6,:).';
        FM_amp_freq_list = 2*abs(FM_freq_list)/N;
        FM_phase_freq_list_deg = rad2deg(angle(FM_freq_list));
        names = {'Fx','Fy','Fz','Mx','My','Mz'};
        if make_plots
            figure('Name','六分量频域重构'); tiledlayout(2,3);
            for j=1:6
                nexttile; plot(freq_list,FM_amp_freq_list(j,:)); grid on;
                title(names{j}); xlabel('Frequency / Hz'); ylabel('Amplitude');
            end
            figure('Name','六分量时域重构'); tiledlayout(2,3);
            for j=1:6
                nexttile; plot(t,fm_time(j,:)); grid on;
                title(names{j}); xlabel('Time / s'); ylabel(names{j});
            end
        end
        dynamic_result = struct('t',t,'fm_time',fm_time,'Fx',Fx,'Fy',Fy, ...
            'Fz',Fz,'Mx',Mx,'My',My,'Mz',Mz,'FM_full',FM_full, ...
            'FM_freq_list',FM_freq_list,'FM_amp_freq_list',FM_amp_freq_list, ...
            'FM_phase_freq_list_deg',FM_phase_freq_list_deg, ...
            'freq_list',freq_list,'selected_channels',selected_channels, ...
            'source_file',dynamic_filename, ...
            'actual_freq_list_dynamic',actual_freq_list_dynamic, ...
            'freq_error_dynamic',freq_error_dynamic,'T_start',T_start,'T_end',T_end);
        dynamic_output = fullfile(output_dir,sprintf('dynamic_q%d_ch_%s.mat',q,channel_tag));
        save(dynamic_output,'dynamic_result','-v7.3');
        fprintf('动态六分量结果已保存：%s\n',dynamic_output);

        %% ===== 5. 与原代码一致的幅值比回代检查 =====
        col_idx = find(cols_to_extract==val_position,1);
        if isempty(col_idx)
            warning('val_position不在标定加载点列表中，跳过幅值比检查。');
        elseif isfield(d,'signal_Flc') || isfield(d,sprintf('signal_Flc%d',val_position))
            if isfield(d,'signal_Flc')
                F_ts = d.signal_Flc;
            else
                F_ts = d.(sprintf('signal_Flc%d',val_position));
            end
            x = F_ts.Data(:);
            if numel(x)>=i2
                x = x(i1:i2); x = x-mean(x);
                XF = fft(x);
                F_ref_freq_list = complex(nan(6,nf));
                for k=1:nf
                    if freq_list(k)>nyquist_freq, continue; end
                    [~,idx_pos] = min(abs(fgrid-freq_list(k)));
                    F_ref_freq_list(:,k) = C(:,col_idx)*XF(idx_pos);
                end
                denom = abs(F_ref_freq_list);
                ratio_amp = abs(FM_freq_list)./denom;
                eps_ref = max(denom(isfinite(denom)),[],'omitnan')*1e-6;
                ratio_amp(denom<=eps_ref | ~isfinite(denom)) = nan;
                dynamic_result.ratio_amp = ratio_amp;
                save(dynamic_output,'dynamic_result','-v7.3');
                if make_plots
                    figure('Name','重构/参考幅值比（标定数据回代检查）');
                    tiledlayout(2,3);
                    for j=1:6
                        nexttile; plot(freq_list,ratio_amp(j,:)); hold on;
                        yline(1,'--'); grid on; ylim([0 2]);
                        xlabel('Frequency / Hz'); ylabel('|F_{est}|/|F_{ref}|'); title(names{j});
                    end
                end
            else
                warning('动态文件中的力参考长度不足，跳过幅值比检查。');
            end
        else
            warning('动态文件缺少力参考通道，跳过幅值比检查。');
        end
        fprintf('提示：如果动态文件属于1~50的标定数据，本图是回代检查而非独立验证。\n');
        fprintf('提示：频谱只在 freq_list 目标频点赋值，IFFT 是相应频带的重构。\n');
    end
end
