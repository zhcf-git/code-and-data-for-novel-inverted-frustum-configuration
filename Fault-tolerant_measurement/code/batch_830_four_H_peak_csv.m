%% ========================================================================
% Batch 830 Hz: 30 files x four H models -> CSV
%
% Based on the original App algorithm:
%   - signal_in1 is the shared reference time axis
%   - no last-sample deletion, no fixed 3-13 s truncation
%   - clean/sort/unique channel time, interpolate if needed
%   - demean then fft(x,N), find nearest calibrated H frequency
%   - six-component SINGLE-SIDED amplitude: 2*abs(H*U)/N
%
% Peak selection: for each file, use the declared loading component from
% the FILENAME in the H24 reconstruction to find the largest peak within
% peak_band. At that SAME physical FFT bin, sample all 6 components for
% H6/H12/H18/H24. Only a common frequency makes comparisons fair.
%
% CSV: exactly one row per source file (30 expected), all 24 amplitudes.
% ========================================================================
clear; clc;

%% ------------------------- USER SETTINGS -------------------------------
data_dir = 'D:\file\code\matlab\app\bd_freqall_260516\data\liner\830';
H_dir = pwd;   % Four H files (change if stored elsewhere)
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
output_csv = fullfile(script_dir,'peak_830_30files_fourH_with_Flc.csv');
peak_band = [800 860];    % Hz; adjust only if the actual 830 Hz peak is outside
strict_30_files = true;   % true: stop if expected files missing
save_mat_details = false; % optional: additionally save per-file complex results

H_files = { ...
    'H_q6_ch_10_11_12_22_23_24.mat', ...
    'H_q12_ch_04_05_06_10_11_12_16_17_18_22_23_24.mat', ...
    'H_q18_ch_04_05_06_07_08_09_10_11_12_16_17_18_19_20_21_22_23_24.mat', ...
    'H_q24_ch_01_02_03_04_05_06_07_08_09_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24.mat' ...
};
expected_channels = { ...
    [10 11 12 22 23 24], ...
    [4 5 6 10 11 12 16 17 18 22 23 24], ...
    [4 5 6 7 8 9 10 11 12 16 17 18 19 20 21 22 23 24], ...
    1:24 ...
};
q_values = [6 12 18 24];
component_names = {'Fx','Fy','Fz','Mx','My','Mz'};
file_components = {'fx','fy','fz','mx','my','mz'};

assert(exist(data_dir,'dir')==7,'Data folder does not exist: %s',data_dir);
assert(all(isfinite(peak_band)) && peak_band(1)>0 && peak_band(2)>peak_band(1), ...
    'Invalid peak_band');

%% ---------------------- Load and check H -------------------------------
models = struct('q',{},'channels',{},'H',{},'freq',{});
for m = 1:4
    fn = fullfile(H_dir,H_files{m});
    assert(exist(fn,'file')==2,'H file not found: %s',fn);
    S = load(fn);
    assert(isfield(S,'H_complex') && isfield(S,'freq_list'), ...
        'H file requires H_complex and freq_list: %s',H_files{m});
    if isfield(S,'selected_channels')
        channels = double(S.selected_channels(:).');
    elseif isfield(S,'rows_to_extract')
        channels = double(S.rows_to_extract(:).');
    else
        error('Missing physical-channel order metadata: %s',H_files{m});
    end
    assert(isequal(channels,expected_channels{m}), ...
        'Channel order mismatch in H file: %s',H_files{m});
    f = double(S.freq_list(:).');
    H = S.H_complex;
    assert(size(H,1)==6 && size(H,2)==q_values(m) && size(H,3)==numel(f), ...
        'H size mismatch: %s',H_files{m});
    assert(all(diff(f)>0),'H frequency list must be increasing: %s',H_files{m});
    models(m).q = q_values(m);
    models(m).channels = channels;
    models(m).H = H;
    models(m).freq = f;
end

%% ----------------------- Prepare 30 filenames --------------------------
expected_files = cell(30,1);
loading_row = zeros(30,1);
loading_level = zeros(30,1);
r = 0;
for c = 1:6
    for n = 1:5
        r = r+1;
        expected_files{r} = sprintf('830hz_%s_%dn.mat',file_components{c},n);
        loading_row(r) = c;
        loading_level(r) = n;
    end
end

% Match case-insensitively to support variations in the spelling of .MAT.
listing = dir(fullfile(data_dir,'*'));
listing = listing(~[listing.isdir]);
available_names = {listing.name};
resolved_files = cell(30,1);
missing = {};
for r = 1:30
    matches = find(strcmpi(available_names,expected_files{r}));
    if isempty(matches)
        missing{end+1} = expected_files{r}; %#ok<SAGROW>
    else
        assert(numel(matches)==1,'Ambiguous file names matching %s',expected_files{r});
        resolved_files{r} = fullfile(data_dir,available_names{matches});
    end
end
if ~isempty(missing)
    warning('Missing %d expected files:\n%s',numel(missing),strjoin(missing,newline));
    if strict_30_files
        error('Missing expected data files. No CSV has been written.');
    end
end

%% --------------------------- Result table ------------------------------
% For every row: metadata + 6 component magnitudes for each q.
file_col = strings(30,1);
load_col = strings(30,1);
level_col = loading_level;
peak_freq = nan(30,1);
H24_driver_peak = nan(30,1);
Flc_amp = nan(30,1); % signal_Flc single-sided FFT amplitude at SAME peak bin
fs_col = nan(30,1);
N_col = nan(30,1);
fft_df_col = nan(30,1);
H_freq_used = nan(30,4);
amp_all = nan(30,4,6);  % file x H group x component
status = strings(30,1);
if save_mat_details
    details = cell(30,1);
end

%% --------------------------- Batch process -----------------------------
for r = 1:30
    file_col(r) = string(expected_files{r});
    load_col(r) = string(component_names{loading_row(r)});
    if isempty(resolved_files{r})
        status(r) = "Missing file";
        continue;
    end
    fprintf('[%02d/30] %s\n',r,expected_files{r});
    D = load(resolved_files{r});

    % Match original App's H24 reference time axis (signal_in1)
    assert(isfield(D,'signal_in1'),'No signal_in1 in %s',expected_files{r});
    t = double(D.signal_in1.Time(:));
    t = t(isfinite(t));
    t = sort(t);
    t = unique(t,'stable');
    assert(numel(t)>=4,'Insufficient time samples: %s',expected_files{r});
    dtVec = diff(t);
    dt = median(dtVec);
    assert(isfinite(dt) && dt>0,'Invalid dt: %s',expected_files{r});
    fs = 1/dt;
    dtError = max(abs(dtVec-dt))/dt;
    if dtError > 1e-3
        t_used = (t(1):dt:t(end)).';
    else
        t_used = t;
    end
    N = numel(t_used);    % deliberately do NOT drop the final sample
    dt = median(diff(t_used));
    fs = 1/dt;
    df = fs/N;
    nyquist = fs/2;
    assert(peak_band(2)<=nyquist, ...
        'Band [%.1f %.1f] Hz exceeds Nyquist %.3f Hz: %s', ...
        peak_band(1),peak_band(2),nyquist,expected_files{r});
    fs_col(r) = fs;
    N_col(r) = N;
    fft_df_col(r) = df;

    % FFT of all original 24 physical channels, same preprocessing as App.
    U24 = complex(zeros(24,N));
    for ch = 1:24
        name = sprintf('signal_in%d',ch);
        assert(isfield(D,name),'Missing %s in %s',name,expected_files{r});
        ts = D.(name);
        tc = double(ts.Time(:));
        xc = double(ts.Data(:));
        assert(numel(tc)==numel(xc),'Time/Data mismatch: %s / %s',expected_files{r},name);
        keep = isfinite(tc) & isfinite(xc);
        tc = tc(keep);
        xc = xc(keep);
        [tc,order] = sort(tc);
        xc = xc(order);
        [tc,unq] = unique(tc,'stable');
        xc = xc(unq);
        assert(numel(tc)>=2,'Too few samples in %s / %s',expected_files{r},name);
        if numel(tc)~=N || max(abs(tc-t_used))>1e-12
            x = interp1(tc,xc,t_used,'linear','extrap');
        else
            x = xc;
        end
        x = x(:);
        x = x-mean(x);
        U24(ch,:) = fft(x,N).';
    end

    % All candidate FFT bins strictly inside user-provided ~830 Hz band.
    positive_bins = 2:(floor(N/2)+1);
    freq_positive = (positive_bins-1)*df;
    mask = freq_positive>=peak_band(1) & freq_positive<=peak_band(2);
    candidate_bins = positive_bins(mask);
    assert(~isempty(candidate_bins),'No FFT bins in peak band: %s',expected_files{r});

    % Choose the peak from H24 in the known applied loading component.
    driver = loading_row(r);
    model24 = models(4);
    amplitudes_driver = nan(numel(candidate_bins),1);
    for b = 1:numel(candidate_bins)
        idx = candidate_bins(b);
        actual_f = (idx-1)*df;
        [~,hi] = min(abs(model24.freq-actual_f));
        Hk = model24.H(:,:,hi);
        if ~all(isfinite(Hk),'all')
            continue;
        end
        y = Hk * U24(:,idx);
        amplitudes_driver(b) = 2*abs(y(driver))/N;
    end
    if ~any(isfinite(amplitudes_driver))
        status(r) = "No valid H24 within peak band";
        warning('No valid H24 bins for %s',expected_files{r});
        continue;
    end
    amplitudes_driver(~isfinite(amplitudes_driver)) = -inf;
    [peak_value,best_pos] = max(amplitudes_driver);
    best_bin = candidate_bins(best_pos);
    best_freq = (best_bin-1)*df;  % raw physical FFT frequency, never rounded
    peak_freq(r) = best_freq;
    H24_driver_peak(r) = peak_value;

    % ---------- signal_Flc: same time axis, preprocessing and FFT ----------
    % The peak bin is determined by H24 above; do NOT locate a second peak.
    assert(isfield(D,'signal_Flc'), ...
        'Missing signal_Flc in %s',expected_files{r});
    ts_Flc = D.signal_Flc;
    t_Flc = double(ts_Flc.Time(:));
    x_Flc = double(ts_Flc.Data(:));
    assert(numel(t_Flc)==numel(x_Flc), ...
        'signal_Flc Time/Data size mismatch in %s',expected_files{r});
    keep_Flc = isfinite(t_Flc) & isfinite(x_Flc);
    t_Flc = t_Flc(keep_Flc);
    x_Flc = x_Flc(keep_Flc);
    [t_Flc,order_Flc] = sort(t_Flc);
    x_Flc = x_Flc(order_Flc);
    [t_Flc,unique_Flc] = unique(t_Flc,'stable');
    x_Flc = x_Flc(unique_Flc);
    assert(numel(t_Flc)>=2, ...
        'Too few valid signal_Flc samples in %s',expected_files{r});
    if numel(t_Flc)~=N || max(abs(t_Flc-t_used))>1e-12
        x_Flc_used = interp1(t_Flc,x_Flc,t_used,'linear','extrap');
    else
        x_Flc_used = x_Flc;
    end
    x_Flc_used = x_Flc_used(:);
    x_Flc_used = x_Flc_used - mean(x_Flc_used);
    Flc_fft = fft(x_Flc_used,N);
    Flc_complex_peak = Flc_fft(best_bin);
    Flc_amp(r) = 2*abs(Flc_complex_peak)/N;

    % Evaluate all models on EXACTLY the same FFT index / physical frequency.
    y_all = complex(nan(6,4));
    for m = 1:4
        [~,hi] = min(abs(models(m).freq-best_freq));
        Hk = models(m).H(:,:,hi);
        H_freq_used(r,m) = models(m).freq(hi);
        if ~all(isfinite(Hk),'all')
            warning('Invalid H%d at %.6f Hz, file %s',models(m).q,best_freq,expected_files{r});
            continue;
        end
        U = U24(models(m).channels,best_bin);
        y = Hk*U;
        y_all(:,m) = y;
        amp_all(r,m,:) = reshape(2*abs(y)/N,1,1,6);
    end
    if all(isfinite(reshape(amp_all(r,:,:),[],1)))
        status(r) = "OK";
    else
        status(r) = "Missing invalid H value(s)";
    end
    if save_mat_details
        details{r} = struct('source_file',resolved_files{r}, ...
            'peak_fft_bin',best_bin,'peak_frequency_Hz',best_freq, ...
            'complex_FM_6x4',y_all, 'amplitude_FM_6x4',2*abs(y_all)/N, ...
            'signal_Flc_complex_peak',Flc_complex_peak, ...
            'signal_Flc_amplitude',Flc_amp(r));
    end
    fprintf('   applied %s: peak %.6f Hz, H24 magnitude %.9g\n', ...
        component_names{driver},best_freq,peak_value);
end

%% ----------------------------- Export CSV ------------------------------
T = table(file_col,load_col,level_col,peak_freq,H24_driver_peak, ...
    fs_col,N_col,fft_df_col, ...
    'VariableNames',{'File','LoadComponent','LoadLevel_N', ...
    'PeakFreq_Hz','H24_LoadPeakAmplitude','SampleRate_Hz', ...
    'FFT_N','FFT_dF_Hz'});
T.signal_Flc_Amp = Flc_amp;
for m = 1:4
    q = models(m).q;
    T.(sprintf('H%d_UsedFreq_Hz',q)) = H_freq_used(:,m);
    for c = 1:6
        T.(sprintf('H%d_%s',q,component_names{c})) = amp_all(:,m,c);
    end
end
T.Status = status;
writetable(T,output_csv,'Encoding','UTF-8');
fprintf('\nCSV saved: %s\n',output_csv);
fprintf('Rows = %d; OK rows = %d; columns = %d\n', ...
    height(T),nnz(status=="OK"),width(T));
if save_mat_details
    output_mat = fullfile(script_dir,'peak_830_30files_fourH_with_Flc_details.mat');
    save(output_mat,'T','details','models','peak_band','-v7.3');
    fprintf('Additional MAT details: %s\n',output_mat);
end
% No plots: designed for unattended 30-file batch processing.
