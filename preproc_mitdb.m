function preproc_mitdb(signal_number)

% Sparse Coding with Anomaly Detection
% Amir Adler, Michael Elad, Yacov Hel-Or, Ehud Rivlin

% clear
% close all
load('data/MITDB.mat');

signal_numbers = [109];
% signal_numbers = [100:109 111:119 121:124 200:203 205 207:210 212:215 217 219:223 228 230:234];

signal_data = data{signal_numbers == signal_number}(:,1);            % signal sample data
signal_Fs = Fs{signal_numbers == signal_number};                     % sampling frequency
signal_ann = ann{signal_numbers == signal_number};                   % vector of annotation positions
signal_type = ann_type{signal_numbers == signal_number};             % vector of annotation values

% signal preprocessing: BP filtering
% [b, a] = butter(5, [1 100]/signal_Fs*2, 'bandpass');
% x = filtfilt(b, a, signal_data);
% plot([signal_data, x])
x = signal_data;

figure; hold on; grid on;
title('Full sample')
plot(signal_data);
plot(signal_ann, signal_data(signal_ann), 'ro');
text(signal_ann, signal_data(signal_ann), signal_type);

% Split signal into segments of 5 minutes
seg_seconds = 300;                % 300 seconds = 5 minutes
seg_len = signal_Fs*seg_seconds;  % number of samples in the segment

% Split sample indices into segments, then use them to segment the data and
% the annotations
signal_ind = 1:numel(signal_data);
signal_ind_matrix = buffer(signal_ind, seg_len);  % Segment with 'buffer()'; Segments = columns

% number of full segments
seg_full_num = floor(size(signal_ind, 2) / seg_len);

% Take each segment
for seg_num = 1:seg_full_num
    if seg_num + 1 == size(signal_ind_matrix, 2)
        seg_idx = [signal_ind_matrix(:, seg_num); signal_ind_matrix(:, seg_num + 1)];
    else
        seg_idx = signal_ind_matrix(:, seg_num); % each segment = each column
    end

    % Ignore padding with 0 at the end of last segment
    seg_idx = seg_idx(seg_idx ~= 0);
    
    % Segment data
    seg_x = x(seg_idx);
    
    % Annotation indices in the segment, indexed from start of segment
    annidx = signal_ann > seg_idx(1) & signal_ann <= seg_idx(end);
    seg_ann = signal_ann(annidx);
    seg_ann = seg_ann - seg_idx(1);
    
    % Annotation values in the segment
    seg_type = signal_type(annidx);
    seg_ann_P = seg_ann(seg_type ~= 'L');
    seg_type_str_P = seg_type(seg_type ~= 'L');

    figure; hold on; grid on;
    title(sprintf('Segment %d', seg_num))
    plot(seg_x);
    plot(seg_ann_P, seg_x(seg_ann_P), 'ro');
    text(seg_ann_P, seg_x(seg_ann_P), seg_type_str_P);
    hold off;

    % Get windows of segment
    winlen = 256;
    win = buffer(seg_x, winlen, winlen-1, 'nodelay');  % buffer with overlap len-1
    P = pca(win');
    PCAsize = 32;
    Pr = P(:,1:PCAsize);
    Cwin = Pr'*win;

    filename = sprintf('data/preproc/preproc_mitdb%d_seg%d.mat', signal_number, seg_num);
    fprintf('Saving %s\n', filename);
    save(filename, 'Cwin', 'Pr', 'win', 'seg_type', 'seg_ann');
end