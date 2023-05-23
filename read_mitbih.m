clc
clear
close all

signal_numbers = [109];
%signal_numbers = [100:109 111:119 121:124 200:203 205 207:210 212:215 217 219:223 228 230:234];

disp('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
wb = waitbar(0, 'Reading samples ECG signal from MIT-BIH Arrhythmia Database');

n_sginals = numel(signal_numbers);
data = cell(1, n_sginals);
Fs = cell(1, n_sginals);
tm = cell(1, n_sginals);
ann = cell(1, n_sginals);
ann_type = cell(1, n_sginals);
ann_subtype = cell(1, n_sginals);
ann_chan = cell(1, n_sginals);
ann_num = cell(1, n_sginals);

for i = 1:numel(signal_numbers)
    signal_name = num2str(signal_numbers(i));
    
    disp(['Reading signal ' signal_name])
    waitbar(i / numel(signal_numbers), wb, ['Reading signal ' signal_name]);
    
    % Read a WFDB record
    [ecg_curr, Fs_curr, tm_curr] = rdsamp(['mitdb/' signal_name]);
    data{i} = ecg_curr; % signal
    Fs{i}   = Fs_curr;  % sampling frequency
    tm{i}   = tm_curr;  % sampling intervals
    
    % Read WFDB annotation
    [ann_curr, type_curr, subtype_curr, chan_curr, num_curr] = rdann(['mitdb/' signal_name],'atr');
    ann{i}     = ann_curr;          % annotation locations
    ann_type{i}    = type_curr;     % annotation types
    ann_subtype{i} = subtype_curr;  % annotation subtype
    ann_chan{i}    = chan_curr;     % annotation channel
    ann_num{i}     = num_curr;      % annotation NUM
end

delete(wb);
save('data/MITDB.mat', 'data', 'Fs', 'tm', 'ann', 'ann_type', 'ann_subtype', 'ann_chan', 'ann_num');
