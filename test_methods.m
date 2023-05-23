clc;
clear;
close all;

% DL params
s = 4;
n_iter = 40;
n_atoms = 96;
n_segments = 6;
signal_number = 109;

% code params
radius = 0.06;
iter_cd = 10;
err_thresh = 1e-7;
round = 6;

% train all DL methods
dictlearn(signal_number, n_segments, n_atoms, n_iter, s);

% classify by DL methods obtained before
classif(signal_number, n_segments, 'VF', n_atoms, n_iter, s, radius, iter_cd, err_thresh, false)
