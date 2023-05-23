function dictlearn(signal_number, segments, n_atoms, n_iter, s)
    data_files = cell(1, segments);
    for seg_num = 1:segments
        data_files{seg_num} = sprintf('data/preproc/preproc_mitdb%d_seg%d.mat', signal_number, seg_num);
    end
    
    for ifile = 1:numel(data_files)
        % Load segment data
        load(data_files{ifile});  % load 'Cwin' matrix
        fprintf('Loaded %s\n', data_files{ifile});

        % Initialize dictionary with random linear combinations, normalized
        Y = Cwin;
        rng(ifile)
        D0 = Y * randn(size(Y,2), n_atoms);
        for i = 1:size(D0, 2)
            D0(:,i) = D0(:,i) / norm(D0(:,i));
        end
    
        % Normal K-SVD with OMP
        [D, errors, D_hist, ~] = K_SVD(Y, D0, n_iter, s);
        savename = [
            'data/dicts/mitdb_' num2str(signal_number)...
            '_seg_' num2str(ifile)...
            '_Ksvd_N_' num2str(n_atoms)...
            '_iter_' num2str(n_iter)...
            '_s_' num2str(s)...
            '.mat'
        ];
        save(savename, 'D', 'D_hist', 'errors', 'n_iter', 's');
    end
end