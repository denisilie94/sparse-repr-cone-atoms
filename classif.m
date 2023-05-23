function classif(signal_number, segments, det_class, n_atoms, n_iter, s, radius, iter_cd, err_thresh, verbose)
    format_values = [signal_number*ones(segments, 1)  (1:segments)'];
    data_files = sprintfc('data/preproc/preproc_mitdb%d_seg%d.mat', format_values);
    
    format_values = [signal_number*ones(segments, 1) (1:segments)' n_atoms*ones(segments,1) n_iter*ones(segments,1) s*ones(segments, 1)];
    param_strings = {...
        'data/dicts/mitdb_%d_seg_%d_Ksvd_N_%d_iter_%d_s_%d.mat', 'OMP', 'OMP + KSVD';
        'data/dicts/mitdb_%d_seg_%d_Ksvd_N_%d_iter_%d_s_%d.mat', 'OMP cone', 'OMP cone + KSVD';
    };

    n_param_strings = size(param_strings, 1);
    X_roc = cell(1, n_param_strings);
    Y_roc = cell(1, n_param_strings);
    T_roc = cell(1, n_param_strings);
    AUC_roc = cell(1, n_param_strings);
    targets_all = cell(1, n_param_strings);
    outputs_all = cell(1, n_param_strings);
    
    fig = figure; hold on; grid on;
    for i = 1:n_param_strings
        dict_files = sprintfc(param_strings{i}, format_values);
        [X_roc{i}, Y_roc{i}, T_roc{i}, AUC_roc{i}, targets_all{i}, outputs_all{i}] = classif_generate_roc(data_files, dict_files, det_class, param_strings{i,2}, radius, iter_cd, err_thresh);
        
        figure(fig);
        plot(X_roc{i}, Y_roc{i});
        fprintf('ROC AUC%d = %f\n', i, AUC_roc{i});
        ip_100 = find(Y_roc{i} == 1, 1);  % first value of 1, all true positives found
        fprintf('FP at 1   : %d\n', length(unique(X_roc{i}(1:ip_100)))-1);
        ip_97 = find(Y_roc{i} > 0.97, 1); % same for 97%
        fprintf('FP at 0.97: %d\n', length(unique(X_roc{i}(1:ip_97)))-1);
        ip_94 = find(Y_roc{i} > 0.94, 1); % same for 94%
        fprintf('FP at 0.94: %d\n', length(unique(X_roc{i}(1:ip_94)))-1);

        ip = find(targets_all{i}==1);
        figure, stem(outputs_all{i}, 'k.'), hold on, stem(ip, outputs_all{i}(ip), 'c.')
        plot([1 length(outputs_all{i})], [1 1]*min(outputs_all{i}(ip)), 'r'), hold off
        axis([0 length(outputs_all{i}) 0 1])
        xlabel('#heartbeat', 'FontSize', 14)
        ylabel('score', 'FontSize', 14)
        grid
        print(['fig_err_109_', num2str(i), '.eps'], '-depsc') 

        if verbose
            file = fopen('logs.txt', 'a');
            fprintf(file, [param_strings{i,3} '\n']);
            fprintf(file, 'ROC AUC%d = %f\n', i, AUC_roc{i});
            fprintf(file, 'FP at 1   : %d\n', length(unique(X_roc{i}(1:ip_100)))-1);
            fprintf(file, 'FP at 0.97: %d\n', length(unique(X_roc{i}(1:ip_97)))-1);
            fprintf(file, 'FP at 0.94: %d\n', length(unique(X_roc{i}(1:ip_94)))-1);
            fprintf(file, '-------------------\n');
            fclose(file);
        end
    end
    legend(param_strings{:, 3})
    
    if verbose
        file = fopen('logs.txt', 'a');
        fprintf(file, ['n_atoms = ' num2str(n_atoms) '\n']);
        fprintf(file, ['n_iter = ' num2str(n_iter) '\n']);
        fprintf(file, ['s = ' num2str(s) '\n']);
        fprintf(file, ['radius = ' num2str(radius) '\n']);
        fprintf(file, ['iter_cd = ' num2str(iter_cd) '\n']);
        fprintf(file, ['err_thresh = ' num2str(err_thresh) '\n']);
        fprintf(file, '===================\n');
        fclose(file);
    end
end
