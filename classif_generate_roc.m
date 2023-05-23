function [x_roc, y_roc, t_roc, auc_roc, targets_all, outputs_all] = classif_generate_roc(data_files, dict_files, det_class, sp_alg, radius, iter_cd, err_thresh)
    targets_all = [];
    outputs_all = [];
    fprintf('%s --- Running classification algorithm %s\n', datetime, sp_alg);
    
    for i_file = 1:numel(data_files)
        fprintf('%s --- Loading file %s\n', datetime, dict_files{i_file});
        load(data_files{i_file})
        load(dict_files{i_file})
    
        % Classification with OMP
        Y = Cwin;

        if strcmp(sp_alg, 'OMP')
            param.L = s;
            X = mexOMP(Y, D, param);
            recerror = sum((Y - D*X).^2);
        elseif strcmp(sp_alg, 'OMP cone')
            % method 1: redistribute radii based on atom utilization
            if 0  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if length(radius) > 1
                param.L = s;
                X = mexOMP(Y, D, param);
                u = sum(X~=0,2);  % utilization of atoms
                [~,ii] = sort(u);
                sort_rad = sort(radius);  % sort increasingly
                radius = sort_rad(ii);
              end
            end

            % method 2: uniform distribution from 0 to radius
            if 0
              if length(radius) == 1
                param.L = s;
                X = mexOMP(Y, D, param);
                u = sum(X~=0,2);  % utilization of atoms
                [~,ii] = sort(u);
                sort_rad = linspace(0, radius, size(X,1)+1);
                sort_rad = sort_rad(2:end)';
                radius = sort_rad(ii);
              end
            end
            
            % method 3: proportional with utilization
            if 0
              if length(radius) == 1
                param.L = s;
                X = mexOMP(Y, D, param);
                u = sum(X~=0,2);  % utilization of atoms
                u = u / max(u);
                radius = u * radius;
              end
            end
            
            % method 0: radius as given
            X = zeros(size(D, 2), size(Y, 2));
            recerror = zeros(1, size(Y,2));
            for i_signal = 1:size(Y, 2)
              [x, Da, support] = omp_cone_single(Y(:, i_signal), D, s, radius, iter_cd, err_thresh);
              X(support, i_signal) = x;
              recerror(i_signal) = norm(Y(:, i_signal) - Da*x)^2;
            end
        else
            error('This case is not covered yet');
        end
    
        % Classify and plot ROC curve
        errors = zeros(size(seg_type));
        targets = ismember(seg_type, det_class)'; % We want 1 when V is detected
        recerror_len = numel(recerror);

        for i = 1:numel(seg_ann)
            segment_i = seg_ann(i) - 127;
            left  = max(segment_i - 100, 1);
            right = min(segment_i + 100, recerror_len);
            if left < right  % it may happen than right < 0 so, skip this 
                errors(i) = median(recerror(left:right));
            else
                errors(i) = 0;  % ignore
            end
        end
        
        outputs = errors' / max(errors);
        targets_all = [targets_all, targets];
        outputs_all = [outputs_all, outputs];
    end

    [x_roc, y_roc, t_roc, auc_roc] = perfcurve(targets_all, outputs_all, 1);
end
