function [D, errors, D_hist, X_hist] = K_SVD(Y, D0, n_iter, s)
    D = D0;
    D_hist = cell(1, n_iter);
    X_hist = cell(1, n_iter);
    errors = zeros(1, n_iter);
    
    for i_iter = 1:n_iter
        fprintf('%s --- K-SVD iter %2d', datetime, i_iter);
    
        param.L = s;
        X = mexOMP(Y, D, param);
        errors(i_iter) = norm(Y - D*X,'fro');
        
        fprintf(', residual = %g\n', errors(i_iter));
        
        % Compute residual matrix
        R = Y - D*X;
        
        % Coordinate descent - update each atom one by one
        for i = 1:size(D,2)
            % Keep contribution of current atom
            Rcurr = R + D(:,i)*X(i,:);
            
            % Restrict R to the data having the current atom in the support
            Rcurr = Rcurr(:, X(i,:) ~= 0);
            
            if (size(Rcurr,2) > 0)
                % SVD
                [U,~,~] = svd(Rcurr,'econ');
                D(:,i) = U(:,1);
            end
        end
    
        % Save history
        D_hist{i_iter} = D;
        X_hist{i_iter} = X;
    end
end