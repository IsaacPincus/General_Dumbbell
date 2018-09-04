function [Q, dQ] = textra(data, iflag, initial_limit)

    x = data(:,1);
    y = data(:,2);
    dy = data(:,3);

    min_err = 100;

    final_coeff = [];
    final_err = [];
    
    limit = initial_limit;
    
    for k=1:6
        for i = 1:length(x)-1
            xi = x(i:end);
            yi = y(i:end);
            dyi = dy(i:end);

            for j = 1:(length(xi)-1)
                if iflag == 1
                    p = j:-1:0;
                elseif iflag == 0
                    p = [j:-1:2,0];
                end
                [coeff, error, chi2val] = pfit2(p, xi, yi, dyi);
                chi_inverse = chi2inv(1-limit, length(xi) - length(p));

                if chi2val < chi_inverse
                    if error(end) < min_err
                        min_err = error(end);

                        final_coeff = coeff;
                        final_err = error;
                    end
                end
            end
        end
        % Do same calculation again with a 10 times smaller limit
        if isempty(final_coeff)
            if k==6
                failed_string = ['Did not converge after %d attempts' ...
                    ', setting Q and dQ to smallest timestep values'];
                sprintf(failed_string, k)
                Q = y(1);
                dQ = dy(1);
                return
            end
            limit = limit/10;
        else
            break;
        end
    end
    
    sprintf('The final limit is %.1d', limit)
    
    Q = final_coeff(end);
    dQ = final_err(end);
end