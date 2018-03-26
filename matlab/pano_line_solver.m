function [wall_d, x_min, f_min] = pano_line_solver(lines, w_res, options)

    pk_loc = lines(1:2:end,1);
    
    x_ini = [0.5; 0.5; 1; pk_loc; w_res];
    x_min = Inf;
    f_min = Inf;
    % solve
    for t = 1:100 %200
        f_old = Inf;
        for k = 1:5
            [x,f] = minFunc(@sampleEijOpt,x_ini,options);
            x_ini = max(x,0.01);
            % iterate until converge
            if abs(f- f_old) < 1e-5
                break
            end
            f_old = f;
        end
        if f < f_min
            f_min = f;
            x_min = x;
        end
        if f < 1e-5
            break
        end
        x_ini = [rand(2,1);1;pk_loc; w_res];
    end
    
    % wall d
    cor = [0 0; 1,0; 1, x_min(3); 0, x_min(3)];
    wall_d = (cor - repmat([x_min(1) x_min(2)], 4, 1)); 
    wall_d = wall_d(:,1).^2 + wall_d(:,2).^2;
    wall_d = sqrt(wall_d); 
