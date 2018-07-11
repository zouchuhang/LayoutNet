function [wall_d, x_min, f_min] = pano_line_solver_6(lines, w_res, options)

    pk_loc = lines(1:2:end,1);
        
    % rotate
    line_label = lines(1:2:end,3);
    line_rot = find(line_label == 2);

    % initialize based on convex/concave
    if line_rot == 1
        
    else if line_rot == 2
            x_ini = [1.5; 0.5; -1; 2; 1; pk_loc; w_res];
        else if line_rot == 3
                x_ini = [0.5; 1.5; 1; 2; 2; pk_loc; w_res];
            else if line_rot == 4
                    x_ini = [0.25; 0.25; 0.5; 0.5; 1; pk_loc; w_res];
                else if line_rot == 5
                        x_ini = [0.75; 0.25; 1; 0.5; 0.5; pk_loc; w_res];
                        else if line_rot == 6
                            end
                    end
                end
            end
        end 
    end
    
    x_min = Inf;
    f_min = Inf;
    % solve
    for t = 1:400
        %while 1
        f_old = Inf;
        for k = 1:5
            [x,f] = minFunc(@sampleEijOpt_6,x_ini,options);
            x_ini = x;
            % iterate until converge
            if abs(f- f_old) < 1e-5
                break
            end
            f_old = f;
        end
        %disp(f)
        if f < f_min
            f_min = f;
            x_min = x;
        end
        if f < 1e-5
            break
        end
        % TODO: fix initialization here
        rng('shuffle');
        if line_rot == 1
        else if line_rot == 2
                xo = rand(3,1);
                x_ini = [xo(1)+1;xo(2);-xo(3);2;1;pk_loc; w_res];
            else if line_rot == 3
                    xo = rand(2,1)*2;
                    x_ini = [xo(1);xo(2);xo(2)/2;2;2;pk_loc; w_res];
                else if line_rot == 4
                        xo = rand(2,1);
                        x_ini = [xo;(xo(2)+1)/2;(1+xo(1))/2;1;pk_loc; w_res];
                    else if line_rot == 5
                            xo = rand(2,1);
                            x_ini = [xo(1);xo(2);1;xo(1)/2;(1+xo(2))/2;pk_loc; w_res];
                        end
                    end
                end
            end
        end
    end
    disp(f_min)
    %keyboard
    
    % wall d
    cor = [0 0; 1,0; 1, x_min(3); x_min(4), x_min(3);x_min(4), x_min(5); 0, x_min(5)];
    wall_d = (cor - repmat([x_min(1) x_min(2)], 6, 1)); 
    wall_d = wall_d(:,1).^2 + wall_d(:,2).^2;
    wall_d = sqrt(wall_d); 
