function [x_est,P_matrix] = Initialise_GNSS_KF(ranges, range_rates, co_std, cd_std, outlier_thres, precision_thres)
    % Runs one epoch GNSS
    % co_std: clock offset std
    % cd_std: clock drift std
    % outlier_thres: threshold for outliers abs(v) > sqrt(C) * thres,
    % precision_thres: prev - current < thres
    % Outputs:
    %   x_est                 Kalman filter estimates:
    %     Rows 1-3            estimated ECEF user position (m)
    %     Rows 4-6            estimated ECEF user velocity (m/s)
    %     Row 7               estimated receiver clock offset (m) 
    %     Row 8               estimated receiver clock drift (m/s)
    %   P_matrix              state estimation error covariance matrix
    
    
    % constants:
    sigma_range = 10;
    sigma_rates = 0.05;
    
    % time_step relaitve info (loop removes from these for outliers)
    i=1;
    satIDs = ranges(1,2:end);
    sat_num = size(ranges,2)-1; % number of sattellites
    time_stamps = ranges(2:end,1); 
    r_as_measured = ranges(i+1,2:end); %first line is satID
    rr_as_measured = range_rates(i+1,2:end);
    
    while true % will break out the loop if there are no outliers
        % get satellite positions
        sat_r_es_es = [];
        sat_v_es_es = [];
        for satID = satIDs
            [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time_stamps(i), satID);
            sat_r_es_es = [sat_r_es_es, sat_r_es_e.'];
            sat_v_es_es = [sat_v_es_es, sat_v_es_e.'];
        end
        % calculate the position and velocity using least squares estimate
        % also get the innovation vector for detecting outliers
        [r_xhat, v_xhat, delt_z_minus, H_measurement] = iterate_GNSS_LS(r_as_measured, rr_as_measured, sat_r_es_es, sat_v_es_es, precision_thres);
        
        % outlier detection
        v = (H_measurement * pinv(H_measurement) - eye(sat_num)) * delt_z_minus; % residuals vector
        Cv = (eye(sat_num) - H_measurement * pinv(H_measurement)) * power(sigma_range,2); % residuals covariance
        
        % exceeding threshold means there is an outlier
        outliers = false;
        for idx = 1:sat_num
            if abs(v(idx)) > sqrt(Cv(idx,idx)) * outlier_thres
                outliers = true;
        %         disp([time_stamps(i), satID(idx), abs(v(idx))]);
            end
        end
        
        % if no outliers, then finish loop
        if outliers == false
            % assign values
            x_est =  [r_xhat(1:end-1); v_xhat(1:end-1);r_xhat(end);v_xhat(end)]; 
            
            P_matrix =  zeros(8);
            P_matrix(1,1) = co_std;
            P_matrix(2,2) = co_std;
            P_matrix(3,3) = co_std;
            P_matrix(4,4) = cd_std;
            P_matrix(5,5) = cd_std;
            P_matrix(6,6) = cd_std;
            P_matrix(7,7) = co_std;
            P_matrix(8,8) = cd_std;
            break;
        end
        disp("remove outlier");
        [~, outlier_sat_idx] = max(abs(v));
        
        disp(outlier_sat_idx);
        satIDs(outlier_sat_idx) = [];
        sat_num = sat_num-1; % number of sattellites
        r_as_measured(outlier_sat_idx) = [];
        rr_as_measured(outlier_sat_idx) = [];
    
    end
end

