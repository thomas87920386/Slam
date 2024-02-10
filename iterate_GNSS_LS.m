function [r_x_hat, v_x_hat, r_delt_z, H] = iterate_GNSS_LS(r_as_measured, rr_as_measured, sat_r_es_es, sat_rr_es_es, thres, initial)
    % iterate the least squares method for GNSS pos and velocity estimates
    % r_as_measured, rr_as_measured: range and range rates
    % sat_r_es_es, sat_rr_es_es: satellite positions and velocity wrt earth
    % thres: iterating threshold
    % initial: initial r_ea estimate


    % Output:
    %   r_x_hat, v_x_hat: estimated user position and velocity
    %   r_delt_z: innovation vector for position
    %   H: measurement matrix

    % constants
    omega_ie = 7.292115E-5;
    omega = Skew_symmetric([0,0,omega_ie]);
    % if no initial estimate, r, v of antenna to earth is 0
    r_ea_minus = [0; 0; 0];
    v_ea_minus = [0; 0; 0];
    delt_rho_est_minus = 0; % predicted receiver offset 
    delt_rho_r_est_minus = 0; % prediceted rho for rates
    sat_num = size(r_as_measured,2); % number of sattellites
    if nargin == 6
        r_ea_minus = initial(1:end-1);
        delt_rho_est_minus = initial(end-1);
    end

    iter = 0;
    while(iter<10000)
        r_ases_est = [];  % r_as
        rr_ases_est = []; % rr_as
        u_ases_est = [];  % line of sight
        for i = 1:sat_num
            % get line of sight, range and rates prediction
            [r_as, c_ie, u_as] = get_range_correction_and_los(sat_r_es_es(:,i), r_ea_minus);
            rr_as = u_as.' * (c_ie*(sat_rr_es_es(:,i) + omega*sat_r_es_es(:,i)) - (v_ea_minus + omega*r_ea_minus)); 
            rr_ases_est = [rr_ases_est, rr_as];
            r_ases_est = [r_ases_est, r_as];
            u_ases_est = [u_ases_est, u_as];
        end
        % get position
        r_x_hat_minus = [r_ea_minus; delt_rho_est_minus];
        r_delt_z_minus = (r_as_measured - r_ases_est - delt_rho_est_minus).';
        H_measurement = [-u_ases_est.' , ones(sat_num,1)];
        r_x_hat_plus = r_x_hat_minus + pinv(H_measurement) * r_delt_z_minus;

        % get velocity
        v_x_hat_minus = [v_ea_minus; delt_rho_r_est_minus];
        v_delt_z_minus = (rr_as_measured - rr_ases_est - delt_rho_r_est_minus).';
        v_x_hat_plus = v_x_hat_minus + pinv(H_measurement) * v_delt_z_minus;
        % stop at threshold
        if sum((r_ea_minus - r_x_hat_plus(1:end-1)).^2) < thres
            break
        end
        % plus is the new old(minus)
        r_ea_minus = r_x_hat_plus(1:end-1);
        delt_rho_est_minus = r_x_hat_plus(end-1);
        v_ea_minus = v_x_hat_plus(1:end-1);
        delt_rho_r_est_minus = v_x_hat_plus(end-1);
        iter = iter + 1;
    end
    H = H_measurement;
    r_delt_z = r_delt_z_minus;
    r_x_hat = r_x_hat_plus;
    v_x_hat = v_x_hat_plus;
end