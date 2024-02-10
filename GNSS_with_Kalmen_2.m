disp("GNSS kalman filter");

ranges = readmatrix('Pseudo_ranges_remove_outlier.csv');
range_rates = readmatrix('Pseudo_range_rates_remove_outlier.csv');

Define_Constants;
deg_to_rad = 0.01745329252;
PSD = 0.01; % acceleration power spectral density
clock_phase = 0.01; % clock_phase PSD
cf_PSD = 0.04; % clock_frequency PSD
rho_range = 10; % range initial std
rho_rate = 0.05; % rate initial std
co_std = 100000;  % clock offset std
cd_std = 200; % clock drift std
precision_thres = 0.05;
outlier_thres = 6;
%% calculate the least squares GNSS solution for the first epoch, getting rid of outliers
% state vector here: ECEF position, ECEF velocity, clock offset, clock drift
[xhat_prev_plus,P_prev_plus] = Initialise_GNSS_KF(ranges, range_rates, co_std, cd_std, outlier_thres, precision_thres);

epochs = size(ranges,1)-1; % how many time stamps to filter
time_stamp = ranges(2:end,:); 
GNSS_table = zeros(epochs,9);
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(xhat_prev_plus(1:3),xhat_prev_plus(4:6)); % only get the pos(NED) and velocity

%% iteration for each epoch
for i = 1:epochs
    if i ~= 1
        t = ranges(i+1,1) - ranges(i,1); % get tau
    else
        t = 0; % if at time = 0, then run kalman filter with no velocity update
    end
    [xhat_now_plus, P_now_plus] = one_step_filter(xhat_prev_plus, P_prev_plus, t, PSD, clock_phase, cf_PSD, i, ranges, range_rates, rho_range, rho_rate);
    % estimate for the next epoch
    xhat_prev_plus = xhat_now_plus;
    P_prev_plus = P_now_plus;
    % convert to NED and radians
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(xhat_now_plus(1:3),xhat_now_plus(4:6));
    GNSS_table(i, 1) = time_stamp(i);
    GNSS_table(i, 2:4) = [L_b/deg_to_rad, lambda_b/deg_to_rad, h_b];
    GNSS_table(i, 5:7) = v_eb_n.';
    GNSS_table(i, 8:9) = xhat_now_plus(7:8);
end

format longg
% the GNSS solution contains more than [Requirements for Deliverable 1: Motion Profile File]
% we need the alttitude and uncertainty for the first epoch in other methods.
% column 1: time
% column 2-4: position in lattitude
% column 5-7: velocity in m/s
% column 8,9: position and velocity uncertainty
writematrix(GNSS_table, "GNSS_sol.csv");

% functions
function  [xhat_now_plus, P_now_plus] = one_step_filter(xhat_prev_plus, P_prev_plus, t, PSD, clock_phase, cf_PSD, time_idx, ranges, range_rates, rho_range, rho_rate)
    % prev_plus: previous filtered estimation(vectors are column)
    % t: transition time tau
    % PSD: acceleration PSD 
    % clock_phase: clock_phase m^2s^(-3)
    % cf_PSD: clock frequency PSD
    % time_idx, ranges, range_rates: to give the range and rates
    % measurement at the current time
    % rho_range, rho_rate: range and rates measurement std
    % delt_x: std in single dimension of position measurement
    % Pos_now: position estimations at this epoch as column vector
    % Output
    %   xhat_now_plus: the updated estimation of state vector at time stamp (now)
    %   P_now_plus: the: the updated estimation of uncertainty matrix at time stamp (now)
    omega_ie = 7.292115E-5;
    omega = Skew_symmetric([0,0,omega_ie]);
    time_stamp = ranges(2:end,:);
    
    % transition matrix
    Phi = [eye(3),t*eye(3),zeros(3,2);
           zeros(3, 3), eye(3), zeros(3,2);
           zeros(2,6), [1, t; 0, 1]];
    % system noise covariance matrix
    Q = [1/3 * PSD*(t^3) * eye(3), 1/2 * PSD*(t^2) * eye(3), zeros(3,2);
         1/2 * PSD*(t^2) * eye(3), PSD*t*eye(3), zeros(3,2);
         zeros(1,6), t*clock_phase + 1/3 * cf_PSD * (t^3), 1/2 * cf_PSD * (t^2);
         zeros(1,6), 1/2 * cf_PSD * (t^2), cf_PSD * t];
    % propagate the state estimates and error covariance matrix
    xhat_now_minus = Phi * xhat_prev_plus;
    P_now_minus = Phi * P_prev_plus * Phi.' + Q;

    % Predict the ranges from the approximate user position to each satellite
    r_ea_minus = xhat_now_minus(1:3);
    v_ea_minus = xhat_now_minus(4:6);
    sat_r_es_es = [];
    sat_v_es_es = [];
    for satID = ranges(1,2:end)
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time_stamp(time_idx), satID);
        sat_r_es_es = [sat_r_es_es, sat_r_es_e.'];
        sat_v_es_es = [sat_v_es_es, sat_v_es_e.'];
    end
    sat_num = size(ranges,2) - 1;
    r_ases_est = [];
    rr_ases_est = [];
    u_ases_est = [];
    for i = 1:sat_num
        % get
        % estimated position
        % Sagnac effect compensation matrix
        % line-of-sight unit vector
        [r_as, c_ie, u_as] = get_range_correction_and_los(sat_r_es_es(:,i), r_ea_minus);
        rr_as = u_as.' * (c_ie*(sat_v_es_es(:,i) + omega*sat_r_es_es(:,i)) - (v_ea_minus + omega*r_ea_minus)); 
        rr_ases_est = [rr_ases_est, rr_as];
        r_ases_est = [r_ases_est, r_as];
        u_ases_est = [u_ases_est, u_as];
    end

    % Compute the measurement matrix
    H = [-u_ases_est.', zeros(sat_num, 3), ones(sat_num, 1), zeros(sat_num, 1);
         zeros(sat_num, 3), -u_ases_est.', zeros(sat_num, 1), ones(sat_num,1)];
    % Compute the measurement noise covariance matrix
    R = [rho_range^2 * eye(sat_num), zeros(sat_num, sat_num);
         zeros(sat_num, sat_num), rho_rate^2 * eye(sat_num)];
    % Kalmain gain matrix
    K = P_now_minus * H.' / (H*P_now_minus*H.' + R);

    % Formulate the measurement innovation vector
    measured_range = ranges(time_idx+1,2:end);
    measured_rates = range_rates(time_idx+1,2:end);
    delt_z = [measured_range.' - r_ases_est.' - xhat_now_minus(7);
              measured_rates.' - rr_ases_est.' - xhat_now_minus(8)];

    % Update the state estimates and error covariance
    xhat_now_plus = xhat_now_minus + K*delt_z;
    P_now_plus = (eye(size(P_now_minus)) - K*H) * P_now_minus;
end