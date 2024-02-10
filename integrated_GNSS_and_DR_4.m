disp("integrated:");

GNSS_sol = readmatrix("GNSS_sol.csv");
DR_table = readmatrix("DR_sol_ih.csv");
dr_reading = readmatrix("Dead_reckoning.csv");

closed = false; % closed or open loop

Define_Constants;
deg_to_rad = 0.01745329252;
PSD = 0.01; % acceleration power spectral density
Sdr = 0.01; % velocity psd 
clock_phase = 0.01; % clock_phase PSD
cf_PSD = 0.04; % clock_frequency PSD
rho_range = 10; % range measurement std
rho_rate = 0.05; % rate measurement std
% sigma_v, sigma_r: velocity position uncertainty in both drections
sigma_v = GNSS_sol(1,9);
sigma_r = GNSS_sol(1,8);

% measurement matrix (velocity first in state vector)
H = [0, 0, -1, 0;
     0, 0, 0, -1;
     -1, 0, 0, 0;
     0, -1, 0, 0];

%% initial epoch
% get first epoch from GNSS solution
h0 = GNSS_sol(1,4);
L0 = GNSS_sol(1,2)*deg_to_rad;

[Rn, Re] = Radii_of_curvature(L0);
Var_n = sigma_r^2/(Rn+h0)^2;
Var_e = sigma_r^2/((Re+h0)^2 * cos(L0)^2);

% x is [delta_vn, delta_ve, delta_L, delta_lambda]
x_last_plus = zeros(4,1); % initially we do not know the DR error
% state estimation error covariance matrix
P_last_plus = [sigma_v^2*eye(2), zeros(2,2);
               zeros(2,2), [Var_n, 0;0,Var_e]];

time_steps = size(DR_table, 1);
integrated_table = zeros(time_steps, 6);
integrated_table(1,:) = [GNSS_sol(1,1), GNSS_sol(1,2), GNSS_sol(1, 3), GNSS_sol(1,5), GNSS_sol(1,6), dr_reading(1, end)];

for i = 2:time_steps
    
    h_last = GNSS_sol(i-1,4);
    L_last = GNSS_sol(i-1,2)*deg_to_rad;
    h_now = GNSS_sol(i,4);
    L_now = GNSS_sol(i,2)*deg_to_rad;

    t_diff = GNSS_sol(i,1) - GNSS_sol(i-1,1);
    % GNSS solution [L, lambda, Vn, Ve]
    G = [GNSS_sol(i, 2)*deg_to_rad; GNSS_sol(i, 3)*deg_to_rad; GNSS_sol(i, 5); GNSS_sol(i, 6)];

    if closed == true
        % closed loop integration
        L_last_inte = integrated_table(i-1,2);
        lambda_last_inte = integrated_table(i-1,3);
        h_inte = GNSS_sol(i-1,4);
        v_eb = integrated_table(i-1,4:5).';

        last = dr_reading(i-1,:);
        now = dr_reading(i, :);
        ave_v = mean(now(4:5)); % mean of four wheels (2:5) or mean of rear two wheels (4:5)
        [L, lambda, Vn, Ve] = one_step_dr(now(7), last(7), ave_v, L_last_inte, lambda_last_inte, h_inte, t_diff, v_eb);
        D = [L*deg_to_rad; lambda*deg_to_rad; Vn; Ve];
        % zero out the kalman filter dr error estimate
        x_last_plus = zeros(4,1);
    else
        D = [DR_table(i, 2)*deg_to_rad; DR_table(i, 3)*deg_to_rad; DR_table(i, 4); DR_table(i, 5)];
    end

    % one iteration of integration
    [x_now_plus, P_now_plus] = one_step_integration(h_last, L_last, t_diff, h_now, L_now, rho_range, rho_rate, Sdr, G, D, H, x_last_plus, P_last_plus);
    % update new epoch 
    x_last_plus = x_now_plus;
    P_last_plus = P_now_plus;
    
    % dead reckoning - error estimates(velocity, position)
    % same as dr + H*state vector
    C = D+H*x_now_plus;
    C(1) = C(1)/deg_to_rad;
    C(2) = C(2)/deg_to_rad;
    integrated_table(i,:) = [GNSS_sol(i,1), C.', dr_reading(i, end)];
end

format longg
if closed==true
    disp("Using closed loop GNSS DR integration: ");
    writematrix(integrated_table, "integrated_sol_closed.csv");
else
    disp("Using open loop GNSS DR integration: ");
    writematrix(integrated_table, "integrated_sol_open.csv");
end

% functions
function [L, lambda, Vn, Ve] = one_step_dr(Psi_now, Psi_last, ave_v_now, L_last, lambda_last, h, t_diff, inst_v_last)
    % L, lambda in degrees
    % last: last epoch
    % t_diff: tau, transition time
    % ave_v_now: average speed for this epoch
    % Psi in degrees
    % Vn, Ve damped instantneous DR velocity
    % inst_v_last is column vector
    % Output
    %   L, lambda: lattitude, longitude
    %   Vn, Ve: instant velocity in north east direction
    deg_to_rad = 0.01745329252;
    % deg2rad
    Psi_now_rad = deg_to_rad*Psi_now;
    Psi_last_rad = deg_to_rad*Psi_last;
    L_last_rad = L_last*deg_to_rad;
    lambda_last_rad = lambda_last*deg_to_rad;
    % speed to velocity
    ave_velocity = 0.5*[cos(Psi_now_rad) + cos(Psi_last_rad); sin(Psi_now_rad) + sin(Psi_last_rad)]*ave_v_now;
    
    % get meridian and transverse radius of curvature
    [Rn, Re] = Radii_of_curvature(L_last_rad);

    % update position
    L_diff = ave_velocity(1)*t_diff/(Rn+h);
    L_rad = L_last_rad + L_diff;

    lambda_diff = ave_velocity(2)*t_diff/(Re+h)/cos(L_rad);
    lambda_rad = lambda_last_rad + lambda_diff;

    L = L_rad/deg_to_rad;
    lambda = lambda_rad/deg_to_rad;
    % instant speed calculation
    inst_v_now = 1.7*ave_velocity - 0.7*inst_v_last;
    Vn = inst_v_now(1);
    Ve = inst_v_now(2);
end

function [x_now_plus, P_now_plus] = one_step_integration(h_last, L_last, t_diff, h_now, L_now, sigma_Gr, sigma_Gv, Sdr, G, D, H, x_last_plus, P_last_plus)
    % h_last, L_last: height and lattitude of last time-step
    % t_diff: propagation time between time-steps
    % h_now, L_now: height and lattitude at this time-step
    % sigma_Gr, sigma_Gv: measurement error std of velocity and position
    % Sdr: DR velocity error ower spectral density
    % G, D: The column vector of GNSS and DR solutions in Pn;Pe;Vn;Ve
    % H: measurement matrix
    % x,P_last_plus: from last estimation of x and P
    
    [Rn, Re] = Radii_of_curvature(L_last);
    
    % transition matrix
    Phi = eye(4);
    Phi(3,1) = t_diff/(Rn+h_last);
    Phi(4,2) = t_diff/(Re+h_last)/cos(L_last);
    % system noise
    Q = [Sdr*t_diff, 0, 0.5*Sdr*t_diff^2/(Rn+h_last), 0;
         0, Sdr*t_diff, 0, 0.5*Sdr*t_diff^2/(Re+h_last)/cos(L_last);
         0.5*Sdr*t_diff^2/(Rn+h_last), 0, 1/3*Sdr*t_diff^3/(Rn+h_last)^2, 0;
         0, 0.5*Sdr*t_diff^2/(Re+h_last)/cos(L_last), 0, 1/3*Sdr*t_diff^3/(Re+h_last)^2/cos(L_last)^2];
    
    % propagation
    x_now_minus = Phi*x_last_plus;
    P_now_minus = Phi*P_last_plus*Phi.' + Q;
    
    % measurement noise covariance
    [Rn, Re] = Radii_of_curvature(L_now);
    R_diag = [sigma_Gr^2/(Rn+h_now)^2, sigma_Gr^2/(Re+h_now)^2/cos(L_now)^2, sigma_Gv^2, sigma_Gv^2];
    R = diag(R_diag);
    
    % kalman gain matrix
    K = P_now_minus*H.' / (H*P_now_minus*H.' + R);
    
    % measurement innovation vector
    delt_z_now_minus = G-D-H*x_now_minus;
    
    % updates
    x_now_plus = x_now_minus + K*delt_z_now_minus;
    P_now_plus = (eye(4) - K*H) * P_now_minus;
end
