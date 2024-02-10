disp("Dead reckoning");

Define_Constants;
deg_to_rad = 0.01745329252;

GNSS = readmatrix('GNSS_sol.csv'); % for initial position and heights
dr_reading = readmatrix("Dead_reckoning.csv");

use_gyro = false; % whether to use gyro-magnetometer integrated heading
if use_gyro == true
    disp("Using gyro integrated heading");
    dr_reading = readmatrix("dr_integrated_heading.csv");
    dr_reading(:, 7) = dr_reading(:,8);
    file_name = "DR_sol_ih.csv";
else
    disp("Using magnetic compass heading");
    file_name = "DR_sol_mh.csv";
end

time_steps = size(dr_reading, 1);
% use GNSS solution as initial position
DR_table = zeros(time_steps, 6);
L_last = GNSS(1,2);
lambda_last = GNSS(1,3);
v_eb = GNSS(1,5:7);
inst_v_last = v_eb(1:2).';
DR_table(1,:) = [GNSS(1,1), L_last, lambda_last, v_eb(1:2), dr_reading(1, end)];

for i = 2:time_steps
    last = dr_reading(i-1,:);
    now = dr_reading(i, :);
    ave_v = mean(now(4:5));  % mean of four wheels (2:5) or mean of rear two wheels (4:5)
    t_diff = now(1) - last(1);
    h = GNSS(i-1,4);
    % dead reackoning 
    [L, lambda, Vn, Ve] = one_step_dr(now(7), last(7), ave_v, L_last, lambda_last, h, t_diff, inst_v_last(1:2));
    DR_table(i,:) = [now(1), L, lambda, Vn, Ve, now(end)];
    % update velocity and positon for next epoch
    inst_v_last = [Vn;Ve];
    L_last = L;
    lambda_last = lambda;
end
format longg
writematrix(DR_table, file_name);

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