Dead_reckoning_datas = readtable("Dead_reckoning.csv");
Dead_reckoning_datas.Properties.VariableNames = {'Time(s)', 'Sensor1', 'Sensor2', 'Sensor3', 'Sensor4', 'Gyroscope', 'Magnetic_compass'};
Integral_gyroscope = table(0, 'VariableNames', {'Gyro integral heading'});
Dead_reckoning_datas_modify_heading_with_Kalmen = table(0, 'VariableNames', {'Gyro integrated heading'});
time_steps = size(Dead_reckoning_datas, 1);

% Define the parameter
sigma_psi = 10; % heading uncertainty
sigma_b = 10;   % bias uncertainty
Gyro_random_noise = 0.0001;   % Gyro random noise with power spectral density
Gyro_bias_variation = 0.0001; % Gyro bias variation with PSD
sigma_M = deg2rad(4); % magnetic heading error std
tau = 0.5; % transition time

% to be estimated
delta_psi_g = 0; % gyro-derived heading error
b_g = 0;         % gyro bais

x = [delta_psi_g; b_g]; % state vector
P = [sigma_psi^2, 0; 0, sigma_b]; % initialize error covariance function
transition_matrix = [1, tau; 0, 1]; % Used to update the previous state to the current state

Q = [Gyro_random_noise*tau + Gyro_bias_variation*(tau^3)/3, 0.5*Gyro_bias_variation*(tau^2); 
     0.5*Gyro_bias_variation*(tau^2), Gyro_bias_variation*tau]; % system_noise_covariance_matrix

H = [-1, 0]; % Measurement matrix

R = sigma_M^2; % magnetic heading noise covariance

% initial gyroscope heading should start from magnetic compass reading
gyroscope = deg2rad(Dead_reckoning_datas{1, "Magnetic_compass"});
gyro_heading = zeros(time_steps);
gyro_heading(1) = rad2deg(gyroscope);

for i = 2:time_steps
    % Gyroscope derived heading
    gyroscope = gyroscope + Dead_reckoning_datas{i, "Gyroscope"} * tau;
    gyro_heading(i) = rad2deg(gyroscope);
    tempRow = table(rad2deg(gyroscope), 'VariableNames', {'Gyro integral heading'});
    Integral_gyroscope = [Integral_gyroscope; tempRow];

    % Update step
    x = transition_matrix * x;
    P_now_minus = transition_matrix*P*transition_matrix' + Q;

    %Kalman gain matrix
    Kalmen_filter_matrix = P_now_minus*H' / (H*P_now_minus*H' + R);
    delta_z_k_minus = deg2rad(Dead_reckoning_datas{i, "Magnetic_compass"}) - gyroscope - H*x;

    % Estimate step
    state_estimated = x + Kalmen_filter_matrix * delta_z_k_minus;
    P = (eye(2) - Kalmen_filter_matrix*H) * P_now_minus;

    heading = gyroscope + H*state_estimated;
    gyroColumn = table(rad2deg(heading), 'VariableNames', {'Gyro integrated heading'});
    Dead_reckoning_datas_modify_heading_with_Kalmen = [Dead_reckoning_datas_modify_heading_with_Kalmen; gyroColumn];
    x = state_estimated;
end

Dead_reckoning_datas_modify_heading_with_Kalmen = [Dead_reckoning_datas, Dead_reckoning_datas_modify_heading_with_Kalmen];

writetable(Dead_reckoning_datas_modify_heading_with_Kalmen, "dr_integrated_heading.csv");

figure
plot(Dead_reckoning_datas_modify_heading_with_Kalmen{:,"Time(s)"}, Dead_reckoning_datas_modify_heading_with_Kalmen{:,"Magnetic_compass"}, '-b', 'DisplayName', 'Magnetic heading');
hold on;
% plot(Dead_reckoning_datas_modify_heading_with_Kalmen{:,"Time(s)"}, Integral_gyroscope{:, "Gyro integral heading"}, '-m', 'DisplayName', 'Integral heading');
legend;
plot(Dead_reckoning_datas_modify_heading_with_Kalmen{:,"Time(s)"}, Dead_reckoning_datas_modify_heading_with_Kalmen{:,"Gyro integrated heading"}, '-r',  'DisplayName', 'Integrated heading');
title('Heading');
xlabel('Time (sec)');
ylabel('Heading (°) ');
hold off;







