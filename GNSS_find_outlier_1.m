format longg
Define_Constants; % To utilize certain parameters that have been defined in the .m file.

s_ranges_data = readtable('Pseudo_ranges.csv'); % The data include the outlier
% s_ranges_data = readtable('Workshop1_Pseudo_ranges_remove_outlier.csv'); % The data remove the outlier

user_position_table = table([], [], [], [], 'VariableName', {'Time(s)', 'Latitude', 'Longitude', 'Height'});
outlier_table = table([], [], [], 'VariableName', {'Epoch index', 'Problem satellite', 'Outlier'});

for state_of_time_index = 2: size(s_ranges_data, 1)
    % 0. show the current step
    state_of_time_index

    % 1. Initial User Positions(assumption)
    r_ea_e = [0; 0; 0];

    % 2. ECEF Satellite Positions
    time = s_ranges_data{state_of_time_index, 1};
    s_position_table = table([], [], [], [], 'VariableNames', {'Satellite', 'x', 'y', 'z'});
    for number_of_satellite = 2: size(s_ranges_data, 2)
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time, s_ranges_data{1, number_of_satellite});
        newRow = table(s_ranges_data{1, number_of_satellite}, sat_r_es_e(1), sat_r_es_e(2), sat_r_es_e(3), 'VariableNames', {'Satellite', 'x', 'y', 'z'});
        s_position_table = [s_position_table; newRow];
    end

    % 3. Iteration for the final result
    % initial condition about x_caret_minus
    epsilon_rou_caret_correction_a_minus = 0;
    x_caret_minus = [r_ea_e; epsilon_rou_caret_correction_a_minus];
    
    for count = 1:100
        epsilon_z_minus = [];
        for number_of_satellite = 1:size(s_ranges_data, 2) - 1
            caret_distance_table_from_a_to_s = prediction_distance_calculator(s_position_table, x_caret_minus(1:3));
            element = s_ranges_data{state_of_time_index, number_of_satellite + 1} - caret_distance_table_from_a_to_s{number_of_satellite, 2} - epsilon_rou_caret_correction_a_minus;
            epsilon_z_minus = [epsilon_z_minus; element];
        end
        H_G_e = unit_vector_calculater(s_ranges_data, s_position_table, x_caret_minus(1:3));

        % Update the state
        x_caret_plus = x_caret_minus + inv(H_G_e' * H_G_e) * H_G_e' * epsilon_z_minus;
        if (x_caret_plus - x_caret_minus) < 0.001
            break;
        end
        x_caret_minus = x_caret_plus;
    end
    x_caret_result = x_caret_minus;

    % 3. Detect the Outlier
    residuals_vector = residuals_vector_calculator(H_G_e, epsilon_z_minus);
    residuals_covariance_matrix = residuals_covariance_matrix_calculator(H_G_e);
    for index_of_satellite_in_table = 2: size(s_ranges_data, 2)
        if abs(residuals_vector(index_of_satellite_in_table-1)) > sqrt(residuals_covariance_matrix(index_of_satellite_in_table-1, index_of_satellite_in_table-1))*T
            epoch_index = s_ranges_data{state_of_time_index, 1};
            problem_satellite = s_ranges_data{1, index_of_satellite_in_table};
            outlier = residuals_vector(index_of_satellite_in_table-1);
            tempRow = table(epoch_index, problem_satellite, outlier, 'VariableName', {'Epoch index', 'Problem satellite', 'Outlier'});
            outlier_table = [outlier_table; tempRow];
        end
    end

    % 4. Transform the ECEF coordinates to NED expression
    r_eb_e = x_caret_plus(1:3);
    v_eb_e = [0; 0; 0];
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_eb_e,v_eb_e); % the unit of return value is radien
    latitude = round(L_b / pi * 180, 8, 'significant');
    longtitude = round(lambda_b / pi * 180, 8, 'significant');
    height = round(h_b, 4, 'significant');
    user_position_table_element = table(s_ranges_data{state_of_time_index, 1}, latitude, longtitude, height, 'VariableNames', {'Time(s)', 'Latitude', 'Longitude', 'Height'});
    user_position_table = [user_position_table; user_position_table_element];
end

%% The prediction value of the distance from user to satellite
function caret_distance_table_from_a_to_s = prediction_distance_calculator(s_position_table, r_ea_e)
Define_Constants;
caret_distance_table_from_a_to_s = table([], [], 'VariableNames', {'Satellite', 'Distances'});
for number_of_satellite = 1:size(s_position_table, 1)
    r_ej_e = [s_position_table{number_of_satellite, 2}; s_position_table{number_of_satellite, 3}; s_position_table{number_of_satellite, 4}];
    r_aj = sqrt((r_ej_e - r_ea_e)' * (r_ej_e - r_ea_e));
    % Sagnac effect
    CI_e = [1, omega_ie*r_aj/c, 0; -omega_ie*r_aj/c, 1, 0; 0, 0, 1];
    r_Is_I = CI_e * r_ej_e;
    range_as = sqrt((r_Is_I - r_ea_e)' * (r_Is_I - r_ea_e));
    distance_table_element = table(s_position_table{number_of_satellite, 1}, range_as, 'VariableNames', {'Satellite', 'Distances'});
    caret_distance_table_from_a_to_s = [caret_distance_table_from_a_to_s; distance_table_element];
end
end

%% Calculate the H martix
function H_G_e = unit_vector_calculater(s_ranges_data, s_position_table, r_ea_e)
Define_Constants;
unit_vector_as_table = table([], [], [], [], 'VariableNames', {'Satellite', 'partial X', 'partial Y', 'partial z'});
for number_of_satellite = 1:size(s_ranges_data, 2)-1
    r_ej_e = [s_position_table{number_of_satellite, "x"}; s_position_table{number_of_satellite, "y"}; s_position_table{number_of_satellite, "z"}];
    vector_as = r_ea_e - r_ej_e;
    r_aj_unmodified = sqrt(vector_as' * vector_as);

    % Sagnac effect modified
    CI_e = [1, omega_ie*r_aj_unmodified/c, 0; -omega_ie*r_aj_unmodified/c, 1, 0; 0, 0, 1];
    r_Is_I = CI_e * r_ej_e;

    % Calculate the unit vector
    vector_as = r_Is_I - r_ea_e;
    u_aj_e = vector_as / sqrt(vector_as' * vector_as);

    % Store the value
    x = round(u_aj_e(1), 6, 'significant');
    y = round(u_aj_e(2), 6, 'significant');
    z = round(u_aj_e(3), 6, 'significant');

    u_vector_element = table(s_ranges_data{1, number_of_satellite+1}, x, y, z, 'VariableNames',{'Satellite', 'partial X', 'partial Y', 'partial z'});
    unit_vector_as_table = [unit_vector_as_table; u_vector_element];
end
H_G_e = -1 * unit_vector_as_table{:, 2:4};
H_G_e = [H_G_e, ones(size(H_G_e,1), 1)];
end

%% Compute the residuals vector 
function residuals_vector = residuals_vector_calculator(H_G_e, epsilon_z_minus)
residuals_vector = (H_G_e * inv(H_G_e'*H_G_e) * H_G_e' - eye(size(H_G_e, 1))) * epsilon_z_minus;
end

%% Compute the residuals covariance matrix
function residuals_covariance_matrix = residuals_covariance_matrix_calculator(H_G_e)
Define_Constants;
residuals_covariance_matrix = (eye(size(H_G_e, 1)) - H_G_e*inv(H_G_e'*H_G_e)*H_G_e')*sigma_rho^2;
end
