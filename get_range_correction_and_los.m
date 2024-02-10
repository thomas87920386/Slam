function [r_as,c_ie, u_as] = get_range_correction_and_los(r_es, r_ea)
    % get as from es and ea
    % r_es, r_ea: column vector, estimated earth to satellite and earth to
    % antenna
    % output:
    % r_as: antenna to satellite (estimated)
    % c_ie: correction matrix
    % u_as: line of sight vector

    c = 299792458; % Speed of light in m/s
    omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
    % first set Sagnax effect = eye
    C_ie = eye(3);
    for i = 1:10
        relative_pos = C_ie * r_es - r_ea;
        r_as = norm(relative_pos);
        C_ie(1, 2) = omega_ie*r_as/c;
        C_ie(2, 1) = -omega_ie*r_as/c;
    end
    c_ie = C_ie;
    u_as = relative_pos / r_as;
end