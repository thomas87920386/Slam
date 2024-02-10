GNSS_sol = readmatrix("GNSS_sol.csv");
DR_table_ih = readmatrix("DR_sol_ih.csv");
DR_table_mh = readmatrix("DR_sol_mh.csv");
integrated_sol_closed = readmatrix("integrated_sol_closed.csv");
integrated_sol_open = readmatrix("integrated_sol_open.csv");

%% GNSS
GNSS_Lattitude = GNSS_sol(:, 2);
GNSS_Longitude = GNSS_sol(:, 3);
GNSS_Vn = GNSS_sol(:, 5);
GNSS_Ve = GNSS_sol(:, 6);

figure;
plot(GNSS_Longitude, GNSS_Lattitude, '-r');
title('GNSS Position');
xlabel('Longitude (°)');
ylabel('Lattitude (°)');
% saveas(gcf,'figures/GNSS.png')

figure;
quiver(GNSS_Longitude, GNSS_Lattitude, GNSS_Ve, GNSS_Vn);
title('GNSS velocity');
xlabel('Longitude (°)');
ylabel('Lattitude (°)');
% saveas(gcf,'figures/GNSS_vel.png')

%% Dead reckoning
DR_Lattitude = DR_table_ih(:, 2); % integrated heading
DR_Longitude = DR_table_ih(:, 3);
DR_Vn = DR_table_ih(:, 4);
DR_Ve = DR_table_ih(:, 5);

DR_Lattitude_mh = DR_table_mh(:, 2); % magnetic heading
DR_Longitude_mh = DR_table_mh(:, 3);
DR_Vn_mh = DR_table_mh(:, 4);
DR_Ve_mh = DR_table_mh(:, 5);

figure;
plot(DR_Longitude, DR_Lattitude, '-k', 'DisplayName', 'DR position using integrated heading');
hold on;
plot(DR_Longitude_mh, DR_Lattitude_mh, '-g', 'DisplayName', 'DR position using magnetic heading');
legend;
title('DR Position');
xlabel('Longitude (°)');
ylabel('Lattitude (°)');
hold off;
% saveas(gcf,'figures/dr.png')

figure;
quiver(DR_Longitude, DR_Lattitude, DR_Ve, DR_Vn, 'DisplayName', 'DR velocity using integrated heading');
hold on;
quiver(DR_Longitude_mh, DR_Lattitude_mh, DR_Ve_mh, DR_Vn_mh, 'DisplayName', 'DR velocity using magnetic heading');
hold off;
legend;
title('DR velocity');
xlabel('Longitude (°)');
ylabel('Lattitude (°)');
% saveas(gcf,'figures/dr_vel.png')

%% integrated

integrated_Lattitude_c = integrated_sol_closed(:, 2);
integrated_Longitude_c = integrated_sol_closed(:, 3);
integrated_Lattitude_o = integrated_sol_open(:, 2);
integrated_Longitude_o = integrated_sol_open(:, 3);
integrated_Vn_c = integrated_sol_closed(:, 4);
integrated_Ve_c = integrated_sol_closed(:, 5);
integrated_Vn_o = integrated_sol_open(:, 4);
integrated_Ve_o = integrated_sol_open(:, 5);

figure;
plot(integrated_Longitude_c, integrated_Lattitude_c, 'DisplayName', 'integrated position using closed loop');
hold on;
plot(integrated_Longitude_o, integrated_Lattitude_o, 'DisplayName', 'integrated position using open loop');
hold off;
legend;
title('Integrated Position');
xlabel('Longitude (°)');
ylabel('Lattitude (°)');
% saveas(gcf,'figures/integrated.png')

figure;
quiver(integrated_Longitude_c, integrated_Lattitude_c, integrated_Ve_c, integrated_Vn_c, 'DisplayName', 'integrated velocity using closed loop');
hold on;
quiver(integrated_Longitude_o, integrated_Lattitude_o, integrated_Ve_o, integrated_Vn_o, 'DisplayName', 'integrated velocity using closed loop');
hold off;
legend;
title('Integrated velocity');
xlabel('Longitude (°)');
ylabel('Lattitude (°)');
% saveas(gcf,'figures/dr_vel.png')



%% todo
% cite
% report
% gyro
% figure for closed loop open loop
% :)