%% Dataset from the stoner

%Load and save the data
dataset=readtable('G1.txt', "VariableNamingRule", 'preserve');

lat = table2array(dataset(1:480, 'Latitudine'));
lon = table2array(dataset(1:480, 'Longitudine'));

%Convertion of the data from degree to utm
[E, N, utm_zone_S] = deg2utm(lat, lon);

ID = (1:length(E))';
tab = table(ID, E, N);

%writing the table in a file
writetable(tab, 'coord_G1_T1.csv')

%% Matlab Data

%Load and save the data
dataset_m = readtable('sensorlog_pos_20221115_101008.csv');

lat_m = table2array(dataset_m(:, 'latitude'));
long_m = table2array(dataset_m(:, 'longitude'));

%Convertion of the data 
id_m = (1:length(lat_m))';
tab_m = table(id_m, lat_m, long_m);

%writing the table in a file
writetable(tab_m, 'coord_sens.csv')

%Convertion of the data from degree to utm
[E_m, N_m, utm_zone] = deg2utm(lat_m, long_m);

stoner_coord=[E, N];
matlab_coord=[E_m, N_m];

%computing the ditances and save the minimum 
distance = zeros(length(matlab_coord), 1);
for i = 1:length(matlab_coord)
    dist_points=sqrt(sum((repmat(matlab_coord(i,:), [length(stoner_coord), 1]) - stoner_coord).^2, 2));
    distance(i, 1)=min(dist_points);
end


%% Plot of the trajectories
plot(E, N, 'b', E_m, N_m, 'r')








