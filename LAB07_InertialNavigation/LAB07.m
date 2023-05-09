
%load data without errors
data = load("Inertial_data.dat");
epoch = data(:, 1);
%accelerations
acc = data(:, 2:3); 
%omega in zed
omegaz = data(:,4);
%call function to compute trajectory without errors
traject= CalcTrajectory1(acc, omegaz, epoch);
plot(traject(:,1), traject(:,2),'-b');

%load data with errors (repeat)
data = load("Inertial_data_ni.dat");
epoch = data(:, 1);
%accelerations
acc = data(:, 2:3); 
%omega in zed
omegaz = data(:,4);
%call function to compute trajectory with errors
traject = CalcTrajectory1(acc, omegaz, epoch);

%plot comparison on the same plot
hold on, plot(traject(:,1), traject(:,2),'-b', traject(:,1), traject(:,2),'-r');
legend('Trajectory without error', 'Trajectory error', 'Location','best');
