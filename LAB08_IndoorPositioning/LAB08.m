%% Read database of control points and relevant vector of measurements
load('user_db.txt');
load("control_points_db.txt");

%% User Trajectory

% for each user point
for i = 1:length(user_db)
    % for each control point
    for j = 1:length(control_points_db)
            % find the distance of measurements of the user from the
            % measurements of the control point: (m_u-m_cp)*(m_u-m_cp)
            % identify the minimum
            d(j, 1:3)=[control_points_db(j, 2:3),... 
                abs(min(((user_db(i, 2:end) - control_points_db(j, 4:end)))))];
            
    end
    % attribute to the user the position of the control point with the
    % nearest measurements
    m_index(i, :)=[control_points_db(d(:,3)==min(d(:,3)), 2:3)];
end

% Plot with the movements of the user
plot(m_index(:,1), m_index(:,2),control_points_db(:,2), control_points_db(:,3), 'x')
legend('User Trajectory', 'Control points')
title('User Trajectory')
rectangle('Position',[0 0.5 20 9])
axis([-2 22 -1 11])

%% All data matrix 
% If we want all the data saved in a matrix we can add the k counter
% initialized to 1 and add the following code in the inside for cycle:
%D(k, :)=[d(j, 1:3)];
%k=k+1;