function[x_y] = CalcTrajectory1(acc, omegaz, epoch)
    %initialize velocities and delta positions in body frame
    a_x = acc(:,1);
    a_y = acc(:,2);
    alpha = 0;
    vx = 0;
    vy = 0;
    x_y = [100 100];
    x = x_y(1);
    y = x_y(2);
    dx = 0;
    dy = dx;
    centr_y=0;

for i = 1 : length(a_x)-1
    %    compute velocities and delta positions in body frame
    %    initialize alpha angle and positions in inertial
    %    clean apparent centrifugal from Y acceleration
    dt = epoch(i+1)-epoch(i);
  
    vx = [vx; vx(i) + a_x(i+1)*dt];
    x = [x; x(i) + vx(i)*dt + 1/2*a_x(i+1)*dt^2];
    dx = [dx; x(i+1)-x(i)];

    vy = [vy; vy(i) + a_y(i+1)*dt];
    y = [y; y(i) + vy(i)*dt + 1/2*a_y(i+1)*dt^2];
    dy = [dy; y(i+1)-y(i)];

    x_y =[x_y; x(i+1) y(i+1)];    

    centr_y = [centr_y; sqrt(vx(i)^2+vy(i)^2) * omegaz(i+1)];
    a_y(i+1) = a_y(i+1) - centr_y(i+1);
    vy(i+1) = vy(i)+a_y(i+1)*dt;
    y(i+1) = y(i) + vy(i)*dt + 1/2*a_y(i+1)*dt^2;
    dy(i+1) = y(i+1)-y(i);
 
    delta_xy=[dx dy];
    x_y(i+1, 2) = y(i+1);

    alpha = [alpha; alpha(i) + omegaz(i+1)*dt];
    R = [cos(alpha(i+1)) sin(alpha(i+1)); -sin(alpha(i+1)) cos(alpha(i+1))];

    x_y(i+1, :) = x_y(i, :) + delta_xy(i+1, :)*R';
    x(i+1) = x_y(i+1, 1);
    y(i+1) = x_y(i+1, 2);
 end
end