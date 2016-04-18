k = 10;
m = 0.01;
c = 0.2;
vi = 0;
xi = 1;
dt = 0.001;
t_steps = 500;

time = zeros(t_steps);
x = zeros(t_steps);
v = zeros(t_steps);

time(1) = 0;
x(1) = xi;
v(1) = vi;

for i=1:t_steps-1
    time(i+1) = (i+1) * dt;
    v(i+1) = v(i) - (k/m) * x(i) * dt - (c/m) * v(i) * dt;
    x(i+1) = x(i) + v(i+1) * dt;
endfor

figure(1);
plot(time,v);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs. Time');

figure(2);
plot(time,x);
xlabel('Position (m)');
ylabel('Time (s)');
title('Position vs. Time');

figure(3);
plot(x,v);
xlabel('Position (m)');
ylabel('Velocity (m/s)');
title('Phase Space');

U = 0.5*k*x.^2;
T = 0.5*m*v.^2;
E = T +U;

figure(4);
plot(time,T,time,U,time,E);
xlabel('Time (s)');
ylabel('Energy (J)');
title('Energy functions');
legend('Kinetic','Potential','Mechanical');