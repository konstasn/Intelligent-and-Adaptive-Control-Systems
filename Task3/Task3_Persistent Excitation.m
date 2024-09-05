clc
clear all
close all

%Reference Model Parameters
a = 11;
b = 30;

%Controller parameters
g = 1e10;
G = g*eye(4); %Matix Γ
l = 100; %Filter Λ(s)
p0 = 1; %Filter P(s)

%Theta* values
t1 = 2 -a;
t2 = 10*a + 2*l - a*l + (a*l^2)/2 - l^2 - 20;
t3 = a - b/2 + l - (a*l)/2 + 8;
c0 = 1/2;

%Reference model input
r = @(t) 0.05*sin(5*t) + 0.05*sin(7*t);

x0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[t,x] = ode45(@(t,x) odefun(x,r,t,a,b,l,p0,G),[0 100],x0);

y = x(:,1);
y_m = x(:,3);

x2 = x(:,2);
x2_m = x(:,4);

theta = x(:,11:14);

figure('Renderer','painters','Position', [100 100 1200 900])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('System response for r = sin(5t) + sin(7t) \n Reference Model: s^2 + %ds + %d, l = %d, p_0 = %d, Γ = %d*I',a,b,l,p0,g))

nexttile
plot(t,y-y_m)
title('Output Error')

nexttile
plot(t,x2-x2_m)
title('x_2 Error')

nexttile
plot(t,theta(:,1))
hold on
yline(t1,'--','θ_1*')
hold off
title('θ_1')

nexttile
plot(t,theta(:,2))
hold on
yline(t2,'--','θ_2*')
hold off
title('θ_2')

nexttile
plot(t,theta(:,3))
hold on
yline(t3,'--','θ_3*')
hold off
title('θ_3')

nexttile
plot(t,theta(:,4))
hold on
yline(c0,'--','c_0*')
hold off
title('c_0')
    


function dx = odefun(x,r,t,a,b,l,p0,G)
    dx = zeros(14,1);
    
    omega = [x(5);x(6);x(1);r(t)];
    phi = x(7:10);
    theta = x(11:14);
    e = x(1) - x(3);
    
    u = theta' * omega - e * phi' * G  * phi;
    
    dx(1) = x(2);                           %x_1
    dx(2) = -20*sin(x(1)) - 2*x(2) + 2*u;   %x_2
    dx(3) = x(4);                           %x_1m
    dx(4) = -b*x(3) - a*x(4) + r(t);        %x_2m
    dx(5) = -l*x(5) +  u;                   %ω_1
    dx(6) = -l*x(6) + x(1);                 %ω_2
    dx(7:10) = -p0*phi + omega;             %φ
    dx(11:14) = -e*G*phi;                   %θ

end