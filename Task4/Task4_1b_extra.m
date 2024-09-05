clc 
clear all
close all

%Reference Model Parameters
a1 = 20;
b1 = 100;

a2 = 60;
b2 = 900;

%Reference input
A = 0.0873;
w = 90;

r1 = @(t) A*(b1-w^2)*sin(w*t) + A*a1*w*cos(w*t);
r2 = @(t) A*(b2-w^2)*sin(w*t) + A*a2*w*cos(w*t);
y_d = @(t) A*sin(w*t);

%Controller Parameters
l = 3;
[t11,t21,t31] = get_theta(a1,b1,l);
[t12,t22,t32] = get_theta(a2,b2,l);
c0 = 1/2;

%Simulations
[ts1,x1] = ode45(@(ts1,x1) odefun(x1,r1,ts1,t11,t21,t31,c0,l,a1,b1),[0 10*2*pi/w], [0 0 0 0 0 0]);
y1 = x1(:,1);
y_d1 = y_d(ts1);

[ts2,x2] = ode45(@(ts2,x2) odefun(x2,r2,ts2,t12,t22,t32,c0,l,a2,b2),[0 10*2*pi/w], [0 0 0 0 0 0]);
y2 = x2(:,1);
y_d2 = y_d(ts2);

%Ploting
plot(ts1,y1-y_d1,ts2,y2-y_d2)
legend('p=10','p=30')
title('Desired output error (y-y_d)')

%Functions
function [t1,t2,t3] = get_theta(a,b,l)
    t1 = 2 -a;
    t2 = 10*a + 2*l - a*l + (a*l^2)/2 - l^2 - 20;
    t3 = a - b/2 + l - (a*l)/2 + 8;
end

function dx = odefun(x,r,t,t1,t2,t3,c0,l,a,b)
    dx = zeros(6,1);

    u = t1*x(3) + t2*x(4) + t3*x(1) + c0*r(t);
    
    dx(1) = x(2);
    dx(2) = -20*sin(x(1)) - 2*x(2) + 2*u;
    dx(3) = -l*x(3) + u;
    dx(4) = -l*x(4) + x(1);
    dx(5) = x(6);
    dx(6) = -b*x(5) - a*x(6) + r(t);
end   