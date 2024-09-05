clc 
clear all
close all

%Reference Model Parameters
a = 7;
b = 12;

%Controller Parameters
l = 3; 
t1 = 2 -a;
t2 = 10*a + 2*l - a*l + (a*l^2)/2 - l^2 - 20;
t3 = a - b/2 + l - (a*l)/2 + 8;
c0 = 1/2;

%Αρχικές συνθήκες
r = @(t) 0*t;

%Ευρωστία ως προς το λ
L = [1 3 5 10 50 100];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of λ \n Output response for initial contidions [0.1745 0] \n Reference Model: s^2 + %ds + %d',a,b))

for i = L
    [t1,t2,t3] = get_theta(a,b,i);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,i,a,b),[0 10], [0.1745 0 0 0 0.1745 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);
    
    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y,t,y_m)
    hold on
    yline(0,'--')
    hold off
    legend('Closed Loop System','Reference Model')
    title(sprintf('λ = %d\nMean Absolute Error %f',i,mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of λ \n Output response for initial contidions [0.8727 0] \n Reference Model: s^2 + %ds + %d',a,b))

for i = L
    [t1,t2,t3] = get_theta(a,b,i);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,i,a,b),[0 10], [0.8727 0 0 0 0.8727 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y,t,y_m)
    hold on
    yline(0,'--')
    hold off
    legend('Closed Loop System','Reference Model')
    title(sprintf('λ = %d\nMean Absolute Error %f',i,mae))
end

%Ευρωστία ως προς το μοντέλο αναφοράς
p = [1 3 5 10 20 30];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output response for initial contidions [0.1745 0], λ = %d',l))

for i = p
    a = 2*i;
    b = i^2;
    
    [t1,t2,t3] = get_theta(a,b,l);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,l,a,b),[0 10], [0.1745 0 0 0 0.1745 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y,t,y_m)
    hold on
    yline(0,'--')
    hold off
    legend('Closed Loop System','Reference Model')
    title(sprintf('s^2 + %ds + %d (p = %d)\nMean Absolute Error %f',a,b,i,mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output response for initial contidions [0.8727 0], λ = %d',l))

for i = p
    a = 2*i;
    b = i^2;
    
    [t1,t2,t3] = get_theta(a,b,l);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,l,a,b),[0 10], [0.8727 0 0 0 0.8727 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);

    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y,t,y_m)
    hold on
    yline(0,'--')
    hold off
    legend('Closed Loop System','Reference Model')
    title(sprintf('s^2 + %ds + %d (p = %d)\nMean Absolute Error %f',a,b,i,mae))
end


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