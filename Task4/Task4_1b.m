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

%Ευρωστία ως προς το λ
L = [1 3 5 10 50 100];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of λ \n Output error tracking 0.0175sin(0.5t) \n Reference Model: s^2 + %ds + %d',a,b))

for i = L
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    [t1,t2,t3] = get_theta(a,b,i);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,i,a,b),[0 10*2*pi/w], [0 0 0 0 0 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('λ = %d\nMean Absolute Error %f',i,mae))
end

%Αλλαγή μοντέλου αναφοράς για την 2η περίπτωση
a = 40;
b = 400;
    
figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of λ \n Output error tracking 0.0873sin(90t) \n Reference Model: s^2 + %ds + %d',a,b))

for i = L
    [t1,t2,t3] = get_theta(a,b,i);
    
    A = 0.0873;
    w = 90;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,i,a,b),[0 10*2*pi/w], [0 0 0 0 0 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('λ = %d\nMean Absolute Error %f',i,mae))
end

%Ευρωστία ως προς το μοντέλο αναφοράς
p = [1 3 5 10 20 30];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output error tracking 0.0175sin(0.5t), λ = %d',l))

for i = p
    a = 2*i;
    b = i^2;
    
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    [t1,t2,t3] = get_theta(a,b,l);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,i,a,b),[0 10*2*pi/w], [0 0 0 0 0 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('s^2 + %ds + %d (p = %d)\nMean Absolute Error %f',a,b,i,mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output error tracking 0.0873sin(90t), λ = %d',l))

for i = p
    a = 2*i;
    b = i^2;
    
    A = 0.0873;
    w = 90;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    [t1,t2,t3] = get_theta(a,b,l);
    
    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,i,a,b),[0 10*2*pi/w], [0 0 0 0 0 0]);
    y = x_c(:,1);
    y_m = x_c(:,5);
    
    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
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