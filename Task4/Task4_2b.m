clc 
clear all
close all

%Reference Model Parameters
a = 7;
b = 12;

A_m = [0 1;-b -a];
B_m = [0;1];

%Solving Lyapunov equation to calutate matrix P
g = 1e5;
Q = g*eye(2);

P = lyap(A_m',Q);

%Ευρωστία ως προς το Q
G = [1 100 1e4 1e6 1e7 1e8];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of Q \n Output error tracking 0.0175sin(0.5t) \n Reference Model: s^2 + %ds + %d',a,b))

for i = G
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    P_i = lyap(A_m',i*eye(2));
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,P_i),[0 10*2*pi/w], [0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('Q = 10^{%d}*I\nMean Absolute Error %f',log10(i),mae))
end

%Αλλαγή μοντέλου αναφοράς
a = 40;
b = 400;

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of Q \n Output error tracking 0.0873sin(90t) \n Reference Model: s^2 + %ds + %d',a,b))

for i = G    
    A_m = [0 1;-b -a];
    P_i = lyap(A_m',i*eye(2));
    
    A = 0.0873;
    w = 90;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,P_i),[0 10*2*pi/w], [0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
    title(sprintf('Q = 10^{%d}*I\nMean Absolute Error %f',log10(i),mae))
end

%Ευρωστία ως προς το μοντέλο αναφοράς
p = [1 3 5 10 20 30];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output error tracking 0.0175sin(0.5t), Q = 10^{%d}*I',log10(g)))

for i = p
    a = 2*i;
    b = i^2;
    
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    A_m = [0 1;-b -a];
    P_i = lyap(A_m',g*eye(2));
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,P_i),[0 10*2*pi/w], [0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('s^2 + %ds + %d (p = %d)\nMean Absolute Error %f',a,b,i,mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output error tracking 0.0873sin(90t), Q = 10^{%d}*I',log10(g)))

for i = p
    a = 2*i;
    b = i^2;
    
    A = 0.0873;
    w = 90;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    A_m = [0 1;-b -a];
    P_i = lyap(A_m',g*eye(2));
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,P_i),[0 10*2*pi/w], [0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
    title(sprintf('s^2 + %ds + %d (p = %d)\nMean Absolute Error %f',a,b,i,mae))
end


function dx = odefun(x,r,t,a,b,P)
    dx = zeros(7,1);

    u = -(x(5)*x(1) + x(6)*x(2)) + x(7)*r(t);
    e1 = x(1) - x(3);
    e2 = x(2) - x(4);
    
    dx(1) = x(2);
    dx(2) = -20*sin(x(1)) - 2*x(2) + 2*u;
    dx(3) = x(4);
    dx(4) = -b*x(3) - a*x(4) + r(t);
    dx(5) = P(1,2)*e1*x(1) + P(2,2)*e2*x(1);
    dx(6) = P(1,2)*e1*x(2) + P(2,2)*e2*x(2);
    dx(7) = -(P(1,2)*e1 + P(2,2)*e2)*r(t);
end