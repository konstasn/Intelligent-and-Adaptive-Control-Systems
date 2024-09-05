clc 
clear all
close all

%Reference Model Parameters
a = 7;
b = 12;
a2 = 40;        %Ορισμός 2 μοντέλων αναφοράς
b2 = 400;       %ένα για την κάθε περίπτωση y_d

%Controller parameters
g = 1e7;
G = g*eye(4); %Matix Γ
l = 10; %Filter Λ(s)
p0 = 1; %Filter P(s)

%Αρχικές συνθήκες
r = @(t) 0*t;

%Ευρωστία ως προς το G
gs = [1 100 1e4 1e6 1e7 1e8];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of G \n Output error tracking 0.0175sin(0.5t) \n Reference Model: s^2 + %ds + %d, λ = %d, p_0 = %d',a,b,l,p0))

for i = gs
    G_i = i*eye(4);
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,l,p0,G_i),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);
    
    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
    title(sprintf('G = 10^{%d}*I\nMean Absolute Error = %f',log10(i),mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of G \n Output error tracking 0.0873sin(90t) \n Reference Model: s^2 + %ds + %d, λ = %d, p_0 = %d',a2,b2,l,p0))

for i = gs   
    G_i = i*eye(4);
    A = 0.0873;
    w = 90;
    r = @(t) A*(b2-w^2)*sin(w*t) + A*a2*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a2,b2,l,p0,G_i),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('G = 10^{%d}*I\nMean Absolute Error = %f',log10(i),mae))
end

%Ευρωστία ως προς το λ
L = [1 3 5 10 50 100];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of λ \n Output error tracking 0.0175sin(0.5t) \n Reference Model: s^2 + %ds + %d, Γ = 10^{%d}*I, p_0 = %d',a,b,log10(g),p0))

for i = L
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,i,p0,G),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);
    
    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
    title(sprintf('λ = %d\nMean Absolute Error = %f',i,mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of λ \n Output error tracking 0.0873sin(90t) \n Reference Model: s^2 + %ds + %d, Γ = 10^{%d}*I, p_0 = %d',a2,b2,log10(g),p0))

for i = L
    A = 0.0873;
    w = 90;
    r = @(t) A*(b2-w^2)*sin(w*t) + A*a2*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a2,b2,i,p0,G),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('λ = %d\nMean Absolute Error = %f',i,mae))
end

%Ευρωστία ως προς το p0
P0 = [0.1 1 3 5 6.9 100];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of p_0 \n Output error tracking 0.0175sin(0.5t) \n Reference Model: s^2 + %ds + %d, Γ = 10^{%d}*I, λ = %d',a,b,log10(g),l))

for i = P0
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,l,i,G),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);
    
    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
    title(sprintf('p_0 = %.1f\nMean Absolute Error = %f',i,mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of p_0 \n Output error tracking 0.0873sin(90t) \n Reference Model: s^2 + %ds + %d, Γ = 10^{%d}*I, λ = %d',a2,b2,log10(g),l))

for i = P0
    A = 0.0873;
    w = 90;
    r = @(t) A*(b2-w^2)*sin(w*t) + A*a2*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a2,b2,l,i,G),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));
    
    nexttile
    plot(t,y-y_m)
    title(sprintf('p_0 = %.1f\nMean Absolute Error = %f',i,mae))
end

%Ευρωστία ως προς το μοντέλο αναφοράς
p = [1 3 5 10 20 30];

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output error tracking 0.0175sin(0.5t), Γ = 10^{%d}*I, λ = %d, p_0 = %d',log10(g),l,p0))

for i = p
    a = 2*i;
    b = i^2;
    A = 0.0175;
    w = 0.5;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,l,p0,G),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);

    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
    title(sprintf('s^2 + %ds + %d (p = %d)\nMean Absolute Error = %f',a,b,i,mae))
end

figure('Renderer','painters','Position', [10 10 1200 1200])
tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
title(tL,sprintf('Robustness in terms of reference model \n Output error tracking 0.0873sin(90t), Γ = 10^{%d}*I, λ = %d, p_0 = %d',log10(g),l,p0))

for i = p
    a = 2*i;
    b = i^2;
    A = 0.0873;
    w = 90;
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,l,p0,G),[0 10*2*pi/w], [0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
    y = x(:,1);
    y_m = x(:,3);
    
    mae = mean(abs(y-y_m));

    nexttile
    plot(t,y-y_m)
    title(sprintf('s^2 + %ds + %d (p = %d)\nMean Absolute Error = %f',a,b,i,mae))
end


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