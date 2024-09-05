clc 
clear all
close all

%Reference Model Parameters
a = 7;
b = 12;

%Controller Parameters
l = 2; 
t1 = 2 -a;
t2 = 10*a + 2*l - a*l + (a*l^2)/2 - l^2 - 20;
t3 = a - b/2 + l - (a*l)/2 + 8;
c0 = 1/2;

%Sinwave parameter values
A_values = [0.0175 0.0873];
w_values = [0.5 90];

for i = 1:2
    A = A_values(i);
    w = w_values(i);
    
    %Reference Model Input
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    %Desired Output
    y_d = @(t) A*sin(w*t);

    [t,x_c] = ode45(@(t,x_c) odefun(x_c,r,t,t1,t2,t3,c0,l,a,b),[0 10*2*pi/w], [0 0 0 0 0 0]);
    
    y = x_c(:,1);
    y_m = x_c(:,5);
    
    x2 = x_c(:,2);
    x2_m = x_c(:,6);
    
    u = t1*x_c(:,3) + t2*x_c(:,4) + t3*y + c0*r(t);
    
    figure('Renderer','painters','Position', [300 300 1200 600])
    tL = tiledlayout(2,1,'TileSpacing','Compact','Padding','compact');
    title(tL,sprintf('System response tracking %.4fsin(%.1ft) \n Reference Model: s^2 + %ds + %d, λ = %d',A,w,a,b,l))
  
    nexttile
    plot(t,y_m,t,y_d(t))
    legend('Reference Model','Desired Output')
    title('Output response')
    
    nexttile
    plot(t,y_m - y_d(t))
    title('Model Output-Desired Output Error')
    
    figure('Renderer','painters','Position', [300 300 1200 600])
    plot(t,y-y_d(t))
    title(sprintf('Desired output error tracking %.4fsin(%.1ft) \n Reference Model: s^2 + %ds + %d, λ = %d',A,w,a,b,l))
    
    figure('Renderer','painters','Position', [100 100 1200 600])
    tL = tiledlayout(2,2,'TileSpacing','Compact','Padding','compact');
    title(tL,sprintf('System response tracking %.4fsin(%.1ft) \n Reference Model: s^2 + %ds + %d, λ = %d',A,w,a,b,l))
    
    nexttile
    plot(t,y,t,y_m)
    legend('Closed Loop System','Reference Model')
    title('Output response')
    
    nexttile
    plot(t,x2,t,x2_m)
    legend('Closed Loop System','Reference Model')
    title('x_2 & x_{m2}')
    
    nexttile
    plot(t,y-y_m)
    title('Output Error')
    
    nexttile
    plot(t,u)
    title('Control Signal')
    
    figure('Renderer','painters','Position', [300 300 1200 600])
    tL = tiledlayout(2,1,'TileSpacing','Compact','Padding','compact');
    title(tL,sprintf('ω_1 & ω_2 for system tracking %.4fsin(%.1ft) \n Reference Model: s^2 + %ds + %d, λ = %d',A,w,a,b,l))
  
    nexttile
    plot(t,x_c(:,3))
    title('ω_1')
    
    nexttile
    plot(t,x_c(:,4))
    title('ω_2')
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