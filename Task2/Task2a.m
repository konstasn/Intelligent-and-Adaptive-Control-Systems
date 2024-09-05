clc
clear all
close all

%Reference Model Parameters
a = 7;
b = 12;

A_m = [0 1;-b -a];
B_m = [0;1];

%Solving Lyapunov equation to calutate matrix P
g = 1e4;
Q = g*eye(2);

P = lyap(A_m',Q);

%Calculate K* & L*
M = 1/2;
G = 10;
C = 1;
A = [0 1;-G/M -C/M];
B = [0;1/M];

K_star = 1/norm(B)^2 * [b/M-G/M^2 a/M-C/M^2];
L_star = 1/(M*norm(B)^2);

%Αρχικές Συνθήκες
q0 = [0.175 0.8275];
r = @(t) 0*t;

for i = q0
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,P),[0 10], [i 0 i 0 0 0 0]);
    
    y = x(:,1);
    y_m = x(:,3);
    
    x2 = x(:,2);
    x2_m = x(:,4);
    
    K = x(:,5:6);
    L = x(:,7);
    
    u = -(x(:,5).*x(:,1) + x(:,6).*x(:,2)) + x(:,7).*r(t);
    
    figure('Renderer','painters','Position', [100 100 1200 300])
    plot(t,y,t,y_m)
    legend('Closed Loop System','Reference Model')
    title(sprintf('Output response for initial contidions [%.4f 0] \n Reference Model: s^2 + %ds + %d, Q = %d*I',i,a,b,g))
    
    figure('Renderer','painters','Position', [100 100 1200 900])
    tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
    title(tL,sprintf('System response for initial contidions [%.4f 0] \n Reference Model: s^2 + %ds + %d, Q = %d*I',i,a,b,g))
    
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
    plot(t,x2-x2_m)
    title('Speed Error')
    
    nexttile
    plot(t,K(:,1),t,K(:,2),t,L(:))
    hold on
    yline(K_star(1),'--','k_1*')
    yline(K_star(2),'--','k_2*')
    yline(L_star,'--','L*')
    hold off
    title('Parameter convergence')
    
    nexttile
    plot(t,u)
    title('Control Signal')
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