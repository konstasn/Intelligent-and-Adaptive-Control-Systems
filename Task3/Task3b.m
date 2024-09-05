clc
clear all
close all

%Reference Model Parameters
a = 7;
b = 12;

%Controller parameters
g = 1e7;
G = g*eye(4); %Matix Γ
l = 10; %Filter Λ(s)
p0 = 1; %Filter P(s)

%Theta* values
t1 = 2 -a;
t2 = 10*a + 2*l - a*l + (a*l^2)/2 - l^2 - 20;
t3 = a - b/2 + l - (a*l)/2 + 8;
c0 = 1/2;

%Sinwave Values
A_values = [0.0175 0.0873];
w_values = [0.5 90];

for i = 1:2
    A = A_values(i);
    w = w_values(i);
    
    %Reference model input
    r = @(t) A*(b-w^2)*sin(w*t) + A*a*w*cos(w*t);
    
    x0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [t,x] = ode45(@(t,x) odefun(x,r,t,a,b,l,p0,G),[0 10*2*pi/w],x0);
    
    y = x(:,1);
    y_m = x(:,3);
    
    x2 = x(:,2);
    x2_m = x(:,4);
    
    theta = x(:,11:14);
    phi = x(:,7:10);
    e = y - y_m;
    
    u = theta(:,1).*x(:,5) + theta(:,2).*x(:,6) + theta(:,3).*x(:,1) + theta(:,4).*r(t) ...
        - e.*(G(1,1)*phi(:,1).^2 + G(2,2)*phi(:,2).^2 + G(3,3)*phi(:,3).^2 + G(4,4)*phi(:,4).^2);
    
    figure('Renderer','painters','Position', [100 100 1200 300])
    plot(t,y,t,y_m)
    title(sprintf('Output response tracking %.4fsin(%.1ft) \n Reference Model: s^2 + %ds + %d, l = %d, p_0 = %d, Γ = %d*I',A,w,a,b,l,p0,g))
    
    figure('Renderer','painters','Position', [100 100 1200 900])
    tL = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
    title(tL,sprintf('System response tracking %.4fsin(%.1ft) \n Reference Model: s^2 + %ds + %d, l = %d, p_0 = %d, Γ = %d*I',A,w,a,b,l,p0,g))
    
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
    plot(t,theta(:,1),t,theta(:,2),t,theta(:,3),t,theta(:,4))
    hold on
    yline(t1,'--','θ_1*')
    yline(t2,'--','θ_2*')
    yline(t3,'--','θ_3*')
    yline(c0,'--','c_0*')
    hold off
    legend('θ_1','θ_2','θ_3','c_0')
    title('Parameter convergence')
    
    nexttile
    plot(t,u)
    title('Control Signal')
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