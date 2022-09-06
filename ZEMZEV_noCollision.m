clc;clear;close all

%% LANDER IC %%
r0 = [-2 0 1.5]'*1000;
v0 = [0 0 -75]';
m_wet = 1905;
t0 = 0;

%% LANDER TC %%
rd = [0 0 0]';
vd = [0 0 0]';
tf = 100;

%% THRUST PARAMS %%
theta = 27*pi/180;
n = 6;
Isp = 226;
Tmax = 3100;
Fmin = 0.3*Tmax;
Fmax = 0.8*Tmax;

%% PLANET PARAMS %%
g = [0 0 -3.7114]';
ge = 9.807;

%% DE solve %%
N = 10000;
tSpan = linspace(t0,tf,N);
init = [r0' v0' m_wet']';

[t,x] = ode45(@odefun, tSpan, init);

%% PLOTS %%

% TRAJECTORY %
figure(1)
plot3(x(:,1),x(:,2),x(:,3), 'LineWidth', 1.5);
hold on
plot3(x(1,1),x(1,2),x(1,3),'ks')
hold on
plot3(x(N,1),x(N,2),x(N,3),'ro')
hold off
grid on

% POSITION %
figure(2)
subplot(3,1,1)
plot(t,x(:,1))
grid on

subplot(3,1,2)
plot(t,x(:,2))
grid on

subplot(3,1,3)
plot(t,x(:,3))
grid on

% VELOCITY %
figure(3)
subplot(3,1,1)
plot(t,x(:,4))
grid on

subplot(3,1,2)
plot(t,x(:,5))
grid on

subplot(3,1,3)
plot(t,x(:,6))
grid on

% MASS %
figure(4)
plot(t,x(:,7))
grid on

% THRUST %
A = [0 0 1]';
del = 1;
phi = del^2/3;
c = 500;
tgo = tf - t;
r = [x(:,1) x(:,2) x(:,3)];
v = [x(:,4) x(:,5) x(:,6)];
ZEM = rd' - (r + v.*tgo + 0.5*g'.*tgo.^2);
ZEV = vd' - (v + g'.*tgo);
a_av = c*(x(:,3).^2 - phi).*(tgo.^2)./(24*(x(:,3).^2 + phi).^2);
T = ((6*ZEM./tgo.^2) - (2*ZEV./tgo) + a_av).*x(:,7);
% 
% for j = 1:1:length(T)
%     for i=1:1:3
%         if abs(T(j,i)) <= 0.3*Tmax
%             T(j,i) = 0.3*Tmax;
%         elseif abs(T(j,i)) >= 0.8*Tmax
%             T(j,i) = 0.8*Tmax;
%         else
%         end
%     end
% end


% a = ((6*ZEM./tgo.^2) - (2*ZEV./tgo) + a_av);
% 
% for i = 1:1:length(a)
%     if norm(a(i,:)) <= Fmin/x(i,7)
%         a(i,:) = Fmin*a(i,:)./(x(i,7).*norm(a(i,:)));
%     elseif (norm(a(i,:)) > Fmin/x(i,7)) || (norm(a(i,:)) < Fmax/x(i,7))
%         a(i,:) = abs(a(i,:));
%     elseif norm(a(i,:)) >= Fmax/x(i,7)
%         a(i,:) = Fmax*a(i,:)./(x(i,7)*norm(a(i,:)));
%     else
%     end
%     T(i,:) = a(i,:).*x(i,7);
% end



figure(5)
subplot(3,1,1)
plot(t,T(:,1))
grid on

subplot(3,1,2)
plot(t,T(:,2))
grid on

subplot(3,1,3)
plot(t,T(:,3))
grid on

%% DE FUNC %%
function dx = odefun(t,x)
    dx = zeros(7,1);
    
    %%% PLANET PARAMS %%%
    g = [0 0 -3.7114]';
    ge = 9.807;
    
    %%% LANDER IC %%%
    r0 = [-2 0 1.5]'*1000;
    v0 = [0 0 -75]';
    
    %%% LANDER TC %%%
    rd = [0 0 0]';
    vd = [0 0 0]';
    tf = 100;
    tgo = tf - t;
    
    %%% THRUST PARAMS %%%
    theta = 27*pi/180;
    n = 6;
    Isp = 226;
    Tmax = 3100;
    Fmin = 0.3*Tmax;
    Fmax = 0.8*Tmax;
    
    %%% SAFETY REGION %%%
    A = [0 0 1]';
    del = 1;
    phi = del^2/3;
    c = 500;
    
    %%% ZEM/ZEV %%%
    r = [x(1) x(2) x(3)]';
    v = [x(4) x(5) x(6)]';
    
    ZEM = rd - (r + v*tgo + 0.5*g*tgo^2);
    ZEV = vd - (v + g*tgo);
    
    %%% THRUST GENERATION %%%
    
    a_av = A*c*(r(3)^2 - phi)*(tgo^2)/(24*(r(3)^2 + phi)^2);
    a = ((6*ZEM/tgo^2) - (2*ZEV/tgo) + a_av);
  
%     if norm(a) <= Fmin/x(7)
%         a = Fmin*a/(x(7)*norm(a));
%     elseif (norm(a) > Fmin/x(7)) || (norm(a) < Fmax/x(7))
%         a = abs(a);
%     elseif norm(a) >= Fmax/x(7)
%         a = Fmax*a/(x(7)*norm(a));
%     else
%     end
    
    T = a*x(7);
    
%     for i=1:1:3
%         if abs(T(i)) == 0
%             T(i) = 0.3*Tmax;
%         elseif (abs(T(i)) <= 0.3*Tmax) && ~(abs(T(i)) == 0)
%             T(i) = 0.3*Tmax;
%         elseif abs(T(i)) >= 0.8*Tmax
%             T(i) = 0.8*Tmax;
%         else
%         end
%     end
    
    %%% ATM PERTURB %%%
%     ap = 0.2*(T./x(7))*sin(pi*t/3);
    

    %%% ODE EQNS %%%
    dx(1) = x(4);%x
    dx(2) = x(5);%y
    dx(3) = x(6);%z
    dx(4) = g(1)  + T(1)/x(7);%vx
    dx(5) = g(2)  + T(2)/x(7);%vy
    dx(6) = g(3)  + T(3)/x(7);%vz
    dx(7) = -norm(T)/(Isp*norm(g));%m

end
