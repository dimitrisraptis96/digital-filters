%% Adaptive Digital Filters - Homework 1
%% Method of Steepest Descend
%% Author: Ioannis Antoniadis 7137, April 2014 

clear

%% Section 1 - Definition of input signals
% Number of discrete time values
n = 200;

% Input signal x(n)
x = zeros(n,1);
% Initial phase of x
phase = pi/2; 
for i=1:1:n
    x(i) = sin(i*pi/4+phase) + cos(i*pi/2+phase);
end

% White noise signal v(n)
% Variance of v(n) = 0.34
% Construction of v(n)
v = rand(n,1);
v = v - mean(v);
v = 0.34*v/norm(v);

% Filter input u(n) as  AR process
u = zeros(n,1);
u(1) = v(1);
for i=2:1:n
    u(i) = -0.25*u(i-1) + v(i);
end

% Desirable output d(n)
d = x + v;

figure(1)
plot([d x u])
legend({'d(n)', 'x(n)', 'u(n)'})

%% Section 2 - Optimum Wiener Solution
% autocorrelation E[u u']
R = [0.3626 -0.09; -0.09 0.3626];
lamda = eig(R);

% cross correlation E[u d]
p = [0.34; 0];

% optimum tap weights
wo = R \ p;                         

%% Section 3 - Steepest Descend Algorithm
% initialization of tap weights
w = [0; 0];
% convergence parameter
mu = 4;           

% Tap weights matrix for all n discrete moments
wt = zeros([2,n]); 
wt(:,1) = w;
% Output initialization
y = zeros(n, 1);

% Steepest Descend loop
s = [0; u];
for i=2:n
  w = w + mu*(p-R*w);
  wt(:,i) = w;
  y(i) = s(i:-1:i-1)' * w;
end

y(1)=[];
y(n)=0;

% Noise free final output
xf = d - y;

% Plot the original x(n) and the filtered output xf(n)
figure(2)
plot([x xf])
legend({'x(n)', 'xf(n)'})

%% Section 4 - Error parameters
figure(3)
we = (wt - wo*ones(1,n)).^2;
e = sqrt(sum(we));

semilogy(e);
xlabel('time step n');
ylabel('Parameter error');
title('Parameter error');

%% Section 5 - Contour Curves and Trajectories
L = 50;
ww = linspace(-2.5,2.5,L);

J = zeros([L,L]);
sigma2d = 0.1;

% Construct the error surface
for i=1:L
  for k=1:L
    wp = [ww(i); ww(k)];
    J(k,i) = sigma2d - 2*p'*wp + wp'*R*wp;
  end
end

min_J = min(J(:));
max_J = max(J(:));

levels = linspace(min_J,max_J,12);

figure(4)
contourf(ww, ww, J, levels); axis square
hold on

plot(wt(1,:), wt(2,:), 'xr--');
hold off
colorbar
xlabel('w(1)');
ylabel('w(2)');
title('Error Surface and Adaptation process');
