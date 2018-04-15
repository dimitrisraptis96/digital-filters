%% Clean-up
clear all;
close all;
clc;

%% Constants values
n = 1000;           % time steps
varV = 0.18;        % white noise's variance
m = 2;              % steepest descent m parameter
epsilon = 1.0e-8;   % steppest descent epsilon parameter

% load('sounds.mat')
%% Initiate vectors
v = zeros(n,1);     % Gaussian white noise
x = zeros(n,1);     % input signal
d = zeros(n,1);     % input signal with white noise
u = zeros(n,1);     % independent noise measurement
e = zeros(n,1);     % clean signal

%% Create functions
v = sqrt(varV) .* randn(n, 1);
v = v - mean(v);

for i=1:n     
    x(i) = cos(pi*i) * sin( (pi/25)*i + pi/3 );
    d(i) = x(i) + v(i);
    
    if i==1;    u(1)=v(1);
    else;       u(i) = -0.78 * u(i-1) + v(i);
    end
    
%     e(i) = d(i) - u(i);
end

%% Compute autocorrelation matrix R of u(n) signal
u_minus_1 = [0; u(1:n-1)];

r = [u'; u_minus_1'];
R = (1/n)*r*(r');

%% Compute cross-correlation vector p between u(n) and d(n) signals
p0 = mean(u.*d);
p1 = mean(u_minus_1.*d);
p = [p0 ; p1];


%% Compute optimal Wiener-Hopf coefficients
w0 = R \ p;
fprintf('Optimal Wiener-Hopf coefficients: \t w0 = %f \t and \t w1 = %f\n', w0(1), w0(2));

%% Calculate the range of the coefficient μ 
min_m = 0;
max_m = 2 / max(eig(R));
fprintf('\nRange of coefficient μ: \t\t\t %f < μ < %f\n', min_m, max_m);


%% Steepest descent method

[w,error,steps,wt] = steepest_descent(m,epsilon,n,R,p);
% w = [0 ; 0];
% w_prev = [0 ; 0];
% m = 2;
% error = 1.0e-8;
% steps = 1;
% 
% while norm(P - R * w, 'fro') > error
%     w_prev=w;
%     w = w + m*(P - R * w_prev);
%     steps = steps + 1;
% end
% 
% steps
% % T = [u [0; u(1:n-1)]]; 
% w0
% w
y = r'*w;

e = d - y;

% J = mean((d-y).^2);
% fprintf('Mean-square error computed: %f\n', J);

%% Plot functions
% figure
% subplot(5,1,1)  
% plot(v)
% title('v(n)')
% 
% subplot(5,1,2)  
% plot(x)
% title('x(n)')
% 
% subplot(5,1,3)  
% plot(d)
% title('d(n)')
% 
% subplot(5,1,4)  
% plot(u)
% title('u(n)')
% 
% subplot(5,1,5)  
% plot(e)
% title('e(n)')

figure(2)
plot([d y])
legend({'d(n)', 'y(n)'})

%% parameter error
figure(3)
we = (wt - w0*ones(1,n)).^2;
e = sqrt(sum(we));

semilogy(e);
xlabel('time step n');
ylabel('Parameter error');
title('Parameter error');

figure(4)
plot([d-y x])
legend({'e', 'x'})
% figure
% subplot(2,1,1)  
% plot(v)
% title('v(n)')
% 
% subplot(2,1,2)  
% plot(u)
% title('u(n)')

%% Function declaration

function [w,error,steps,wt] = steepest_descent(m,epsilon,n_max,R,p)
    w = [0 ; 0];
    error = zeros(n_max,1);
    wt = zeros([2,n_max]); wt(:,1) = w;
    
    steps = 1;
    error(1) = norm(p - R * w, 'fro');
    
    while error(steps) > epsilon && steps < n_max
        w_prev=w;
        w = w + m * (p - R * w_prev);
        wt(:,steps) = w;
        
        steps = steps + 1;              
        error(steps) = norm(p - R * w, 'fro');   %calculate new error
    end    

    if steps >= n_max; fprintf('\nNot converge'); end;

end
