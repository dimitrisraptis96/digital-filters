%% Clean-up
clear all;
close all;
clc;

%% Constants values
n = 320000;           % time steps
varV = 0.18;          % white noise's variance

load('sounds.mat')

% soundsc(d,Fs)
%% Initiate vectors
% v = zeros(n,1);     % Gaussian white noise
% x = zeros(n,1);     % input signal
% d = zeros(n,1);     % input signal with white noise
% u = zeros(n,1);     % independent noise measurement
% e = zeros(n,1);     % clean signal
% 
% %% Create functions
% v = sqrt(varV) .* randn(n, 1);
% v = v - mean(v);
% 
% for i=1:n     
%     x(i) = cos(pi*i) * sin( (pi/25)*i + pi/3 );
%     d(i) = x(i) + v(i);
%     
%     if i==1;    u(1)=v(1);
%     else;       u(i) = -0.78 * u(i-1) + v(i);
%     end
%     
% %     e(i) = d(i) - u(i);
% end
% 
%% Compute autocorrelation matrix R of u(n) signal

tmp = zeros(320000,500);
tmp(1:n,1) = u';
for i = 2:500
%     if i==1
%         tmp(i)= [ u'];
%     else
        tmp(1:i-1,i) = ones(1,i-1)'; 
%         size (u(1:n-i+1))
        tmp(i:320000,i) = u(1:n-i+1)';
        
%         tmp(:,i) = [zeros(1,i-1)' u(1:n-i)']; 
%     end
end

R = (1/n) *(tmp')*(tmp);
%     
% r = [u'; 0 u(1:n-1)'];
% R = (1/n)*r*(r');

%% Compute cross-correlation vector p between u(n) and d(n) signals
p0 = mean(u .* d);
% p1 = 0;
% P = [p0 ; p1];
P = [p0 zeros(1,499)]';

w=R\P;

%% Compute optimal Wiener-Hopf coefficients
% w = R \ P;
% fprintf('Optimal Wiener-Hopf coefficients: \t w0 = %f \t and \t w1 = %f\n', w(1), w(2));

%% Calculate the range of the coefficient μ 
min_m = 0;
max_m = 2 / max(eig(R));
fprintf('\nRange of coefficient μ: \t\t\t %f < μ < %f\n', min_m, max_m);


%% Steepest descent method
    
wsm = 50*ones(500,1);
w_old = zeros(500,1);
m = 1;

while max(abs(wsm - w_old)) > 1.0e-9
    w_old=wsm;
    wsm = wsm + m*(P-R*wsm);
end
% T = [u [0; u(1:n-1)]]; 
y = tmp*wsm;

e = d - y;
% 
% % J = mean((d-y).^2);
% % fprintf('Mean-square error computed: %f\n', J);
% 
% %% Plot functions
% % figure
% % subplot(5,1,1)  
% % plot(v)
% % title('v(n)')
% % 
% % subplot(5,1,2)  
% % plot(x)
% % title('x(n)')
% % 
% % subplot(5,1,3)  
% % plot(d)
% % title('d(n)')
% % 
% % subplot(5,1,4)  
% % plot(u)
% % title('u(n)')
% % 
% % subplot(5,1,5)  
% % plot(e)
% % title('e(n)')
% 
figure
subplot(2,1,1)  
plot(d)
title('d(n)')

subplot(2,1,2)  
plot(e)
title('e(n)')
