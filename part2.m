%% Clean-up
clear all;
close all;
clc;

%% Load input signals
load('sounds.mat'); 

%% Constants values
n = length(d);      % time steps
coeff = 500;        % filter's coefficients number

%% Wiener-Hopf equations

% generate matrix U = [ u(n) u(n-1) ... u(n-k) ]'
U = zeros(coeff,n);
U( 1, 1:n ) = u';

for k = 2:coeff
        U( k, 1:k-1 ) = zeros( 1, k-1 );  %add zeros at the front of the vector
        U( k, k:n ) = u( 1:n-k+1 )';      %add u(n-k)
end

% Compute auto-correlation matrix R of u(n)
R = (1/n) * U * (U');

% Compute cross-correlation vector p between u(n) and d(n) 
p = (1/n) * U * d;

% Compute optimal Wiener-Hopf coefficients
w=R\p;

% compute desired and noise free signals
e = d-U'*w;

%% Plot resutls

figure(1);
xlabel('time steps n');

ax1 = subplot(2,1,1);  
plot(d)
title('Desired signal')
ylabel(ax1,'d(n)')

ax2 = subplot(2,1,2);  
plot(e);
title('Noise free signal');
ylabel(ax2, 'e(n)');

% soundsc(e,Fs);