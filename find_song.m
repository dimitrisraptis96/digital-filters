%% Clean-up
clear all;
close all;
clc;

%% Constants values
n = 320000;           % time steps
coeff = 500;

load('sounds.mat');                                                                                                                                                                                        


%% Wiener-Hopf equations

% generate matrix U = [ u(n) u(n-1) ... u(n-k) ]'
U = zeros(coeff,n);
U( 1, 1:n ) = u';

for k = 2:coeff
        U( k, 1:k-1 ) = zeros( 1, k-1 );  %add zeros at the front of the vector
        U( k, k:n ) = u( 1:n-k+1 )';      %add u(n-k)
end

% Compute auto-correlation matrix R of u(n)
R = (1/n) * (U) * (U');

% Compute cross-correlation vector p between u(n) and d(n) 
p = (1/n) * U * d;

% Compute optimal Wiener-Hopf coefficients
w=R\p;

% compute desired and noise free signals
y = U'*w;
e = d-y;

%% Plot resutls

figure(1);
subplot(2,1,1);  
plot(d);
title('d(n)');


subplot(2,1,2);  
plote(e);
title('e(n)');
