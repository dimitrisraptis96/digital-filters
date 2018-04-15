%% Clean-up
clear all;
close all;
clc;

%% Constants values
n = 1000;               % time steps
coeff = 2;              % filter's coefficients number
varV = 0.18;            % white noise's variance
mu = [0.01 1 2.5 4];	% steepest descent m parameter
epsilon = 1.0e-8;       % steppest descent epsilon parameter

%% Initiate vectors
x = zeros(n,1);     % input signal
d = zeros(n,1);     % input signal with white noise
u = zeros(n,1);     % independent noise measurement
e = zeros(n,1);     % clean signal

%% Create functions

v = sqrt(varV) .* randn(n, 1); % Gaussian white noise
v = v - mean(v);

for i=1:n     
    x(i) = cos(pi*i) * sin( (pi/25)*i + pi/3 );
    d(i) = x(i) + v(i);
    
    if i==1;    u(1)=v(1);
    else;       u(i) = -0.78 * u(i-1) + v(i);
    end
    
end

%% Wiener-Hopf equations

% Compute autocorrelation matrix R of u(n) 
u_minus_1 = [0; u(1:n-1)];

r = [u'; u_minus_1'];
R = (1/n)*r*(r');

% Compute cross-correlation vector p between u(n) and d(n) signals
p0 = mean(u.*d);
p1 = mean(u_minus_1.*d);
p = [p0 ; p1];


% Compute optimal Wiener-Hopf coefficients
w0 = R \ p;
fprintf('Optimal Wiener-Hopf coefficients: \t w0 = %f \t and \t w1 = %f\n', w0(1), w0(2));



%% Steepest Descent algorithm

% Calculate the range of the coefficient μ 
min_m = 0;
max_m = 2 / max(eig(R));
fprintf('\nRange of coefficient μ: \t\t %f < μ < %f\n', min_m, max_m);


% Apply SD aalgorithm for varing mu parameter
for i=1:length(mu)
    [wt,isConverged] = steepest_descent(mu(i),epsilon,n,R,p);

    fprintf('\nSD algorithm coefficients: \t\t w0 = %f \t and \t w1 = %f\n', wt(1,n), wt(2,n));

    figure(i)
    we = (wt - w0*ones(1,n)).^2;
    e = sqrt(sum(we));

    semilogy(e);
    xlabel('time steps n');
    ylabel('Parameter error');
    title('Parameter error');
    
    %% contour curves and trajectories
    
    L = 50;
    ww = linspace(-2.5,2.5,L);

    J = zeros([L,L]);
    sigma2d = 0.1;

    % Construct the error surface
    for j=1:L
      for k=1:L
        wp = [ww(j); ww(k)];
        J(k,j) = sigma2d - 2*p'*wp + wp'*R*wp;
      end
    end

    min_J = min(J(:));
    max_J = max(J(:));

    levels = linspace(min_J,max_J,12);

    figure(5+i)
    contourf(ww, ww, J, levels); axis square
    hold on

    plot(wt(1,:), wt(2,:), '.r-');
    plot(w0(1), w0(2), 'ob');
    hold off
    colorbar
    xlabel('w(1)');
    ylabel('w(2)');
    title('Error Surface and Adaptation process');
    %%
end

y = r'*wt;

e = d - y;

% J = mean((d-y).^2);
% fprintf('Mean-square error computed: %f\n', J);

%%

%%

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

% figure(2)
% plot([d y])
% legend({'d(n)', 'y(n)'})

%% parameter error
% figure(3)
% we = (wt - w0*ones(1,n)).^2;
% e = sqrt(sum(we));
% 
% semilogy(e);
% xlabel('time steps n');
% ylabel('Parameter error');
% title('Parameter error');

% figure(4)
% plot([d-y x])
% legend({'e', 'x'})
% figure
% subplot(2,1,1)  
% plot(v)
% title('v(n)')
% 
% subplot(2,1,2)  
% plot(u)
% title('u(n)')

%% Function declaration

function [Wt,isConverged] = steepest_descent(mu,epsilon,n,R,p)
    W = [0 ; 0]; 
    Wt = zeros(2,n);
    Wt(:,1) = W;

    for k=2:n
         W = W + mu * (p - R*W);
         Wt(:,k) = W;
    end
    
    if norm(mu * (p - R*W),'fro') > epsilon 
        fprintf('\nNot converged\n');
        isConverged=true;
    else 
        fprintf('\nConverged\n');
        isConverged=false;
    end

end
