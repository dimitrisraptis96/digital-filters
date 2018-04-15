%% Clean-up
clear all;
close all;
clc;

%% 1. Define parameters and signals

% Constants values
n = 1000;               % time steps
coeff = 2;              % filter's coefficients number
varV = 0.18;            % white noise's variance
mu = [0.01 1 2 4 10];	% steepest descent m parameter
epsilon = 1.0e-8;       % steppest descent epsilon parameter

% Initiate vectors
x = zeros(n,1);     % input signal
d = zeros(n,1);     % input signal with white noise
u = zeros(n,1);     % independent noise measurement
e = zeros(n,1);     % clean signal

% Create functions
v = sqrt(varV) .* randn(n, 1); % Gaussian white noise
v = v - mean(v);

for i=1:n     
    x(i) = cos(pi*i) * sin( (pi/25)*i + pi/3 );
    d(i) = x(i) + v(i);
    
    if i==1;    u(1)=v(1);
    else;       u(i) = -0.78 * u(i-1) + v(i);
    end
    
end

%% 2. Wiener-Hopf equations

U = [u'; [0; u(1:n-1)]'];

% autocorrelation matrix R of u(n) 
R = (1/n) * U * (U');

% cross-correlation vector p between u(n) and d(n) signals
p = (1/n) * U * d;

% optimal Wiener-Hopf coefficients
w0 = R \ p;

fprintf('Optimal Wiener-Hopf coefficients:'); display(w0);
fprintf('\n=================================================\n');


%% 3. Steepest Descent algorithm

% Calculate the range of the coefficient μ 
min_m = 0;
max_m = 2 / max(eig(R));
fprintf('\nRange of coefficient μ: \t%f < μ < %f\n', min_m, max_m);
fprintf('\n=================================================\n');


% Apply SD algorithm for varing mu parameter
for i=1:length(mu)
    fprintf('\nPARAMETER mu = %f\n',mu(i));
    % Compute Wiener coefficients
    [w, wt] = steepestDescent(mu(i),n,R,p); 

    fprintf('\nSD algorithm coefficients:'); display(w);

    % Error parameter
    we = (wt - w0*ones(1,n)).^2;
    e = sqrt(sum(we));
    
    figure(i)
    subplot(1,2,1);
    semilogy(e);
    xlabel('time steps n');
    ylabel('Parameter error');
    title('Parameter error');
    
    % Contour curves and trajectories
    L = 50;
    ww = linspace(-2.5,2.5,L);

    J = zeros([L,L]);
    sigma2d = var(d);

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

    subplot(1,2,2);
    contourf(ww, ww, J, levels); axis square
    hold on

    plot(wt(1,:), wt(2,:), '.r-');
    plot(w0(1), w0(2), 'ob');
    hold off
    colorbar
    xlabel('w(1)');
    ylabel('w(2)');
    title('Error Surface and Adaptation process');
    
    % Compute final noise-free signal
    y = U'*w;
    e = d-y;
    
    
    % Plot system's signals
    figure(i+length(mu));
    subplot(5,1,1); plot(v); title('v(n)');

    subplot(5,1,2); plot(x); title('x(n)');

    subplot(5,1,3); plot(d); title('d(n)');

    subplot(5,1,4); plot(u); title('u(n)');

    subplot(5,1,5); plot(e); title('e(n)');
    
    fprintf('\n=================================================\n');
end
