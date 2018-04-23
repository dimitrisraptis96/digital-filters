%% Clean-up
clear all;
close all;
clc;

k = 14;

time = zeros(k,2);
error = zeros(k,1);

for i=1:k
    % Define n and m  (m <= n)
    n=2^i;
    m=2^(i-1);
    
    % Define input vectors
    u = randn (n, 1);
    w = randn (m, 1);
    
    a =  u(m:n);
    b =  u(m:-1:1);
    
    % Toeplitz matrix-vector Fast Multiplication nlog(n) complexity
    tic
    
    n = length(a);
    c = [a; 0; fliplr(b(2:end).').'];
    p = ifft( fft(c) .* fft([w; zeros(n,1)]));
    yf = p(1:n);
    
    time(i,1) = toc;

    
    % Product n^2 complexity
    tic
    
    T = toeplitz ( a, b); % large memory
    y = T * w;
    
    time(i,2) = toc;
    
    
    error(i) = norm ( yf - y ) / norm(y);        
end

fprintf('Error: '); display(error);
fprintf('Time: '); display(time');

figure(1);
plot(time);
xlabel( 'Time step 2^n' );  
ylabel( 'Time elapsed (seconds)' );
title( 'y =Tw multiplication' );
legend( { 'n log(n)', 'n^2' } )

figure(2);
plot(error);
xlabel( 'Time step 2^n' );  
ylabel( ' norm( yf - y ) / norm(y)' );
title( 'Relative Error' );