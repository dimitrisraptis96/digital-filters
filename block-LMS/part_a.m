%% Clean-up
clear all;
close all;
clc;

time = zeros(16,2);
error = zeros(16,1);

for i=1:14
    
    fprintf('\n=================================================\n');
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
    F = fft(c);
    p = ifft(F.*fft([w; zeros(n,1)]));
    yf = p(1:n);
    time(i,1) = toc
    fprintf('n * log (n)');

    
    % Toeplitz matrix-vector product n^2 complexity
    
    
    tic
    T = toeplitz ( a, b); % large memory
    y = T * w;
    time(i,2) = toc
    fprintf('n^2 complexity:'); 
    
    
    error(i) = norm ( yf - y )/norm(y);
    fprintf('Error: '); display(error);
        
end

figure(1);
plot((1:16), time(:,1));
hold on
plot((1:16), time(:,2));
hold off