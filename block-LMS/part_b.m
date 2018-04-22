close all;
clc;

n = 4096; % x lenth
m = 4096; % y length  


%% Generate random signals
x = randn(n,1) + 1i*randn(n,1);
y = randn(m,1) + 1i*randn(m,1);

%% 1. Convolution with build-in conv()
fprintf('\n\n\nConvolution with build-in conv()\n');
fprintf('======================================\n');

tic
w1=conv(x,y);
t1 = toc;
fprintf('Elapsed time: \t%f seconds\n', t1)

%% 2. Convolution with Toeplitz matrix
fprintf('\nConvolution with Toeplitz matrix\n');
fprintf('======================================\n');

tic
Y=toeplitz([y;zeros(m-1,1)],[y(1);zeros(m-1,1)]);

w2=Y*x;
t2 = toc;
fprintf('Elapsed time: \t%f seconds\n', t2)
fprintf('Error : \t%e\n', norm(w1 - w2))

%% 3. Convolution with Circulant matrix
fprintf('\nConvolution with Circulant matrix\n');
fprintf('======================================\n');

tic
yr=[y(1) zeros(1,m-1) y(end:-1:2).'];
C=gallery('circul',yr);

w3=C*[x ; zeros(n-1,1)];
t3 = toc;
fprintf('Elapsed time: \t%f seconds\n', t3)
fprintf('Error: \t\t%e\n', norm(w1 - w3))

%% 4. Convolution with FT
fprintf('\nConvolution with FT\n');
fprintf('======================================\n');

tic
yp=[y; zeros(m-1,1)];
xp=[x; zeros(n-1,1)];

w4=ifft(fft(yp).*fft(xp));
t4 = toc;
fprintf('Elapsed time: \t%f seconds\n', t4)
fprintf('Error : \t%e\n', norm(w1 - w4))