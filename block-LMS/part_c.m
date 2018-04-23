clear all;
close all;

%% Parameters
n = 2^14;       % time steps
T = 40  ;        % number of independent trials
L = 2^10;       % adaptive filter order
mu = 0.0004;

sigma = 0.57;
a = 0.10;

J = zeros(n, 4);

time = zeros(4,1);
%% Block LMS algorithms

for t=1:T    
    % Generate input signal for each trial
    v = zeros(n,1); v = sqrt(sigma) .* randn(n, 1); v = v - mean(v);    % gaussian white noise
    
    u = zeros(n,1); u(1) = v(1);                                        % noisy channel input                            
    for k=2:1:n         
        u(k) = -a * u(k-1) + v(k);
    end

    d = plant(u')';                                                     % desired signal
    
    
    %% 1. Two nested loops
    tic;
    
    y = zeros(n,1);
    w = zeros(L,1);
    e = zeros(n,1);
    
    
    for k=1:floor(n/L)-1
        
        phi = zeros(L,1);  
        
        for i=1:L
            y(k*L+i) = w' * u(k*L + i: -1: (k-1)*L + i + 1);

            e(k*L+i) = d(k*L+i)-y(k*L+i); 
            
            phi = phi + mu * e(k*L + i) * u(k*L + i: -1: (k-1)*L + i + 1);     
        end  
         
        w = w + phi;
    end 
    
    J(:,1) = J(:,1) + e.^2;
    
    time(1) = time(1) + toc;
    
   %% 2. One nested loop and matric calculations
    tic;
    
    y = zeros(n,1);
    w = zeros(L,1);
    e = zeros(n,1);
    
    for k=1:floor(n/L)-1        
        %Set up input signal matrix, dim. MxM 
        U = toeplitz(u(k*L:1:(k+1)*L-1),u(k*L:-1:(k-1)*L+1));

        %Set up vector with desired signal
        dvec=d(k*L:1:(k+1)*L-1);

        %calculate output signal 
        yvec=U*w;

        %calculate error vector 
        evec=dvec-yvec;

        %log error
        e(k*L:1:(k+1)*L-1)=evec;

        %calculate gradient estimate
        phi=U.'*evec;

        %update filter coefficients
        w=w+mu*phi;	
    end
    
    J(:,2) = J(:,2) + e.^2;
    
    time(2) = time(2) + toc;
    
    %% 3. Contstrained FFT
    tic;
    
    y = zeros(n,1);
    e = zeros(n,1);
    W = zeros(2*L,1);

    for k = 2*L:L:n

        U = fft( u(k-2*L+1: k));
        
        C = ifft( W .* U);
        
        y (k-L+1:k) = C(L+1:end);
        
        e (k-L+1:k) = d(k-L+1:k) - y(k-L+1:k);
        
        phi = ifft(fft([zeros(L,1); e(k-L+1:k)]) .* conj(U));
        
        PHI = fft( [ phi(1:end-L); zeros(L,1) ] );
        
        W = W + mu * PHI;
    end
    
    J(:,3) = J(:,3) + e.^2;
    
    time(3) = time(3) + toc;
    
    %% 4. Uncontstrained FFT
    tic
    
    y = zeros(n,1);
    e = zeros(n,1);
    W = zeros(2*L,1);
    
    for k = 2*L:L:n
        U = (fft(u(k-2*L+1:k)));
        
        C = ifft(W.* U);
        
        y(k-L+1:k) = C(L+1:end);
        
        e(k-L+1:k) = d(k-L+1:k) - y(k-L+1:k);
        
        phi = fft( [ zeros(L,1); e(k-L+1:k) ] ) .* conj(U);
        
        W = W + mu * phi;
    end

    J(:,4) = J(:,4) + e.^2;
    
    time(4) = time(4) + toc;
end

J = J / T;
time = time /T


%% Display BlockLMS learning curves
figure(1)
colors = ['r','g','b','y'];
for i=1:4
    subplot(4,1,i); 
    plot(J(L+1:end,i),colors(i));
    xlabel('Time step n');  
    ylabel('Ee^{2}(n)');
    title( sprintf( 'Block LMS - Learning Curves of algorithm %d', i ) );
end

figure(5)
plot(J(L+1:end,:));
xlabel('Time step n');  
ylabel('Ee^{2}(n)');
legend({'Two nested loops','One loop','Constrained FFT','Unconstrained FFT'});
title('Block LMS - Learning Curves');