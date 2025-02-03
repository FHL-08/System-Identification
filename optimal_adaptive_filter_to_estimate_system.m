close all;

% the more samples passed, the higher the accuracy
N = 2048; % length of input/output
umin = -1;
umax = 1 ;
% Generate white noise between umax and umin as input to the system
% Could technically pass any input to the system.
u = umin + (umax-umin).*rand(N,1);

d_old = 0;
d = zeros(length(u), 1); % output of unknown system when us is passed as the input

% first order system model y[n] = a*y[n-1] + b*x[n]
b = 0.1;
a = 0.9;
for i=1:length(u)
    if i == 1
        d(i) = a*d_old + b*u(i);
    else
        d(i) = a*d(i-1) + b*u(i);
    end
end

% the higher the filter length, the higher the accuracy
M = 1024; % length of FIR (must be less than N)

A = toeplitz([u(1); zeros(M-1, 1)], u');
B = toeplitz(u, [u(1) zeros(1, M-1)]);
R_c = A*d; % cross correlation
R_a = A*B; % auto correlation


% Singular Value Decomposition to deal with ill-conditioned R_a
[left_singular, Singular, right_singular] = svd(R_a);
inv_S = diag(1./diag(Singular));
inv_R_a = right_singular*inv_S*(left_singular');

% QP to obtain optimal FIR filter that minimises mean squared error
fir = inv_R_a*R_c;
% half of the mean squared error (makes optimisation easier)
E = 0.5*(fir')*R_a*fir - (R_c')*fir + 0.5*(d')*d;
NE = (2*E)/N;
percentage_accuracy = 100*(1-NE);
percentage_accuracy = min(max(percentage_accuracy, 0), 100);

% Frequency response of FIR filter
[H, w] = freqz(fir, 1,"ctf",2048);

% Desired IIR filter order
order_num = 0; % Numerator order (i.e., number of zeros)
order_den = 1; % Denominator order (i.e., number of poles)

% Use invfreqz to approximate the IIR filter. This is what gives the
% transfer function of the system. b_iir = numerator, a_iir = denominator.
[b_iir, a_iir] = invfreqz(H, w, order_num, order_den);
sys_d = tf(b_iir, a_iir, 1); % discrete time transfer function.

% compare FIR approximation with IIR approximation
figure (1);
freqz(fir, 1, 2048); hold on;
freqz(b_iir, a_iir, 2048);
legend('FIR Filter', 'IIR Approximation');
title('Comparison of FIR and IIR Filters');