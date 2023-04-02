
%clear all
close all

K = 3;
M= 10;

H = wgn(K,M,0)+1i*wgn(K,M,0);
s = (randsrc(K,1)+1i*randsrc(K,1))/sqrt(2);

A = [real(H) -imag(H);imag(H) real(H)];

S_real = [diag(real(s)) diag(imag(s))];
Mvect = S_real*A;
Ik = [eye(K) 1i*eye(K)];
% Ik_inv_left = [eye(K) ;eye(K)/1i]/2;
% Ik_inv_right = [eye(K) zeros(K) ;zeros(K) eye(K)/1i];
% IM_inv_left = [eye(M) ;eye(M)/1i]/2;

rad = rand(M,1)*2*pi;
x = cos(rad) + 1i*sin(rad);
X_real = [real(x); imag(x)];

%x = conj(H).'*s;
%X_real = [real(x); imag(x)];
plus = relu(Mvect*X_real);

tau = 0.00001;
tmax = 10000;
[x_quant1, X_out,summary_U, xi] = c2po(s, H ,K ,M ,tau ,tmax,1);

% disp(sum(summary_U(:)==1))
% disp(K*(tmax-1))

X_final = X_out(:,end);
x_out = X_final(1:M)+1i*X_final(M+1:2*M);

x_quant = xi(1)*(sign(X_out(1:M,end)) + 1i*sign(X_out(M+1:2*M,end)));

r = H*x_quant;

figure()
plot(s,'o')
figure()
plot(r,'o')

r_quant = xi(1)*((sign(real(r))) + 1i*sign(imag(r)));

figure()
plot(r_quant,'o')


disp(r_quant/norm(r_quant) - s/norm(s))


%% Original C2PO
tau = 0.001;
tmax = 1000;
[x_quant,X_out, xi] = original_c2po(s, H ,K ,M ,tau ,tmax);

r = H*x_quant;
r_quant = unit_quant(r);
disp(r_quant/norm(r_quant) - s/norm(s))


%%

x = bussgang(s,H);

x_quant = unit_quant(x);
r = H*x_quant;

r_quant = unit_quant(r);
s = unit_quant(s);


disp(r_quant/norm(r_quant) - s/norm(s))

s-r_quant

abs(s-r_quant)

ser = SER(s,r_quant)


%% c2po sigmoid
lambda = 1000;
tau = 0.00001;
tmax = 1000;

[x_quant, X_out,summary_U, xi] = c2po_sigm(s, H ,K ,M ,tau ,tmax,1, lambda);

r = H*x_quant;
r_quant = xi(1)*((sign(real(r))) + 1i*sign(imag(r)));

disp(r_quant/norm(r_quant) - s/norm(s))


%% C2PO shared phase
clear all
K = 5;
M= 10;

H = wgn(K,M,0)+1i*wgn(K,M,0);
s = (randsrc(K,1)+1i*randsrc(K,1))/sqrt(2);

tau = 0.001;
tmax = 10;
delta = 1;
[x_quant, X_out_real, X_out_imag, xi] = c2po_shared_phase(s, H ,K ,M ,tau ,tmax,delta);

r = H*x_quant;
r_quant = unit_quant(r);

disp(r_quant/norm(r_quant) - s/norm(s))
