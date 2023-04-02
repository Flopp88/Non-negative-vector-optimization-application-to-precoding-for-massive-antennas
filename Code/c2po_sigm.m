function[x_quant,summary_X_real,summary_U, xi] = c2po_sigm(s,H,K,M,tau,tmax,delta,lambda)

ampl = abs(real(s(1)));
A = [real(H), -imag(H);imag(H), real(H)];
S_real = [diag(real(s)) diag(imag(s))];
Mvect = S_real*A;
Ik = [eye(K) 1i*eye(K)];
x = conj(H).'*s;
X_real = [real(x); imag(x)];
rho = 1/(1-delta*tau);
xi = repmat(sqrt((norm(ampl*unit_quant(x))^2)/(2*M)),2*M,1);  %on a posé P=M pour obtenir une sortie sur le cercle unité
summary_U = zeros(K, tmax-1);
summary_X_real=zeros(2*M, tmax);
summary_X_real(:,1) = X_real;
for i=2:tmax
    i;
    
    step = sigmoid(Mvect*summary_X_real(:,i-1),lambda);
    summary_U(:,i-1) = step;
    
    
    deriv = (diag(s)*diag(step)* (eye(K) + lambda*diag(Mvect*summary_X_real(:,i-1))*(eye(K) - diag(step)))...
           *Mvect - Ik*A).';
    fct = diag(s) * diag(step)*Mvect*summary_X_real(:,i-1) - Ik*A*summary_X_real(:,i-1);
    
    z = summary_X_real(:,i-1) - tau * 2 * real(deriv*conj(fct));
    %z = summary_X_real(:,i-1) - tau * 2 * real(deriv.'*conj(deriv))*summary_X_real(:,i-1);
    x_it = min_val(max_val(rho*z,-xi),xi);
    
    summary_X_real(:,i) = x_it;
end
%Quantization and mapping to complex space
x_quant = ampl*(sign(summary_X_real(1:M,end)) + 1i*sign(summary_X_real(M+1:2*M,end)));

% if b~=0
%     a
%     b
% end
%relu(Mvect*x_it)

%Mvect*x_it



end