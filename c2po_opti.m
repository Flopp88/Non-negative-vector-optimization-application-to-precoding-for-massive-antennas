function[x_quant,summary_X_real,summary_U, xi]=c2po_opti(s,H,K,M,tau,tmax,delta)
[r_quant,x_quant,x] = zeroforcing_quant(s,H);
if r_quant == unit_quant(s)
    summary_X_real = 0;
    summary_U = 0;
    xi = 0;
    return
end

ampl = abs(real(s(1)));
A = [real(H), -imag(H);imag(H), real(H)];
S_real = [diag(real(s)) diag(imag(s))];
Mvect = S_real*A;
Ik = [eye(K) 1i*eye(K)];
x = conj(H).'*s;
% x = x_quant;
X_real = [real(x); imag(x)];
rho = 1/(1-delta*tau);
xi = repmat(sqrt((norm(ampl*unit_quant(x))^2)/(M)),2*M,1);  %on a posé P=M pour obtenir une sortie sur le cercle unité
summary_U = zeros(K, tmax-1);
summary_X_real=zeros(2*M, tmax);
summary_X_real(:,1) = X_real;
a=0;
b=0;
for i=2:tmax
    i;
    step = heaviside(Mvect*summary_X_real(:,i-1));
    U = repmat(step,1,2*M);
    summary_U(:,i-1) = U(:,1);
    
    for k=1:length(step)
       if step(k)==0
            a=a+1;
       end
    end
    if a~=0
        b=b+1;
    end
    
    deriv = diag(s)*(Mvect.*U)-Ik*A;
    z = summary_X_real(:,i-1) - tau * 2 * real(deriv.'*conj(diag(s)*...
        ((Mvect*summary_X_real(:,i-1)).*U(:,1))-Ik*A*summary_X_real(:,i-1)));
    %z = summary_X_real(:,i-1) - tau * 2 * real(deriv.'*conj(deriv))*summary_X_real(:,i-1);
    x_it = min_val(max_val(rho*z,-xi),xi);
    
    summary_X_real(:,i) = x_it;
end
%Quantization and mapping to complex space
x_quant_c2po = xi(1)*(sign(summary_X_real(1:M,end)) + 1i*sign(summary_X_real(M+1:2*M,end)));
r_c2po = unit_quant(H*x_quant_c2po);

if SER(s,r_c2po)<SER(s,r_quant)
    x_quant = x_quant_c2po;
end
% if b~=0
%     a
%     b
% end
%relu(Mvect*x_it)
end