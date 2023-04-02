function[x_quant,summary_X, xi]=original_c2po(s,H,K,M,tau,tmax)
ampl = abs(real(s(1)));
A = (eye(K)-s*conj(s).'/(norm(s)^2))*H;
x = conj(H).'*s;
rho = 1/(1-tau);
xi = repmat(sqrt((norm(ampl*unit_quant(x))^2)/(2*M)),2*M,1); %on a posé P=M pour obtenir une sortie sur le cercle unité
 
summary_X=zeros(M, tmax);
summary_X(:,1) = x;
for i=2:tmax
    %i
    z = summary_X(:,i-1) - tau * (conj(A).')*A*summary_X(:,i-1);
    x_it = min_val(max_val(rho*real(z),-xi),xi)+ 1j*min_val(max_val(rho*imag(z),-xi),xi);
    
    summary_X(:,i) = x_it;
end
%Quantization and mapping to complex space
x_quant = xi(1)*(sign(real(summary_X(1:M,end))) + 1i*sign(imag(summary_X(1:M,end))));

sH = conj(s).';
alpha = (sH*H*summary_X(:,end))/(norm(s)^2);
if -pi/4>=angle(alpha) || angle(alpha)>=pi/4 || abs(alpha)==0
    disp(alpha);
end 
    
%H*summary_X(:,end) - alpha*s
end