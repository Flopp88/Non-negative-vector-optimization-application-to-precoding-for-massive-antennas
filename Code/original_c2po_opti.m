function[x_quant,phi, alpha,summary_X, xi]=original_c2po_opti(s,H,K,M,tau,tmax)
[r_quant,x_quant,x] = zeroforcing_quant(s,H);
if r_quant == unit_quant(s)
    summary_X = 0;
    xi = 0;
    phi = 0;
    alpha = 0;
    return
end

ampl = abs(real(s(1)));
A = (eye(K)-s*conj(s).'/(norm(s)^2))*H;
x = conj(H).'*s;
% x = x_quant;
% x = 2*(randi(2,size(x))-1.5)+1i*2*(randi(2,size(x))-1.5);
rho = 1/(1-tau);
xi = repmat(sqrt((norm(ampl*unit_quant(x))^2)/(M)),2*M,1); %on a posé P=M pour obtenir une sortie sur le cercle unité
 
summary_X=zeros(M, tmax);
summary_X(:,1) = x;
for i=2:tmax
    %i
    z = summary_X(:,i-1) - tau * (conj(A).')*A*summary_X(:,i-1);
    x_it = min_val(max_val(rho*real(z),-xi),xi)+ 1j*min_val(max_val(rho*imag(z),-xi),xi);
    
    summary_X(:,i) = x_it;
end
%Quantization and mapping to complex space
x_quant_c2po = xi(1)*(sign(real(summary_X(1:M,end))) + 1i*sign(imag(summary_X(1:M,end))));
r_c2po = unit_quant(H*x_quant_c2po);

if SER(s,r_c2po)<SER(s,r_quant)
    x_quant = x_quant_c2po;
end

sH = conj(s).';
alpha = (sH*H*summary_X(:,end))/(norm(s)^2);
phi = angle(alpha);
alpha = abs(alpha);
% disp("orig C2PO:");
% disp(phi);
% disp(abs(alpha));

% if -pi/4>=angle(alpha) || angle(alpha)>=pi/4 || abs(alpha)==0
%     disp(alpha);
% end 
%H*summary_X(:,end) - alpha*s
end