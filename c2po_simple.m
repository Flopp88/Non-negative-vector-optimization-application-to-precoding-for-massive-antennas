function[x_quant]=c2po_simple(s,H,K,M,tau,tmax,delta)
[r_quant,x_quant,x] = zeroforcing_quant(s,H);
if r_quant == unit_quant(s)
    summary_X_real = 0;
    summary_U = 0;
    xi = 0;
    return
end

ampl = abs(real(s(1)));
x = conj(H).'*s;
rho = 1/(1-delta*tau);
xi = repmat(sqrt((norm(ampl*unit_quant(x))^2)/(M)),2*M,1);  %on a posé P=M pour obtenir une sortie sur le cercle unité
summary_X=zeros(M, tmax);

for i=2:tmax
   x = summary_X_real(:,i-1);
   z = x - tau(conj(H).'*H*x - conj(s).'*H);
   x_it = min_val(max_val(rho*real(z),-xi),xi)+ 1j*min_val(max_val(rho*imag(z),-xi),xi);
   
   summary_X(:,i) = x_it;
end

x_quant_c2po = xi(1)*(sign(real(summary_X(1:M,end))) + 1i*sign(imag(summary_X(1:M,end))));
r_c2po = unit_quant(H*x_quant_c2po);

if SER(s,r_c2po)<SER(s,r_quant)
    x_quant = x_quant_c2po;
end

end