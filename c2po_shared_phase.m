function[x_quant,phi, d, summary_X_real,summary_X_imag, xi]=c2po_shared_phase(s,H,K,M,tau,tmax,delta)
[r_quant,x_quant,x] = zeroforcing_quant(s,H);
if r_quant == unit_quant(s)
    summary_X_real = 0;
    summary_X_imag = 0;
    xi = 0;
    phi = 0;
    d = 0;
    return
end

ampl = abs(real(s(1)));
x = conj(H).'*s;
% x = x_quant;
rho = 1/(1-delta*tau);
xi = repmat(sqrt((norm(ampl*unit_quant(x))^2)/(M)),M,1); %on a posé P=M pour obtenir une sortie sur le cercle unité
 
summary_X_real=zeros(M, tmax);
summary_X_imag=zeros(M, tmax);

summary_X_real(:,1) = real(x);
summary_X_imag(:,1) = imag(x);


for k=2:tmax
%     k
    %rho = 1/(1-delta*tau(k));
    z_real = zeros(M,1);
    z_imag = zeros(M,1);
    
    x = summary_X_real(:,k-1)+1i*summary_X_imag(:,k-1);
    R = H*x;
    a = real(s).*real(R) + imag(s).*imag(R);
    b = real(s).*imag(R) - imag(s).*real(R);
    denom = 2*sum( real(R).*imag(R).*real(s).*imag(s) );
    syst = sum(imag(s).*real(s).*( imag(R).^2 - real(R).^2)) / denom;
     
    phi = atan(syst)/2;
%     real(s*exp(1i*phi)).*real(R) + imag(s*exp(1i*phi)).*imag(R)
    
    RMSE = real(s).*( (cos(phi)^2)*a + cos(phi)*sin(phi)*b ) ...
        - imag(s).*( cos(phi)*sin(phi)*a + (sin(phi)^2)*b) - real(R);
    IMSE = imag(s).*( (cos(phi)^2)*a + cos(phi)*sin(phi)*b ) ...
        + real(s).*( cos(phi)*sin(phi)*a + (sin(phi)^2)*b) - imag(R);
    
    for l=1:M
        da_r = real(s).*real(H(:,l)) + imag(s).*imag(H(:,l));
        db_r = real(s).*imag(H(:,l)) - imag(s).*real(H(:,l));
        dphi_r =  sum( imag(s).*real(s).*( (imag(H(:,l)).*imag(R) - real(H(:,l)).*real(R))*denom ...
            + ((real(R).^2) - (imag(R).^2))*sum( imag(s).*real(s).*(real(R).*imag(H(:,l)) + ...
            imag(R).*real(H(:,l))) ))) /((denom^2)*(1+syst^2));
        A_r = -2*cos(phi)*sin(phi)*dphi_r*a + (cos(phi)^2)*da_r + ...
            cos(phi)*sin(phi)*db_r + (cos(phi)^2)*dphi_r*b - (sin(phi)^2)*dphi_r*b;
        B_r = cos(phi)*sin(phi)*da_r + (cos(phi)^2)*dphi_r*a - (sin(phi)^2)*dphi_r*a ...
            + 2*sin(phi)*cos(phi)*dphi_r*b + (sin(phi)^2)*db_r ;
        
        dRMSE_r = real(s).*A_r - imag(s).*B_r - real(H(:,l)) ;
        dIMSE_r = imag(s).*A_r + real(s).*B_r - imag(H(:,l)) ;
        dMSE_r = 2 * sum( RMSE.*dRMSE_r + IMSE.*dIMSE_r);
        z_real(l) = summary_X_real(l,k-1) - tau * dMSE_r;
        
        da_i = -real(s).*imag(H(:,l)) + imag(s).*real(H(:,l));
        db_i = da_r;
        dphi_i = sum( real(s).*imag(s).*( denom*( imag(R).*real(H(:,l)) + real(R).*imag(H(:,l)))...
            - sum( imag(s).*real(s).*(real(R).*real(H(:,l)) - imag(R).*imag(H(:,l))))* ...
            (imag(R).^2 - real(R).^2))) / ( (denom^2)*(1+syst^2)); 
        A_i = -2*cos(phi)*sin(phi)*dphi_i*a + (cos(phi)^2)*da_i + ...
            cos(phi)*sin(phi)*db_i + (cos(phi)^2)*dphi_i*b - (sin(phi)^2)*dphi_i*b;
        B_i = cos(phi)*sin(phi)*da_i + (cos(phi)^2)*dphi_i*a - (sin(phi)^2)*dphi_i*a ...
            + 2*sin(phi)*cos(phi)*dphi_i*b + (sin(phi)^2)*db_i;
        
        dRMSE_i = real(s).*A_i - imag(s).*B_i + imag(H(:,l));
        dIMSE_i = imag(s).*A_i + real(s).*B_i - real(H(:,l));
        dMSE_i = 2* sum(RMSE.*dRMSE_i + IMSE.*dIMSE_i);
        z_imag(l) = summary_X_imag(l,k-1) - tau * dMSE_i ;
    end
    
    x_real = min_val(max_val(rho*z_real,-xi),xi);
    x_imag = min_val(max_val(rho*z_imag,-xi),xi);
    summary_X_real(:,k) = x_real;
    summary_X_imag(:,k) = x_imag;
end

x_quant_c2po = xi(1)*(sign(summary_X_real(1:M,end)) + 1i*sign(summary_X_imag(1:M,end)));
r_c2po = unit_quant(H*x_quant_c2po);

if SER(s,r_c2po)<SER(s,r_quant)
    x_quant = x_quant_c2po;
end

x = summary_X_real(:,end)+1i*summary_X_imag(:,end);
R = H*x;
denom = 2*sum( real(R).*imag(R).*real(s).*imag(s) );
syst = sum(imag(s).*real(s).*( imag(R).^2 - real(R).^2)) / denom;

phi = atan(syst)/2;
% disp("shared phase:");
% disp(phi)
d =  real(s*exp(1i*phi)).*real(R) + imag(s*exp(1i*phi)).*imag(R);
% disp(d)

end