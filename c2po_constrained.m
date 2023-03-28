function[x_quant,summary_X_real,summary_X_imag, xi]=c2po_constrained(s,H,K,M,tau,tmax,delta)
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
    
    phi = max(-pi/4,min(pi/4,atan(b./a)));
    r = max(0,cos(phi).*a+sin(phi).*b);
    
%     RMSE = real(s).*( (cos(phi)^2)*a + cos(phi)*sin(phi)*b ) ...
%         - imag(s).*( cos(phi)*sin(phi)*a + (sin(phi)^2)*b) - real(R);
%     IMSE = imag(s).*( (cos(phi)^2)*a + cos(phi)*sin(phi)*b ) ...
%         + real(s).*( cos(phi)*sin(phi)*a + (sin(phi)^2)*b) - imag(R);
    RMSE = r.*(cos(phi).*real(s) - sin(phi).*imag(s)) - real(R);
    IMSE = r.*(sin(phi).*real(s) + cos(phi).*imag(s)) - imag(R);
    
    for l=1:M
        da_r = real(s).*real(H(:,l)) + imag(s).*imag(H(:,l));
        db_r = real(s).*imag(H(:,l)) - imag(s).*real(H(:,l));
        darctan_r = ( db_r.*a - b.*da_r )./((1+(b./a).^2).*(a.^2));
        
        dphi_r =  (darctan_r/4) - sign(atan(b./a) - (pi/4)).*darctan_r/4 ...
            - sign( (atan(b./a) - abs(atan(b./a)-(pi/4)))/2 +(pi/2))...
            .*(darctan_r/2 -sign(atan(b./a)-(pi/4)).*darctan_r/2);
        
        dr_r = heaviside(cos(phi).*a+sin(phi).*b).*(cos(phi).*da_r - ...
            a.*sin(phi).*dphi_r + sin(phi).*db_r + b.*cos(phi).*dphi_r);
        
        dRMSE_r = real(s).*(cos(phi).*dr_r -r.*sin(phi).*dphi_r) - ...
            imag(s).*(r.*cos(phi).*dphi_r + sin(phi).*dr_r) - real(H(:,l));
        dIMSE_r = real(s).*(r.*cos(phi).*dphi_r + sin(phi).*dr_r) + ...
            imag(s).*(cos(phi).*dr_r - r.*sin(phi).*dphi_r) - imag(H(:,l));
        
        dMSE_r = 2 * sum( RMSE.*dRMSE_r + IMSE.*dIMSE_r);
        z_real(l) = summary_X_real(l,k-1) - tau * dMSE_r;
        
        da_i = -real(s).*imag(H(:,l)) + imag(s).*real(H(:,l));
        db_i = da_r;
        darctan_i = ( db_i.*a - b.*da_i )./((1+(b./a).^2).*(a.^2));
        
        dphi_i =  (darctan_i/4) - sign(atan(b./a) - (pi/4)).*darctan_i/4 ...
            - sign( (atan(b./a) - abs(atan(b./a)-(pi/4)))/2 +(pi/2))...
            .*(darctan_i/2 -sign(atan(b./a)-(pi/4)).*darctan_i/2);
        
        dr_i = heaviside(cos(phi).*a+sin(phi).*b).*(cos(phi).*da_i - ...
            a.*sin(phi).*dphi_i + sin(phi).*db_i + b.*cos(phi).*dphi_i);
        
        dRMSE_i = real(s).*(cos(phi).*dr_i -r.*sin(phi).*dphi_i) - ...
            imag(s).*(r.*cos(phi).*dphi_i + sin(phi).*dr_i) + imag(H(:,l));
        dIMSE_i = real(s).*(r.*cos(phi).*dphi_i + sin(phi).*dr_i) + ...
            imag(s).*(cos(phi).*dr_i - r.*sin(phi).*dphi_i) - real(H(:,l));
        
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


end