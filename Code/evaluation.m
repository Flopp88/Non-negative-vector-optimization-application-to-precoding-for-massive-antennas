clear all
% close all

path = "./Datasets";
data_struct = dataloader(path);

index =5;

data1 = data_struct(index).data;
K = data_struct(index).K;
M = data_struct(index).M;

% K=5;
% M=25;

nb_sets = 1000;
%SNR_db = cat(2,-30:10:-10,-5:0.25:5,10:5:30);
SNR_db = -30:2:50;

SNR = 10.^(SNR_db/10);
rho = 1;
rho_0 = rho*M;

SNR_results_ZF = zeros(length(SNR),nb_sets);
SNR_results_Buss = zeros(length(SNR),nb_sets);
SNR_results_C2PO = zeros(length(SNR),nb_sets);
SNR_results_C2PO_orig = zeros(length(SNR),nb_sets);
SNR_results_C2PO_opti = zeros(length(SNR),nb_sets);
SNR_results_C2PO_orig_opti = zeros(length(SNR),nb_sets);
SNR_results_C2PO_sigm = zeros(length(SNR),nb_sets);
SNR_results_C2PO_sigm_opti = zeros(length(SNR),nb_sets);

SNR_results_C2PO_shared = zeros(length(SNR),nb_sets);
SNR_results_C2PO_constrained = zeros(length(SNR),nb_sets);
SNR_results_C2PO_simple = zeros(length(SNR),nb_sets);

d_vect = zeros(K, 1);
alpha_vect = zeros(1, 1);

for i=1:nb_sets
    if mod(i,100)==0
        disp(i)
    end
    set = data1( (1:K)+K*(i-1), :); 
    s = set(1:K,1);
    H = set(1:K,2:end);
%     H = (wgn(K,M,0)+1i*wgn(K,M,0));
%     s = (randsrc(K,1)+1i*randsrc(K,1))/sqrt(2);

    [~,x_zero] = zeroforcing_quant(s,H);
    x_quant_zero = unit_quant(x_zero);
    

%     x_buss = saxena(s,H);
%     x_quant_buss = unit_quant(x_buss);
    

    %[x_quant_c2po,~,~,~] = c2po(s,H,K,M,0.000000001,100,1);
    
    [x_quant_c2po_opti,~,~,~] = c2po_opti(s,H,K,M,0.000000001,100,1);
    
 
    %[x_quant_orig_c2po,~,~] = original_c2po(s,H,K,M,0.001,100);
    
    [x_quant_orig_c2po_opti,phi_orig,alpha,~,~] = original_c2po_opti(s,H,K,M,0.001,100);
    
    
    %[x_quant_c2po_sigm,~,~] = c2po_sigm(s,H,K,M,0.000000001,100,100,1000);
    
    [x_quant_c2po_sigm_opti,~,~] = c2po_sigm_opti(s,H,K,M,0.000000001,100,100,100);
    
    
    tau = fliplr(1:1000);
    tau = tau*0.01/1000;
    [x_quant_c2po_shared,phi_shared,d,~,~] = c2po_shared_phase(s,H,K,M,0.01,100,1);
    %d_vect(:,i) =  d;
    if (-pi/4)>=phi_shared || phi_shared>=(pi/4)
        disp(phi_shared)
    end
    
    [x_quant_c2po_constrained,~,~] = c2po_constrained(s,H,K,M,0.1,100,1);
    
    [x_quant_c2po_simple] = c2po_constrained(s,H,K,M,0.01,100,1);
    
    for k=1:length(SNR) %1: length(SNR)
        sigma_sq = 2*rho_0/SNR(k);
        
        noise = sqrt(sigma_sq)*(wgn(K,1,0)+1i*wgn(K,1,0));
        %noise = zeros(K,1);
        
        r_zero = H*x_quant_zero + noise;
        r_zero = unit_quant(r_zero);
        ser_zero = SER(s,r_zero);
        SNR_results_ZF(k,i) = ser_zero;
        
%         r_buss1 = H*x_quant_buss + noise;
%         r_buss = unit_quant(r_buss1);
%         ser_buss = SER(s,r_buss);
%         SNR_results_Buss(k,i) = ser_buss;
        
%         r_c2po = H*x_quant_c2po + noise;
%         r_c2po = unit_quant(r_c2po);
%         ser_c2po = SER(s,r_c2po);
%         SNR_results_C2PO(k,i) = ser_c2po;
        
        r_c2po_opti = H*x_quant_c2po_opti + noise;
        r_c2po_opti = unit_quant(r_c2po_opti);
        ser_c2po_opti = SER(s,r_c2po_opti);
        SNR_results_C2PO_opti(k,i) = ser_c2po_opti;
        
%         r_orig_c2po = H*x_quant_orig_c2po + noise;
%         r_orig_c2po = unit_quant(r_orig_c2po);
%         ser_orig_c2po = SER(s,r_orig_c2po);
%         SNR_results_C2PO_orig(k,i) = ser_orig_c2po;
        
        
        r_orig_c2po_opti = H*x_quant_orig_c2po_opti + noise;
        r_orig_c2po_opti_quant = unit_quant(r_orig_c2po_opti);
        ser_orig_c2po_opti = SER(s,r_orig_c2po_opti_quant);
        SNR_results_C2PO_orig_opti(k,i) = ser_orig_c2po_opti;
        
%         r_c2po_sigm = H*x_quant_c2po_sigm + noise;
%         r_c2po_sigm = unit_quant(r_c2po_sigm);
%         ser_c2po_sigm = SER(s,r_c2po_sigm);
%         SNR_results_C2PO_sigm(k,i) = ser_c2po_sigm;
        
        r_c2po_sigm_opti = H*x_quant_c2po_sigm_opti + noise;
        r_c2po_sigm_opti = unit_quant(r_c2po_sigm_opti);
        ser_c2po_sigm_opti = SER(s,r_c2po_sigm_opti);
        SNR_results_C2PO_sigm_opti(k,i) = ser_c2po_sigm_opti;
                
        r_c2po_shared = H*x_quant_c2po_shared + noise;      
        r_c2po_shared_quant = unit_quant(r_c2po_shared);        
        ser_c2po_shared = SER(s,r_c2po_shared_quant);        
        SNR_results_C2PO_shared(k,i) = ser_c2po_shared;
        
        if ser_c2po_shared~=0 %|| ser_orig_c2po_opti~=0
%             disp("orig")
%             disp("r")
%             disp(r_orig_c2po_opti)
%             disp("r_quant")
%             disp(r_orig_c2po_opti_quant)
%             disp("SER")
%             disp(ser_orig_c2po_opti)
%             disp("Phi orig")
%             disp(phi_orig)
%             disp("alpha")
%             disp(alpha)
%             
%             disp("shared")
%             disp("r")
%             disp(r_c2po_shared)
%             disp("r_quant")
%             disp(r_c2po_shared_quant)
%             disp("SER")
%             disp(ser_c2po_shared)
%             disp("phi shared")
%             disp(phi_shared)
%             disp("d")
%             disp(d)
            d_vect = [d_vect d];
            alpha_vect = [alpha_vect alpha];
        end
        
        r_c2po_constrained = H*x_quant_c2po_constrained + noise;      
        r_c2po_shared_constrained = unit_quant(r_c2po_constrained);        
        ser_c2po_constrained = SER(s,r_c2po_shared_constrained);        
        SNR_results_C2PO_constrained(k,i) = ser_c2po_constrained;
        
        r_c2po_simple = H*x_quant_c2po_simple + noise;      
        r_c2po_shared_simple = unit_quant(r_c2po_simple);        
        ser_c2po_simple = SER(s,r_c2po_shared_simple);        
        SNR_results_C2PO_simple(k,i) = ser_c2po_simple;
    end
end



mean_ser_ZF = mean(SNR_results_ZF,2);
% mean_ser_Buss = mean(SNR_results_Buss,2);
% mean_ser_C2PO = mean(SNR_results_C2PO,2);
mean_ser_C2PO_opti = mean(SNR_results_C2PO_opti,2);
% mean_ser_C2PO_orig = mean(SNR_results_C2PO_orig,2);
mean_ser_C2PO_orig_opti = mean(SNR_results_C2PO_orig_opti,2);
% mean_ser_C2PO_sigm = mean(SNR_results_C2PO_sigm,2);
mean_ser_C2PO_sigm_opti = mean(SNR_results_C2PO_sigm_opti,2);

mean_ser_C2PO_shared = mean(SNR_results_C2PO_shared,2);
mean_ser_C2PO_constrained = mean(SNR_results_C2PO_constrained,2);
mean_ser_C2PO_simple = mean(SNR_results_C2PO_simple,2);

figure()
%loglog(SNR_db,mean_ser_ZF,SNR_db,mean_ser_Buss,SNR_db,mean_ser_C2PO,SNR_db,mean_ser_C2PO_orig)
% ,SNR_db,mean_ser_Buss ,"Saxena"
semilogy(SNR_db,mean_ser_ZF,SNR_db,mean_ser_C2PO_opti...
    ,SNR_db,mean_ser_C2PO_orig_opti,SNR_db,mean_ser_C2PO_sigm_opti...
    ,SNR_db, mean_ser_C2PO_shared,SNR_db,mean_ser_C2PO_constrained,SNR_db,mean_ser_C2PO_simple)
% ,SNR_db,mean_ser_C2PO ,SNR_db,mean_ser_C2PO_orig ,SNR_db,mean_ser_C2PO_sigm
legend("Quantized ZF ","Positive Vector C2PO","Original C2PO","Sigmoid C2PO"...
    , "Shared phase C2PO", "Constrained C2PO","simple C2PO")
% ,"C2PO","C2PO original","C2PO sigmoid"
grid on
title(['Simulation results for K=',num2str(K),' and M=',num2str(M)])
ylim([10^-2 1])
ylabel("Average SER")
xlabel("SNR")

% d_vect= d_vect(:,2:end);
% alpha_vect = alpha_vect(2:end);
% d_vect = sort(d_vect,1,'descend');%/mean(alpha_vect);
% var_d = var(d_vect,1,2);
% %d_vect = mean(d_vect,2);
% figure()
% hold on
% %errorbar(d_vect,var_d)
% 
% stem(d_vect(:,4)/alpha_vect(4))
% yline(1,'r','LineWidth',2)
% hold off
% title("Mean values of the Sorted amplification vector compared to alpha for K=5 and M=25")
% ylabel("Amplitude of the amplification")
% xlabel("Elements of the amplification sorted by amplitude")
% legend("Elements of d","Amplitude of alpha")





