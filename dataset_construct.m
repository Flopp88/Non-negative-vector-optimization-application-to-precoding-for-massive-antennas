clear all
close all

K = 5;
M = 20;
size_data = 1000;

Dataset = zeros(size_data*K,M+1);
count_set_size=0;
count_iteration = 0;

while count_set_size<size_data
    H = wgn(K,M,0)+1i*wgn(K,M,0);
    s = (randsrc(K,1)+1i*randsrc(K,1))/sqrt(2);
    [r_quant,x_quant] = zeroforcing_quant(s,H);
    if ~isequal(s,r_quant)
        append = 1;
        for i=0:count_set_size
            if (isequal(Dataset((i)*K + (1:K),1),s)) && (isequal(Dataset((i)*K + (1:K),2:end),H))
                append = 0;
                break
            end
        end
        if append
            Dataset((count_set_size)*K + (1:K),1) = s;
            Dataset((count_set_size)*K + (1:K),2:end) = H;
            count_set_size = count_set_size + 1
        end
    end
    count_iteration = count_iteration + 1;
end

save('./Datasets/Dataset1.mat', 'Dataset');

%%

H_test = Dataset((200)*K + (1:K),2:end);
s_test = Dataset((200)*K + (1:K),1);
[r_quant_test,x_quant_test] = zeroforcing(s_test,H_test);

disp(r_quant_test/norm(r_quant_test) - s_test/norm(s_test))

tau = 0.001;
tmax = 1000;
[x_quant_test,X_out_test, xi_test] = original_c2po(s_test, H_test ,K ,M ,tau ,tmax);
r_test = H_test*x_quant_test;
r_quant_test = unit_quant(r_test);
disp(r_quant_test/norm(r_quant_test) - s_test/norm(s_test))




