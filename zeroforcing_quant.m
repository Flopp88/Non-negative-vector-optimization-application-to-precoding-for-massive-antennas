function[r_quant,x_quant,x]=zeroforcing_quant(s,H)
s = unit_quant(s);
H_dag = conj(H).';
H_croix = H_dag*inv(H*H_dag);
P= H_croix;
x = P*s;
x_quant = unit_quant(x);
r_quant = unit_quant(H*x_quant);
end