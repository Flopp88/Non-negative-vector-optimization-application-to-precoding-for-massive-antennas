function[x] = saxena(s,H)
s = unit_quant(s);
M = size(H,2);
H_dag = conj(H).';
H_croix = H_dag*inv(H*H_dag);
T = H_croix;
P = T;
x = P*s;
x_quant = unit_quant(x);
q = unit_quant(H*x_quant);
if ~isequal(q,s)
   %disp("zeroforcing not sufficient");
   T_tild = abs(T).^2;
   d_square = inv(conj(T_tild).'*T_tild)*conj(T_tild).'*ones(M,1);
   D = diag(sqrt(d_square));
   P_tild = H_croix*D;
   
   x_tild = P_tild*s;
   x_quant_tild = unit_quant(x_tild);
   q_tild = unit_quant(H*x_quant_tild);
   vect = ceil((q - s)/max(q-s));
   vect_tild = ceil((q_tild - s)/max(q_tild-s));
   if norm(vect,1)>norm(vect_tild,1)
       P = P_tild;
       disp("P_tild");
   end
end
x = P*s;
end