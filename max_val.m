function[Y]=max_val(X1,X2)
Y = zeros(size(X1));
for i=1:size(X1)
   if X1(i)>=X2(i)
      Y(i)=X1(i);
   else
      Y(i)=X2(i); 
   end
end
end