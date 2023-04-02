function[Y]=heaviside(X)
Y = zeros(size(X));
for i=1:size(X)
    if X(i)>0
        Y(i)=1;
    else
        Y(i)=0;
    end
end
end