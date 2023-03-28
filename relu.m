function[Y]=relu(X)
Y = zeros(size(X));
for i=1:length(X)
    if X(i)>0
        Y(i)=X(i);
    end
end
end