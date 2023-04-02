function[y]=sigmoid(x,lambda)
y = zeros(size(x));
for i = 1:length(x)
    y(i) = 1/(1+exp(-lambda*x(i))); 
end
end