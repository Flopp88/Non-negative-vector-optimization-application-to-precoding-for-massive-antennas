function[ser] = SER(s,r)
err = abs(r-s);
idx = err==0;
nb_true = sum(idx(:));
ser = (length(s) - nb_true)/length(s);
end