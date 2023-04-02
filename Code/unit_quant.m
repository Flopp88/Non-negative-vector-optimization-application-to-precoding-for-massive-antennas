function[x_quant] = unit_quant(x)
x_quant = (sign(real(x)) + 1i*sign(imag(x)))/sqrt(2);
end