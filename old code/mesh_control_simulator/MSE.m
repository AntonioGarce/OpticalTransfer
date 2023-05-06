function [mse]=MSE(freq, H_master, H_nom)

MSE_i = 0;

for i=1:length(freq)
    MSE_i = MSE_i + (abs(H_master(i) - H_nom(i)))^2;
end

mse = MSE_i/length(freq);