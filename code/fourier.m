function fouriered = fourier (indata, x, k)

n_x = length(x);
n_k = length(k);

fouriered = zeros(1,n_k);

for k_it=1:n_k
    for x_it=1:n_x 
        fouriered(k_it) = fouriered(k_it)+ exp(1i*k(k_it)*x(x_it))*indata(x_it);  
    end
end

end