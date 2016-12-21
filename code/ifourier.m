function fouriered = ifourier (indata, x, k)

n_x = length(x);
n_k = length(k);

fouriered = zeros(1,n_x);

for k_it=1:n_k
    for x_it=1:n_x 
        fouriered(x_it) = fouriered(x_it)+ 1/(2*pi)*exp(1i*k(k_it)*x(x_it))*indata(k_it);  
    end
end

end