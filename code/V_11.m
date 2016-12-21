function v = V_11(x,a,b,c,d) 
if (x > 0) 
    v =a*(2-exp(-x/b)); 
else
    v =a*exp(x/b);
end
end




