%% Newton-Raphson method (finding roots of function)
function x = newt_raph(x_init,T,p,q_v)
    x_old = x_init;    
    err = 1;
    tol = abs(F(x_init,T,p,q_v))*1e-8;
    
    while err > tol
        x_new = x_old - F(x_old,T,p,q_v)/F_prim(x_old,p);
        x_old = x_new;
        
        err = abs(F(x_old,T,p,q_v));
    end
    
    x = x_new;
end