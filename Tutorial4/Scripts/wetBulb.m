%% Main Function
function T_w = wetBulb(T,p,q_v)
    %% Default function parameter values    
    if nargin < 3
        T = 273.15+30;
        p = 1000;
        q_v = 0.025;
    end
    %% Constants
    c_pd = 1;           %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R_v = 0.4615;       %[kJ/kgK]
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    ep = 0.622;
    %% Calculating c_p 
    c_p = c_pd*(1+0.87*q_v);
    c_p = c_pd;
    
        
    %% Finding wet bulb temperature T_w 
    T_w = newt_raph(T,T,p,q_v);
    
    %% Newton-Raphson function
    function x = newt_raph(x_init,T,p,q_v)
        x_old = x_init;
        x_new = 1;
        err = 1;
        tol = abs(F(x_init,T,p,q_v))*1e-8;
        while err > tol
            x_new = x_old - F(x_old,T,p,q_v)/F_prim(x_old,p);
            x_old = x_new;
            
            err = abs(F(x_old,T,p,q_v));
        end        
        x = x_new;
    end
    %% F = 0 Newton-Raphson Element    
    function F = F(x,T,p,q_v)
        F = x - T + L/c_p *(q_s(x,p) - q_v);
    end
    %% dF/dx
    function F_prim = F_prim(x,p)
        F_prim = 1 + L/c_p *q_s_prim(x,p);
    end
    %% Saturated specific humidity q_s
    function q_s = q_s(T_w,p)
        e_s = e_s0*exp(L/R_v * (T_0^-1 - T_w^-1));
        q_s = ep/p * e_s;
    end
    %% dq_s/dx or equivalently dq_s/dT_w
    function q_s_prim = q_s_prim(T_w,p)
        e_s_prim = L/R_v*T_w^-2 * e_s0*exp(L/R_v * (T_0^-1 - T_w^-1));
        q_s_prim = ep/p * e_s_prim;
    end
end