%% Main Function
function T_w = wetBulb()
    %% Constants
    c_pd = 1;           %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R_v = 0.4615;       %[kJ/kgK]
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    ep = 0.622;
    p = 1000            %hPa
    
    %% Generating domains
    t = linspace(0,40,100) + 273.15;    %[K]
    Q_v = linspace(0,30,100) /1000;     %[kg/kg]
    
    %% Calculating c_p 
    c_p = c_pd*(1+0.87*Q_v);
    
    %% Finding wet bulb temperature T_w 
    T_w = newt_raph();
    
    %% Newton-Raphson function
    function f = newt_raph()
        for i = 1:numel(t)
            for j = 1:numel(Q_v)            
                T = t(i);
                q_v = Q_v(j);
                
                x_init = T;
                x_old = x_init;
                err = 1;
                tol = abs(F(x_init,T,p,q_v))*1e-8;
                while err > tol
                    x_new = x_old - F(x_old,T,p,q_v)/F_prim(x_old,p);
                    x_old = x_new;

                    err = abs(F(x_old,T,p,q_v));
                end
                f(j,i) = x_old;
            end
        end
    end
    %% F = 0 Newton-Raphson Element    
    function F = F(x,T,p,q_v)
        F = x - T + L/c_pd *(q_s(x,p) - q_v);
    end
    %% dF/dx
    function F_prim = F_prim(x,p)
        F_prim = 1 + L/c_pd *q_s_prim(x,p);
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
    %% Plotting
    close all;
    [C,h] = contourf(t-273.15,Q_v*1000,T_w-273.15);
    clabel(C,h);
    
    
    %% patch
    X1 = t - 273.15;
    for i = 1:numel(t)
        Y1(i) = q_s(X1(i)+273.15,1000);
    end
    X2 = fliplr(X1);
    Y2 = zeros(1,numel(t)) + Y1(end);
    X3 = zeros(1,numel(t)) + t(1) - 273.15;
    Y3 = fliplr(Y1);
    
    X = horzcat(X1,X2,X3);
    Y = horzcat(Y1,Y2,Y3)*1000;
    
    patch(X,Y,'w');
    ylim([0 30]);
    title("Wet-bulb temperature for Different Initial Properties [^oC]")
    xlabel("Temperature T [^oC]")
    ylabel("Saturation Specific Humidity q_s [g/kg]")
    set(gcf,'position',[50,100,550,350])
    
end