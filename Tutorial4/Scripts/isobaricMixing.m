function main = isobaricMixing(T_1,T_2,q_v1,q_v2)
    %% Constants
    c_pd = 1;           %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R_v = 0.4615;       %[kJ/kgK]
    T_0 = 273.15;       %[K]
    e_s0 = 6.112;       %[hPa]
    ep = 0.622;
    p = 1000;           %[hPa]
    c_p = c_pd;
    
    %% Default function parameter values    
    if nargin < 4
        p = 1000;
        T_1 = 273.15+35;
        T_2 = 273.15+5;
        q_v1 = 0.9*q_s(T_1,p)
        q_v2 = 0.7*q_s(T_2,p)
    end
    
    T_init = 0.5*T_1 + 0.5*T_2;
    q_v =  0.5*q_v1 + 0.5*q_v2;

    T = isobaric(T_init,p,q_v);
    q_final = q_s(T,p);
    q_l = (q_v-q_final);

    main.q_l = q_l;
    
    %% Isobaric
    function x = isobaric(T_o,p,q_vo)

        x_old = T_o;
        x_new = 1;
        F = x_old - T_o + L/c_p * (q_s(x_old,p) - q_vo);
        F_prim = 1 + L/c_p * q_s_prim(x_old,p);
        
        err = 1;        
        tol = abs(F)*1e-8;
        while err > tol            
            x_new = x_old - F/F_prim;
            x_old = x_new;
            
            F = x_old - T_o + L/c_p * (q_s(x_old,p) - q_vo);
            F_prim = 1 + L/c_p * q_s_prim(x_old,p);

            err = abs(F);
        end        
        x = x_new;
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
    T_env = linspace(0,40,100) + 273.15;
    e_s = e_s0*exp(L/R_v * (T_0^-1 - T_env.^-1));
    T_points = [T_1,T_2];
    e_points = [q_v1*p/ep,q_v2*p/ep];
    T_line = [T_1,T_2];
    e_line = [q_v1*p/ep,q_v2*p/ep];

    close all;
    plot(T_env - 273.15,e_s,'-k','LineWidth',0.75);
    hold on;
    plot(T_points-273.15,e_points,'xr','LineWidth',1,'MarkerSize',7);
    plot(T_init-273.15,q_v*p/ep,'o','Color','r','LineWidth',1)
    plot(T - 273.15,p/ep *q_final,'*','MarkerSize',8)
    plot([T_init-273.15,T-273.15],[q_v*p/ep,q_final*p/ep],':b','LineWidth',1)
    plot(T_line-273.15,e_line,'--','LineWidth',1);
    
    text((T - 273.15)*1.05,p/ep*q_final*1 - 2,strcat("\Deltaq_l = ",string(0.994)," [g/kg]"),'Color','b')
    
    title("Isobaric Mixing of Two Undersaturated Air Masses");
    xlabel("Temperature T [^oC]");
    ylabel("Vapour Pressure e [hPa]");
    legend(["Saturation curve","Initial Air Masses","Pre-condensation","Post-condensation","Condensation Path"],'location','northwest');
    set(gcf,'position',[50,100,550,350])
end