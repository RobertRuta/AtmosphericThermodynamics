function y = FindLCL()
    
    %% Constants
    c_p = 1;            %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R = 0.2869;         %[kJ/kgK]
    R_v = 0.461;         %[kJ/kgK]
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    ep = 0.622;
    p_0 = 1000;         %hPa
    G_env = 8;          %K/km
    G_d = 10;           %K/km
    g = 9.81;           %m/s^2
    %% Ground conditions
    T_env_0 = 20;
    T_parc_0 = 20;
    T_o = linspace(-10,30,100) + 273.15;
    q_v = 7.48;   %g/kg
    
    z = linspace(0,10,1000);
    T_o = T_env(z,T_env_0);
    f = (q_v) ./ q_s(z,T_env_0) * 100;
    one = min(abs(100 - f));
    idx_LCL = (abs(100 - f) == one);
    z_LCL = z(idx_LCL);
    p_LCL = p(z_LCL);
    T_LCL = T_env(z_LCL,T_env_0);
    
    T_parcel = T_parc(z, T_parc_0);
    one = (abs(T_parcel - T_o)) < 0.03;    
    z_peaks = z(one);
    z_LFC = z_peaks(3);
    p_LFC = p(z_LFC);
    T_LFC = T_env(z_LFC,T_env_0);
    idx_LFC = find(z == z_LFC);
    
    z_LNB = z_peaks(4);
    p_LNB = p(z_LNB);
    T_LNB = T_env(z_LNB,T_env_0);
    idx_LNB = find(z == z_LNB);
    
    Z = z(idx_LFC:idx_LNB);
    lnp = log(p(Z));
    Q1 = sum(diff(lnp).*T_parcel(idx_LFC+1:idx_LNB));
    Q2 = sum(diff(lnp).*T_o(idx_LFC+1:idx_LNB));
    CAPE = -R*(Q1-Q2);
    
    
    y.LCL = [z_LCL,p_LCL,T_LCL];
    y.LFC = [z_LFC,p_LFC,T_LFC];
    y.LNB = [z_LNB,p_LNB,T_LNB];
    y.CAPE = CAPE;
    
    
    
    %% Nested functions
    function F = F(x,T_o,f_0)
        F = log(f_0) + 1/R * (c_p * log(x/T_o) + ep*L * (1/x - 1/T_o));
    end

    function F_prim = F_prim(x)
        F_prim = 1/R * (c_p/x - ep*L/(x^2));
    end

    function X = newt_raph(T,f)
        x_old = T;
        err = 1;
        tol = abs(F(x_old,T,f))*1e-8;
        
        while err > tol
            x_new = x_old - F(x_old,T,f)/F_prim(x_old);
            x_old = x_new;
            
            err = abs(F(x_old,T,f));
        end
        X = x_old;
    end

    function q_s = q_s(z,T_env_0)
        T = T_env_0 + 273.15 - G_env*z;
        t = T - 273.15;
        e_s = e_s0*exp(17.27*t./(t+237.7));
        q_s = ep*e_s./p(z) * 1000;
    end

    function p = p(z)
        p = p_0 * (1 - G_env/T_0 * z).^(g/(R*G_env));
    end

    function T_env = T_env(z,T_env_0)
        T_env = T_env_0 - G_env*z;
    end

    function T_parc = T_parc(z,T_parc_0)
        logic_d = z < z_LCL;
        logic_s = z >= z_LCL;
        z_1 = z .* logic_d;
        z_2 = (z - z_LCL) .* logic_s;
        
        index = find(abs(diff(logic_d)) == 1);
        
        T_1 = (T_parc_0 - G_d*z_1).* logic_d;
        T_parc_LCL = T_1(index);
        To = T_parc_LCL;
        for i = 1:numel(z_2)
            q = ep*e_s(To+273.15)/200;
            G = G_d*(1+q*L./(R*(To+273.15)))./(1+q*L^2./(c_p*R_v*(To+273.15).^2));
            
            gamma(i) = (1+q*L./(R*(To+273.15)))./(1+q*L^2./(c_p*R_v*(To+273.15).^2));
            
            T_2(i) = ( T_parc_LCL - G*z_2(i) )*logic_s(i);
            To = T_2(i);
        end
        
        T_parc = T_1 + T_2;        
    end

    function G_s = G_s(z,T)
        T = T + 273.15;    
        %lambda =(1 + L*q_s(z,T)/1000/(R*T)) / (1 + L^2*q_s(z,T)/1000*ep/(R*T^2*c_p));
        gamma_ps = (1+q_s*L./(R*T))./(1+q_s*L^2./(c_p*R_v*T.^2))
        G_s = lambda*G_d;
    end

    function e_s = e_s(T)
        t = T - 273.15;
        e_s = e_s0*exp(17.27*t./(t+237.7));
    end
end
