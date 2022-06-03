function main = Equivalent()
    %% Constants
    c_pd = 1;           %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R_v = 0.4615;       %[kJ/kgK]
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    ep = 0.622;
    g = 9.81;      %m/s^2
    R_d = 8.3145/0.02897;   %m^2/(s^2 K)
    p_0 = 1000;         %hPa
    k = 0.287;

    c_p = c_pd;

    z = linspace(0,10,100);
    p_env = p_0*(1 - 6/T_0 * z).^(g/(R_d*6e-3));
    T_env = T_0 - 6*z;
    env(:, 1) = T_env;
    env(:, 2) = p_env;
    %% Data Generation
    T = linspace(-40,40,100) + 273.15;
    P = linspace(100,1000,100);

    main.TH_e_pot = TH_pot_e(T,P);
    main.TH_pot_e = TH_e_pot(T,P);
    
    %% Plotting
    close all;
    v = [250,300,350,400,500,600,800,1200,2000]
    [C,h] = contour(P,T-273.15,TH_e_pot(T,P),v);
    clabel(C,h,v);
    xlabel("Pressure p [hPa]");
    ylabel("Temperature T [^oC]")
    title("Equivalent Potential Temperature")
    set(gcf,'position',[50,50,650,250])
    
    figure;
    [C,h] = contour(P,T-273.15,TH_pot_e(T,P),v);
    clabel(C,h,v);
    xlabel("Pressure p [hPa]");
    ylabel("Temperature T [^oC]")
    title("Potential Equivalent Temperature")
    set(gcf,'position',[50,50,650,250])
    
    figure;
    difference = [0.01,0.1,1,10,100,1000,10000,100000]
    [C,h] = contour(P,T-273.15,-TH_pot_e(T,P) + TH_e_pot(T,P),difference);
    clabel(C,h,difference);
    xlabel("Pressure p [hPa]");
    ylabel("Temperature T [^oC]")
    title("Difference")
    set(gcf,'position',[50,50,650,250])
    
    %% Nested functions
    function TH_pot_e = TH_pot_e(T,P)
        Qs = q_s(T,P);
        T_e = T.' + L/c_p * Qs;
        TH_pot_e = bsxfun(@times,T_e,(p_0./P).^k);
    end
    function TH_e_pot = TH_e_pot(T,P)
        Qs = q_s(T,P);
        T = T.';           
        TH = bsxfun(@times,T,(p_0./P).^k);
        TH_e_pot =  TH .* exp(bsxfun(@times,Qs,L./(c_p*T)));
    end

    %% Saturated specific humidity q_s
    function q_s = q_s(T,P)
        e_s = e_s0*exp(L/R_v * (T_0^-1 - T.^-1));
        e_s = e_s.';
        q_s = bsxfun(@times,ep./P,e_s);
    end
end