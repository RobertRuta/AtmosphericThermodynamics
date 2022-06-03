function y = wetLapse()
    
    c_pd = 1;            %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R_d = 0.2869;         %[kJ/kgK]
    R_v = 0.461;         %[kJ/kgK]
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    ep = 0.622;
    G = 10;              %K/km
    g = 9.81/1000;       %km/s^2
    p = 1000;   %hPa
    q_d = 0.2;
    


    t = linspace(-40,40,100);
    T = t + 273.15;
    e_s = e_s0*exp(17.27*t./(t+237.7));
    q_s = ep*e_s/p;
    R = q_s*R_v + q_d*R_d;
    c_p = 1;
    
    %gamma_ad = c_pd/c_p * (1+q_s*B_T*R_v/R)/(1+q_s*B_T*L/(c_p*T));
    gamma_ps = (1+q_s*L./(R_d*T))./(1+q_s*L^2./(c_pd*R_v*T.^2))
    gamma = gamma_ps + 0.2*(cos(2*pi*T) + 1);
    
    z = linspace(0,30,100);
    T_env = T_0 - 10*z;
    
    close all;
    figure;
    plot(gamma_ps,t,'LineWidth',1)
    ylabel("Temperature [^oC]")
    xlabel("\gamma_{pseudo}")
    set(gca,'YDir','reverse')
    set(gcf,'position',[50 100 600 350])
    
        figure;
    plot(T_env,z,'LineWidth',1)
    ylabel("Temperature [^oC]")
    xlabel("\gamma_{pseudo}")
    set(gca,'YDir','reverse')
    set(gcf,'position',[50 100 600 350])
end