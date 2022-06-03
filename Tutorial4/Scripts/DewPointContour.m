function fun = DewPoint()

%%Constants
    L_lv0 = 2501;       %kJ/kg
    L_iv0 = 2834;       %kJ/kg
    c_pl = 4.218;       %kJ/kgK
    c_pi = 2.106;       %kJ/kgK
    c_pv = 1.870;       %kJ/kgK
    e_s0 = 6.112;       %hPa
    R_v = 461.5/1000;   %kJ/kgK
    ep = 0.622;         %Unitless

    T = linspace(-40,40,100);              %[K]
    T_d = linspace(0,40,100).' + 273.15;    %[K]
    f = exp(-L_lv0/R_v * bsxfun(@times,(T - T_d),1./bsxfun(@times,T,T_d)));
    deficit = T - T_d;
    
%%Plotting
    close all;
    figure;
    contour(T,T_d,f);
    %xlim([0 40])
    %ylim([0 40])
    title('Dew Point Deficit with T');
    xlabel('Temperature [^oC]');
    ylabel('Dew Point Deficit [^oC]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"Lin.L(T)","Const.L(T)","Magnus-Tentens"},'location','northwest');
    set(gcf,'position',[50,550,550,300]);

end