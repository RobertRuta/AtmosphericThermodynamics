function y = latentHeats()

    %Real data
    T_v = [-50:10:60];
    T_s = [-60:10:0]
    L_v = [629.3,621.7,615.0,608.9,603,597.3,591.7,586,580.4,574.7,569,563.2] * 4.1868; %kJ/kg
    L_s = [677.5,677.9,678,678,677.9,677.5,677] * 4.1868; %kJ/kg

    %Constants
    c_pl = 4.218;   %kJ/kgK
    c_pi = 2.106;   %kJ/kgK
    c_pv = 1.870;   %kJ/kgK
    L_lv0 = 2501;   %kJ/kg
    L_iv0 = 2834;   %kJ/kg
    T_0 = 273.15;   %K

    %Data Generation
    T_iv = linspace(-60,0,100) + 273.15;    %K
    T_lv = linspace(-60,60,100) + 273.15;    %K
    L_lv = ( L_lv0 - (c_pl - c_pv)*(T_lv - T_0) ) / 1000;
    L_iv = ( L_iv0 - (c_pi - c_pv)*(T_iv - T_0) ) / 1000;

    %Const Error
    ERR_V = (L_lv0/1000 - L_lv) ./ L_lv;
    ERR_S = (L_iv0/1000 - L_iv) ./ L_iv;

    %True Error
    %L_V = interp1(T_v,L_v,T_lv-273.15);
    %L_S = interp1(T_s,L_s,T_iv-273.15);
    L_V = ( L_lv0 - (c_pl - c_pv)*(T_v+273.15 - T_0) ) / 1000;
    L_S = ( L_iv0 - (c_pi - c_pv)*(T_s+273.15 - T_0) ) / 1000;
    err_V = (L_v/1000 - L_V) ./ L_V;
    err_S = (L_s/1000 - L_S) ./ L_S;

    %Plotting
    close all;
    plot(T_lv-273.15,L_lv,'--b', 'linewidth',1.5);
    hold on;
    plot(T_iv-273.15,L_iv,':r', 'linewidth',1.5);

    xlim([-40 40]);
    %ylim([2.3 3]);
    title('Latent Heat as a Function of Temperature');
    xlabel('Temperature [^oC]');
    ylabel('Latent Heat [10^6 J/kg]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"Vap","Sub"});
    set(gcf,'position',[50,550,500,320]);

    figure;
    plot(T_lv-273.15,ERR_V*100,'--b', 'linewidth',1.5);
    hold on;
    plot(T_iv-273.15,ERR_S*100,':r', 'linewidth',1.5);

    xlim([-40 40]);
    %ylim([2.3 3]);
    title('Graph of Errors - Constant Latent Heat');
    xlabel('Temperature [^oC]');
    ylabel('Error [%]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"Vap","Sub"},'location','northwest');
    set(gcf,'position',[50,550,500,250]);

    figure;
    plot(T_v,err_V*100,'ob', 'linewidth',1.5);
    hold on;
    plot(T_s,err_S*100,'xr', 'linewidth',1.5);

    %xlim([-60 60]);
    %ylim([2.3 3]);
    title('Graph of Errors - SMC Latent Heat');
    xlabel('Temperature [^oC]');
    ylabel('Error [%]');
    set(gca,'XMinorTick','on','YMinorTick','on'); 
    legend({"Vap","Sub"});
    set(gcf,'position',[50,550,500,250]);

    y.lv = L_lv;
    y.iv = L_iv;
end