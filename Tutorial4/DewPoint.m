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

    f = linspace(0.1,0.9,9).';              %[\%/100]
    T = linspace(0,40,100) + 273.15;    %[K]
    T_d = bsxfun(@times,L_lv0*T, 1./(L_lv0 - R_v*T.*log(f)));  %[K]
    deficit = T - T_d;
    
    T2 = linspace(5,40,8).' + 273.15;                %[K]
    T_d2 = linspace(-20,40,100) + 273.15;                 %[K]
    f2 = exp(-L_lv0/R_v * bsxfun(@times,T2-T_d2,1./bsxfun(@times,T2,T_d2)));
    
%%Plotting
    close all;
    
    figure;
    for i = 1:numel(f)
        h1 = plot(T-273.15,deficit(i, :),'LineWidth',.75);
        hold on;
        label(h1,strcat({string(f(i)*100) + "%"}),'location','right','FontSize',7,'verticalalignment','top')
    end
    xlim([-1 41])
    ylim([-1 41])
    title('Dew Point Deficit');
    xlabel('Temperature T [^oC]');
    ylabel('T - T_d [^oC]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    set(gcf,'position',[50,550,500,300]);
    
    figure;
    for i = 1:numel(T2)
        h = plot(T2(i) - T_d2,f2(i,:)*100,'LineWidth',.5);
        hold on;
        label(h,strcat({string(T2(i)-273.15) + "^oC"}),'location','right','FontSize',7,'verticalalignment','top')
    end
    xlim([0 60])
    ylim([-10 100])
    title('Relative Humidity');
    xlabel('T - T_d [^oC]');
    ylabel('f [%]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    set(gcf,'position',[50,550,500,300]);
    
    
    

end