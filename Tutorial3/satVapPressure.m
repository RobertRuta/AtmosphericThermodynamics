function f = satVapPressure()
%%Constants
    L_lv0 = 2501;       %kJ/kg
    L_iv0 = 2834;       %kJ/kg
    c_pl = 4.218;       %kJ/kgK
    c_pi = 2.106;       %kJ/kgK
    c_pv = 1.870;       %kJ/kgK
    T_0 = 273.15;       %K
        e_s0 = 6.112;   %hPa
    R_v = 461.5/1000;   %kJ/kgK
    
%%Real data
    e_i = [12.85,22.36,38.02,63.3,103.28,165.32,259.92,401.78,611.15];
    t_i = [-40,-35,-30,-25,-20,-15,-10,-5,0]

%%Functions
    T = linspace(-40,40,500) + 273.15; %K
    t = T - 273.15; %oC
    e_s = e_s0*exp( (L_lv0/R_v) * (1/T_0 - 1./T) );
    e_sL = e_s0*exp( (L_lv0 + T_0*(c_pl - c_pv))/R_v * (1/T_0 - 1./T) - (c_pl-c_pv)/R_v *log(T/T_0));
    e_st = e_s0*exp(17.27*t./(t+237.7));

    Ti = linspace(-40,0,250) + 273.15; %K
    ti = Ti - 273.15; %oC
    e_si = e_s0*exp( (L_iv0/R_v) * (1/T_0 - 1./Ti) );
    e_sLi = e_s0*exp( (L_iv0 + T_0*(c_pi - c_pv))/R_v * (1/T_0 - 1./Ti) - (c_pi-c_pv)/R_v *log(Ti/T_0));
    
%%Goodness of Fit Testing
    chi2_1 = sum((e_s - e_st).^2 ./ e_st);
    chi2_2 = sum((e_sL - e_st).^2 ./ e_st);
    f.chi2_1 = chi2_1;
    f.chi2_2 = chi2_2;
    
    SSE1 = sum( (e_s - e_st).^2 );
    SST1 = sum( (e_s - mean(e_st)).^2 );
    R_square1 = 1 - SSE1/SST1;
    f.R1 = R_square1;
    SSE2 = sum( (e_sL - e_st).^2 );
    SST2 = sum( (e_s - mean(e_st)).^2 );
    R_square2 = 1 - SSE2/SST2;
    f.R2 = R_square2;
    
    t_i = t_i + 273.15;
    e_si = e_s0*exp( (L_iv0/R_v) * (1/T_0 - 1./(t_i)) );
    e_sLi = e_s0*exp( (L_iv0 + T_0*(c_pi - c_pv))/R_v * (1/T_0 - 1./t_i) - (c_pi-c_pv)/R_v *log(t_i/T_0));
    SSE3 = sum( (e_si - e_i/100).^2 );
    SST3 = sum( (e_si - mean(e_i/100)).^2 );
    R_square3 = 1 - SSE3/SST3;
    f.R3 = R_square3;
    
    SSE4 = sum( (e_sLi - e_i/100).^2 );
    SST4 = sum( (e_sLi - mean(e_i/100)).^2 );
    R_square4 = 1 - SSE4/SST4;
    f.R4 = R_square4;

%%Plotting
    close all;

    plot(t,e_sL,'-k','LineWidth',0.75);
    hold on;
    plot(t,e_s,'-b','LineWidth',0.75);
    plot(t,e_st,'--r','LineWidth',2);
    xlim([-40 40]);
    %ylim([0 80]);
    title('Saturation Vapour Pressure for Various Latent Heat Models');
    xlabel('Temperature [^oC]');
    ylabel('Saturation Vapour Pressure [hPa]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"Lin.L(T)","Const.L(T)","Magnus-Tentens"},'location','northwest');
    set(gcf,'position',[50,550,550,300]);
    
    figure;
    plot(t,e_sL,'-k','LineWidth',0.75);
    hold on;
    plot(t,e_s,'-b','LineWidth',0.75);
    plot(t,e_st,'--r','LineWidth',2);
    xlim([7 7.05]);
    ylim([10.01 10.07]);
    title('Zoom');
    %xlabel('Temperature [^oC]');
    %ylabel('Saturation Vapour Pressure [hPa]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    %legend({"Lin.L(T)","Const.L(T)","Magnus-Tentens"},'location','northwest');
    set(gcf,'position',[750,550,600,360]);
    close;

    %{
    figure;
    plot(ti,e_si,':r','LineWidth',4);
    hold on;
    plot(ti,e_sLi,'-k','LineWidth',1);
    plot(t_i,e_i/100,'ob','LineWidth',0.75,'MarkerSize',8);
    xlim([-40 0]);
    ylim([0 6.5]);
    title('Saturation Vapour Pressure Over Ice');
    xlabel('Temperature [^oC]');
    ylabel('Saturation Vapour Pressure [hPa]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"Const.L(T)","Lin.L(T)","Real data"},'location','northwest');
    set(gcf,'position',[50,550,550,300]);
    %}
    
    
end