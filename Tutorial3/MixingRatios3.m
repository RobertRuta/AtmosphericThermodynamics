function f = MixingRatios3()
    %%Constants
    L_lv0 = 2501;       %kJ/kg
    L_iv0 = 2834;       %kJ/kg
    c_pl = 4.218;       %kJ/kgK
    c_pi = 2.106;       %kJ/kgK
    c_pv = 1.870;       %kJ/kgK
    e_s0 = 6.112;       %hPa
    R_v = 461.5/1000;   %kJ/kgK
    ep = 0.622;

    
    lapse = [10,6];
    T_0 = [30,10,-10] + 273.15;
    p_o = 1000; %hPa
    g = 9.8e-3; %km*s^-2
    R = 8.3145/0.02897;
    
    
        %w = ep*e_s ./ (p - e_s);
    linstyle = {'-','--','-.'};
    colour = {'r','b','k'};
    graphShift = 0;
    z = linspace(0,12,1000);
    close all;
    for j = 1:numel(lapse)
        figure;
        for i = 1:numel(T_0)
            T = T_0(i) - lapse(j)*z; %K
            p = p_o*(1-lapse(j)*z/T_0(i)).^(g/(R*lapse(j)));
            %e_s = e_s0*exp( (L_lv0 + T_0(i)*(c_pl - c_pv))/R_v * (1/T_0(i) - 1./T) - (c_pl-c_pv)/R_v *log(T/T_0(i)));
            e_s =  e_s0*exp(17.27*(T-273.15)./(T-273.15+237.7));
            w{i,j} = ep*e_s ./ (p - e_s);
            q{i,j} = ep*e_s ./ (p - (1-ep)*e_s);
            
            
            plot(w{i,j}*1000,z,strcat(linstyle{i},colour{i}),'linewidth',0.75);
            hold on;
            plot(q{i,j}*1000,z,strcat(linstyle{i},colour{i}),'linewidth',0.75);
        end
        %Pressure curve plotting parameters
        %xlim([200 1000]);
        ylim([0 12]);
        title('Atmospheric Mixing Ratio and Humidity for \Gamma = ' + string(lapse(j)) + ' [K/km]');
        xlabel('Mixing Ratio [g/kg]');
        ylabel('Altitude z [km]');
        set(gca,'XMinorTick','on','YMinorTick','on');    
        legend(["T_0 = 30^oC","T_0 = 10^oC","T_0 = -10^oC"]);
        set(gcf,'position',[50+graphShift,300,500,250]);
        graphShift = 700 + graphShift;
    end

end
