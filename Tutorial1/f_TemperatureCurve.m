function T = f_TemperatureCurve()
    lats = linspace(-90,90,1800);   %Generating latitude domain
    T = 20*cos(2*lats*2*pi/360)+10; %Generating ground temperature(latitude) curve

    plot(lats,T,'-k');  %Plotting ground temp against latitude
    xlabel('latitude');
    ylabel('Temperature [^oC]');
end