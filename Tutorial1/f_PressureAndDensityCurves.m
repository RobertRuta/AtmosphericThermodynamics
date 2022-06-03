function S = PressureAndDensityCurves()
%Constants
    G = [10,6,4];
    T = [30,10,-10] + 273.15;
    p_o = 1000; %hPa
    g = 9.8; %m*s^-2
    R = 8.3145/0.02897;

%Data Generation
    linstyle = {'-','--','-.'};
    colour = {'r','b','k'};
    z = linspace(0,20,2000);
    graphShift = 0;
    close all;
    for i = 1:numel(T)
        for j = 1:numel(G)
            T_z{i,j} = T(i)-G(j).*z;
            P{i,j} = p_o*(1-(G(j)/T(i)).*z).^(g*1000/(R*G(j)));
            Z_tr{i,j} = ( 1-(200/p_o)^(R*G(j)/(g*1000)) )*T(i)/G(j);
            RHO{i,j} = 100*(p_o/(R*T(i)))*(1-(G(j)/T(i)).*z).^(g*1000/(R*G(j)) - 1);
        end 
    end

    for i = 1:numel(T)
        fig = figure;
        for j = 1:numel(G)          
            plot(P{i,j},z,strcat(linstyle{j},colour{j}),'linewidth',1.5);
            hold on;
        end
        %Pressure curve plotting parameters
        xlim([200 1000]);
        ylim([0 14]);
        title('Pressure Curve at T(z=0) = ' + string(T(i)-273.15) + '^oC');
        xlabel('Pressure p [hPa]');
        ylabel('Altitude z [km]');
        set(gca,'XMinorTick','on','YMinorTick','on');    
        legend({"\Gamma_1 = 1 [K/100m]" + newline + "z_{tr1} = " + Z_tr{i,1} + " [km]","\Gamma_2 = 0.6 [K/100m]" + newline + "z_{tr2} = " + Z_tr{i,2} + " [km]","\Gamma_3 = 0.4 [K/100m]" + newline + "z_{tr3} = " + Z_tr{i,3} + " [km]"});
        set(gcf,'position',[50+graphShift,50,300,400]);
        graphShift = 400 + graphShift;
    end
    
    graphShift = 0;

    for i = 1:numel(T)
        figure;
        for j = 1:numel(G)            
            plot(RHO{i,j},z,strcat(linstyle{j},colour{j}),'linewidth',1.5);
            hold on;
        end
        %Density curve plotting parameters
        xlim([0 1.5]);
        title('Density Curve at T(z=0) = ' + string(T(i)-273.15) + '^oC');
        xlabel('Density \rho [kg/m^3]');
        ylabel('Altitude z [km]');
        set(gca,'XMinorTick','on','YMinorTick','on');    
        legend({"\Gamma_1 = 1 [K/100m]","\Gamma_2 = 0.6 [K/100m]","\Gamma_3 = 0.4 [K/100m]"});
        set(gcf,'position',[50+graphShift,550,300,400]);
        graphShift = 400 + graphShift;
    end  
    
end
    
    
    
