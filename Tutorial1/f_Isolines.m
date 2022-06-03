function f = Isolines()
p_o = 1000;
g = 9.8; %[m*s^-2]
R = 8.3145/0.02897; %[kg*m^2*K^-1*mol^-1*s-2] / [kg*mol^-1)] = [m^2/(K*s^2)]
k = 0.287;
th = [300,320,340,360,380,400]; %Potential Temperature [K]
G = [6,4]; %lapse rate [K/km]
P = [200,300,400,600,800];

T = f_TemperatureCurve +273.15; %Groundtemp(latitude) [K]
lats = linspace(-90,90,1800);   %Generating latitude domain

close all;
plot(lats,T-273.15,'k','linewidth',1.5)
xlim([-90 90]);
ylim([-15 35])
title('Surface Temperature Curve along Latitude Line');
xlabel('Latitude \phi^o');
ylabel('Temperature T_o [^oC]');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gcf,'position',[200,200,500,300]);

%%Pressure isolines
figure;
for i = 1:numel(P)
    for j = 1:numel(G)
        z_isoP{i,j} = T/G(j) .* ( 1 - (P(i)/p_o).^(R*G(j)/(g*1000)) );
        y = z_isoP{i,j};
        idx2 = find(lats>65,1);
        text_x2(i) = 65;
        text_y2(i) = y(idx2) + 0.5;
        labels2{i} = strcat(string(P(i)),'hPa');
        plot(lats,z_isoP{i,j},'k','linewidth',1.5);
        hold on;
    end
end
labelingP = text(text_x2, text_y2, labels2, 'Color', 'k');
xlim([-90 90]);
ylim([0 15]);
title('Pressure Isolines along Latitude Lines');
xlabel('Latitude \phi^o');
ylabel('Altitude z [km]');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gcf,'position',[200,200,500,300]);

%%Temperature isolines
n = 9;
T_n = linspace(-50,30,n) + 273.15;
for j = 1:numel(G)
    figure;
    for i = 1:n
        z_isoT{i,j} = 1/G(j) * (T - T_n(i));
        y = z_isoT{i,j};
        idx3 = find(lats>-5,1);
        text_x3(i) = -5;
        shift = [0.8, 1.2];
        text_y3(i) = y(idx3) + shift(j);
        labels3{i} = strcat(string(floor(T_n(i) - 273.15)),'^oC');
        plot(lats,z_isoT{i,j},'k','linewidth',1.5);
        hold on;
    end
    labelingT = text(text_x3, text_y3, labels3, 'Color', 'k');
    xlim([-90 90]);
    celing = [15, 22]
    ylim([0 celing(j)]);
    title({['Temperature Isolines along Latitude Lines']
        ['\Gamma = ' + string((G(j)/10)) + ' [K/100m]']});
    xlabel('Latitude \phi^o');
    ylabel('Altitude z [km]');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gcf,'position',[200,200,500,350]);
end


%%Potential Temperature isolines
figure;
for i = 1:numel(th)
    for j = 1:numel(G)
        z_isoTH{i,j} = T/G(j) .* ( 1 - (T/th(i)).^(G(j)/g) );
        y = z_isoTH{i,j};
        idx1 = find(lats>-5,1);
        text_x1(i) = -5;
        text_y1(i) = y(idx1) + 0.5;
        labels1{i} = strcat(string(th(i)),'K');
        plot(lats,z_isoTH{i,j},'k','linewidth',1.5);
        hold on;
    end
end
labelingTH = text(text_x1, text_y1, labels1, 'Color', 'k');
xlim([-90 90]);
ylim([-1 11]);
title('Potential Temperature Isolines along Latitude Lines');
xlabel('Latitude \phi^o');
ylabel('Altitude z [km]');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gcf,'position',[1000,200,500,350]);

end



