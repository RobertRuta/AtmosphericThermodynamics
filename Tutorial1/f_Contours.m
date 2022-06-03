function P = f_Contours()
    
    %constants
    p_o = 1000; %hPa
    g = 9.8; %[m*s^-2]
    R = 8.3145/0.02897; %[kg*m^2*K^-1*mol^-1*s-2] / [kg*mol^-1)] = [m^2/(K*s^2)]
    k = 0.287;
    th = [300,320,340,360,380,400]; %Potential Temperature [K]
    G = 6; %lapse rate [K/km]

    z = linspace(0,30,1000);
    lats = linspace(-90,90,1000); %deg
    X = lats;   %deg
    T_o = (20*cos(2*X*2*pi/360) + 10) + 273.15; %K
    A = G./T_o; % 1/km
    Y = z.';  %km    
    B = Y*A;

    T = bsxfun(@minus,T_o,G*z.');
    P = p_o*(1-B).^(g*1000/(R*G));    
    TH = T.*(p_o ./ P).^k;
    
%%Plotting
    close all;
    %Temperature Isolines
    figure;
    %v = [300,320,340,360,380,400];
    [C,h] = contourf(X,Y,T-273.15);
    clabel(C,h);
    
    xlim([-90 90]);
    ylim([0 25]);
    title('Temperature Isolines [^oC] for \Gamma = 6 [K/km]');
    xlabel('Latitude \phi^o');
    ylabel('Altitude z [km]');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gcf,'position',[600,200,500,280]);

    %Pressure Isolines
    figure;
    [C,h] = contourf(X,Y,P);
    v = [200,400,600,800];
    clabel(C,h,v);
    
    xlim([-90 90]);
    ylim([0 18]);
    title('Pressure Isolines [hPa] for \Gamma = 6 [K/km]');
    xlabel('Latitude \phi^o');
    ylabel('Altitude z [km]');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gcf,'position',[200,200,500,280]);
    
    %Potential Temperature Isolines
    figure;
    v = [300,320,340,360,380,400,420];
    [C,h] = contourf(X,Y,TH,v);
    clabel(C,h,v);
    
    xlim([-90 90]);
    ylim([0 25]);
    title('Potential Temperature Isolines [K] for \Gamma = 6 [K/km]');
    xlabel('Latitude \phi^o');
    ylabel('Altitude z [km]');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gcf,'position',[200,200,500,280]);
end