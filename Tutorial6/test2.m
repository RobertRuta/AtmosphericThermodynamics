function f = SkewT()
    %% Constants
    p_0 = 1000;     %hPa
    T_0 = 273.15;
    Gamma = 6;     %K/km
    k = 0.287;
    e_s0 = 6.112;   %hPa
    L_lv = 2501;    %kJ/kg
    c_p = 1;        %kJ/kg
    ep = 0.622;     %Unitless
    
    %% Domains
    p = linspace(100,1000,10).';  %[hPa]
    T = linspace(-100,40,1000) + 273.15; %[K]
    b = 10^(1/9);
    p_log = 100*b.^(0:9).';
    
    %% Results    
    ZT = 40 + T*
    P = p.*(zeros(1,numel(T))+1);
    %% Plotting
    close all;
    r_vals = [0.01,0.5,4,8,12,20]/1000;
    T_vals = linspace(-100,40,15);
    TH_vals = linspace(250,400,10);
    THe_vals = TH_vals;
    alpha = 0.7;
    
    p_vals = linspace(100,1000,10);
    [C1,h1] = contour(T-273.15,p,P,p_vals,'Color','k','LineWidth',0.8);    
    hold on;
    [C2,h2] = contour(T-273.15,p_log,ZT,'Color',[0 0 0]+alpha,'LineWidth',0.8);
    clabel(C2,h2);
    %[C3,h3] = contour(T-273.15,TH,T_iso,'Color',[0 0 0]+alpha,'LineWidth',0.8)
    %[C4,h4] = contour(T-273.15,TH,r_s,r_vals,'Color','g','LineWidth',0.5)
    %[C5,h5] = contour(T-273.15,TH,THe,THe_vals,'c--','LineWidth',0.8)
    
    set(gca, 'YDir', 'reverse')    
    set(gca, 'YScale', 'log')

    xlim([-40, 40]);
    %ylim([200,500]);
    
    %% Nested functions
    function e_s = e_s(T)
        t = T - 273.15;
        e_s = e_s0*exp(17.27*t./(t+237.7));
    end

end