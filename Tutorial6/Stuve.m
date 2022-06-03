function f = Stuve()
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
    p = linspace(100,1000,100).';  %[hPa]
    T = linspace(-90,40,1000) + 273.15; %[K]
    
    %% Results
    exner = (p/p_0).^k;
    TH = T.*exner.^-1;
    
    A = p-e_s(T);
    r = ep*bsxfun(@times,e_s(T),1./A)*1000;
    
    THe = TH.*exp(L_lv*r/1000./(c_p*T));
    P = p.*(zeros(1,numel(T))+1);
    T_iso = T.*(1+zeros(numel(p),1));
    
    z = T_0/Gamma .* (1 - exner);
    T_env = T_0 - Gamma*z;
    
    %% Plotting
    labelsize = 6;
    close all;
    r_vals = [0.1,0.2,0.6,1,2,3,5,10,20,30];
    TH_vals = linspace(250,400,10);
    THe_vals = TH_vals;
    alpha = 0.7;
    
    [C5,h5] = contourf(T-273.15,p,T_iso-273.15,'Color',[0 0 0]+1,'LineStyle','--','LineWidth',0.1,'HandleVisibility','off');
    colormap(parula);
    h = colorbar;
    caxis([-150 40]);
    hold on;
    [C4,h4] = contour(T-273.15,p,P,'Color',[0 0 0]+0,'LineStyle','--','LineWidth',0.1);    
    [C1,h1] = contour(T-273.15,p,TH,TH_vals,'Color','b','LineWidth',0.8);
    [C2,h2] = contour(T-273.15,p,THe,THe_vals,'Color','c','LineStyle','-','LineWidth',0.8);
    
    [C3,h3] = contour(T-273.15,p,r,r_vals,'Color','r','LineStyle','-','LineWidth',0.8);
    
    clabel(C3,h3,'FontSize',labelsize,'Color','r');
    clabel(C2,h2,'FontSize',labelsize,'Color','c')
    
    legend('p [hPa]','\theta [K]','\theta_e [K]','r_s [g/kg]','location','southwest','NumColumns',1);
    set(get(h,'label'),'string','Temperature [^oC]');
    
    set(gca,'YDir','reverse');
    
    xlabel("Temperature T [^oC]");
    ylabel("Pressure p [hPa]")
    
    %% Nested functions
    function e_s = e_s(T)
        t = T - 273.15;
        e_s = e_s0*exp(17.27*t./(t+237.7));
    end


end