function f = Tephigram()
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
    T = linspace(-150,40,1000) + 273.15; %[K]
    TH = linspace(200,500,100).';
    
    %% Results    
    P = p_0*(T./TH).^(1/k);
    TH_iso = TH .* (1+zeros(1,numel(T)));
    T_iso = T .* (1+zeros(numel(TH),1));
    A = P-e_s(T);
    r_s = ep*bsxfun(@times,e_s(T),1./A)*1000;
    THe = TH.*exp(L_lv*r_s/1000./(c_p*T));
    
    %% Plotting
    close all;
    labelsize = 6;
    r_vals = [0.01,0.5,4,8,12,20];
    TH_vals = linspace(250,400,10);
    THe_vals = TH_vals;
    alpha = 0.7;
    
    p_vals = linspace(200,1000,9);
    [C3,h3] = contourf(T-273.15,TH,T_iso-273.15,'Color','w','LineWidth',0.8,'HandleVisibility','off')
    colormap(parula);
    h = colorbar;
    caxis([-150 40]);
    hold on;
    [C2,h2] = contour(T-273.15,TH,TH_iso,'Color',[0 0 0]+1,'LineWidth',0.8,'HandleVisibility','off');
    [C5,h5] = contour(T-273.15,TH,THe,THe_vals,'c-','LineWidth',0.8);
    [C1,h1] = contour(T-273.15,TH,P,p_vals,'-k','LineWidth',0.8);    
    [C4,h4] = contour(T-273.15,TH,r_s,r_vals,'r-','LineWidth',0.5);
    
    legend('\theta_e [K]','p [hPa]','r_s [g/kg]','location','northwest','NumColumns',3);
    set(get(h,'label'),'string','Temperature [^oC]');
    
    clabel(C1,h1,p_vals,'FontSize',labelsize);
    clabel(C4,h4,'FontSize',labelsize)
    %clabel(c4,h4,'FontSize',labelsize)
    
    xlabel("Temperature T [^oC]");
    ylabel("Potential Temperature \theta [K]")
    
    
    set(gca, 'YScale', 'log')
    xlim([-120, 40]);
    ylim([200,500]);
    
    %% Nested functions
    function e_s = e_s(T)
        t = T - 273.15;
        e_s = e_s0*exp(17.27*t./(t+237.7));
    end

end