function y = SkewT()
    %% Constants
    p_0 = 1;     %bar
    T_0 = 273.15;
    Gamma = 6;     %K/km
    g = 9.81/1000;  %km/s/s
    R = 8.3145/0.02897; 
    k = 0.287;
    e_s0 = 6.112/1000;   %bar
    L_lv = 2501;    %kJ/kg
    c_p = 1;        %kJ/kg
    ep = 0.622;     %Unitless

    n = 100;
    b = 10^(1/(n-1));
    p = linspace(0.1,1,n);
    p_log = 0.1*b.^(0:n-1);
    T = linspace(-40,40,n);
    [x,Y] = meshgrid(T,p_log);
    T_mesh = x + x' - 40;
    
    TH = (T_mesh+273.15).*(p_0./Y).^k;    
    r_s = ep*e_s(T_mesh+273.15) ./ (Y - e_s(T_mesh+273.15))*1000;
    %THe = TH.*exp(L_lv.*r_s*1e-3 ./ (c_p*(T_mesh+273.15)));
    
    %% Air Parcel
    G_d = 10;   %K/km
    G = 8;      %K/km
    T_ground = 20+273.15; %K
    T_parc_0 = 20+273.15; %K
    %T_LCL = 1./( 1./(T_vec.' - 55) - log(f_vec)/2840) + 55;
    %z_LCL =     
    p_parc = Y(:,1);
    z = T_ground/G * (1 - (p_parc./p_0).^(R*G/1000000/g));
    T_parc = T_parc_0 - G_d*z(:,1);
    %T_parc_skew = T_parc +
    
    %% Plotting
    labelsize = 6;
    r_vals = [0.01,0.1,0.5,2,4,6,8,10,12,16,20,24];
    T_vals = linspace(-150,40,20);
    
    Y = Y * 1000;
    close all;
    [c,h] = contourf(x,Y,T_mesh,T_vals,'Color','w','LineWidth',0.1,'HandleVisibility','off');
    colormap(parula);
    h = colorbar;
    caxis([-150 40]);
    hold on;
    [c1,h1] = contour(x,Y,Y,':k','LineWidth',0.75,'HandleVisibility','off');
    %[c2,h2] = contour(x,Y,TH,'b','LineWidth',0.75);
    [c3,h3] = contour(x,Y,r_s,r_vals,'r','LineWidth',0.75);
    %[c4,h4] = contour(x,Y,THe,'c','LineWidth',0.75);
    set(gca,'YDir','reverse')
    %clabel(c2,h2,'FontSize',labelsize)
    clabel(c3,h3,'FontSize',labelsize)
    %clabel(c4,h4,'FontSize',labelsize)
    set(gca,'YScale','log')
    xtickangle(45);
    yticks((0:10)*100);
    title("Skew T")
    %yticklabels({'100','200','300','0','\pi','2\pi','3\pi'})
    %legend('\theta [K]','r_s [g/kg]','\theta_e [K]','location','northwest','NumColumns',3);
    set(get(h,'label'),'string','Temperature [^oC]');
    
    xlabel("Temperature T [^oC]");
    ylabel("Pressure p [hPa]")
    
    
    function y = log_b(b,x);
        y = log(x)/log(b);
    end
    function e_s = e_s(T)
        t = T - 273.15;
        e_s = e_s0*exp(17.27*t./(t+237.7));
    end



end