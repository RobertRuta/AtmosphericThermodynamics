function f = test3()
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
    n = 100;
    TH = linspace(200,500,n).';
    b = (TH(end)/TH(1))^(1/(n-1));
    TH_log = TH(1)*b.^(0:n-1);    
    T = linspace(-150,40,n); %[oC]
    [x,y] = meshgrid(T,TH_log);
    X = x - x.' -150;
    Y = y;
    
    %% Results    
    P = p_0*(X./Y).^(1/k);
    r_s = ep*e_s(X) ./ (P - e_s(X));
    THe = TH.*exp(L_lv*r_s./(c_p*X));
    
    %% Plotting
    close all;
    
    [c,h] = contourf(x,y,X)
    colormap(parula);
    h = colorbar;
    %caxis([-150 40]);
    hold on;
    [c1,h1] = contour(x,y,Y)
    
    set(gca, 'YScale', 'log')
    xlim([-120, 40]);
    ylim([200,500]);
    
    %% Nested functions
    function e_s = e_s(T)
        t = T - 273.15;
        e_s = e_s0*exp(17.27*t./(t+237.7));
    end
    function y = log_b(b,x);
        y = log(x)/log(b);
    end

end