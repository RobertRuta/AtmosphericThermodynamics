function f = mixingRatios()
    
    %%Constants
    L_lv0 = 2501;       %kJ/kg
    L_iv0 = 2834;       %kJ/kg
    c_pl = 4.218;       %kJ/kgK
    c_pi = 2.106;       %kJ/kgK
    c_pv = 1.870;       %kJ/kgK
    e_s0 = 6.112;       %hPa
    R_v = 461.5/1000;   %kJ/kgK
    ep = 0.622;    
    
    t = linspace(-40,40,100);           %oC    
    p = linspace(100,1000,100).';       %hPa
    e_s = e_s0*exp(17.27*t./(t+237.7)); %hPa


    %w = ep*e_s ./ (p - e_s);
        A = bsxfun(@minus,p,e_s);
        B = 1./A;
    w = ep*bsxfun(@times,e_s, B);
    %q = ep*e_s ./ (p - (1-ep)*e_s);
        C = p - (1-ep)*e_s;
        D = 1./C;
    q = ep*bsxfun(@times,e_s,D);
    
    difference = w - q;
    
%%Checking validity
    w_simple = ep*bsxfun(@times,e_s,1./p);
    
    eps1 = (w_simple - w).^2;
    SSE1 = sum(sum(eps1.^2));
        %R_square1 = SSE1 / w_simple;
    f.SSE1 = SSE1;
    
    eps2 = (w_simple - q).^2;
    SSE2 = sum(sum(eps2.^2));
        %R_square1 = SSE1 / w_simple;
    f.SSE2 = SSE2;

    %%Plotting
    close all;
    figure;
    w_vals = 100*[0,1e-3,2e-3,4e-3,8e-3,16e-3,32e-3,64e-3,128e-3,256e-3];
    [C,h] = contour(t,p,100*w,w_vals);
    clabel(C,h,w_vals,'FontSize',8,'Color','k');
    title('Vapour Mixing Ratio Contour Plot');
    xlabel('Temperature [^oC]');
    ylabel('Pressure [hPa]');
    legend({"w [%]"},'location','northwest');
    set(gca,'XMinorTick','on','YMinorTick','on');  
    set(gcf,'position',[50,550,550,200]);
    
    figure;
    q_vals = 100*[0,1e-3,2e-3,4e-3,8e-3,16e-3,32e-3,64e-3,128e-3,256e-3];
    [C,h] = contour(t,p,100*q,q_vals);
    clabel(C,h,q_vals,'FontSize',8,'Color','k');
    title('Specific Humidity Contour Plot');
    xlabel('Temperature [^oC]');
    ylabel('Pressure [hPa]');
    legend({"q [%]"},'location','northwest');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    set(gcf,'position',[750,550,550,200]);
    
    figure;
    difference_vals = 100*[0,0.00000001*logspace(1,5,5),0.01,0.1];
    [C,h] = contour(t,p,100*difference,difference_vals);
    clabel(C,h,difference_vals,'FontSize',8,'Color','k');
    title('Difference Contour Plot');
    xlabel('Temperature [^oC]');
    ylabel('Pressure [hPa]');
    legend({"d [%]"},'location','northwest');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    set(gcf,'position',[1350,550,550,200]);
end