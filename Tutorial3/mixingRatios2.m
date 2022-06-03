function f = mixingRatios2()
    
    %%Constants
    L_lv0 = 2501;       %kJ/kg
    L_iv0 = 2834;       %kJ/kg
    c_pl = 4.218;       %kJ/kgK
    c_pi = 2.106;       %kJ/kgK
    c_pv = 1.870;       %kJ/kgK
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    R_v = 461.5/1000;   %kJ/kgK
    ep = 0.622;
    
   
    %e_s = e_s0*exp( (L_lv0/R_v) * (1/T_0 - 1./T) );
    %e_sL = e_s0*exp( (L_lv0 + T_0*(c_pl - c_pv))/R_v * (1/T_0 - 1./T) - (c_pl-c_pv)/R_v *log(T/T_0));
    
    t = linspace(-40,40,100); %oC    
    e_s = e_s0*exp(17.27*t./(t+237.7)); %hPa
    p = linspace(100,1000,10).';   %hPa

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
    SSE1 = sum(eps1.^2);
        %R_square1 = SSE1 / w_simple;
    f.SSE1 = SSE1;
    
    SSE2 = sum( (w_simple - q).^2 );
    SST2 = sum( (w_simple - mean(w_simple)).^2 );
    R_square2 = 1 - SSE2/SST2;
    f.R2 = R_square2;

    %%Plotting
    close all;

    figure;
    surf(t,p,w);
    %figure;
    hold on;
    surf(t,p,q);
    %xlim([800 1000]);
    %ylim([-40 -30]);
    
    hold on;
    surf(t,p,q);
    %xlim([800 1000]);
    %ylim([-40 -30]);
    
    figure;
    surf(t,p,difference);
    %xlim([0 1000]);
    %ylim([-40 40]);
end