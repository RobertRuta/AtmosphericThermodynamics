function y = FindLCL()
    
    %% Constants
    c_p = 1;            %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R = 0.2869;         %[kJ/kgK]
    R_v = 0.461;         %[kJ/kgK]
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    ep = 0.622;
    p = 1000;           %hPa
    G = 10;              %K/km
    g = 9.81/1000;       %km/s^2
    %% Ground conditions
    T_vec = linspace(0,40,100) + 273.15;    %[K]
    f_vec = linspace(1,100,100) / 100;
    T_d = 1./ (1./T_vec.' - R_v/L * log(f_vec));%[K]
    
    %% Disgusting non-vectorised calculation of LCL temp
    for i = 1:numel(T_vec)
        T = T_vec(i);
        for j = 1:numel(f_vec)        
            f = f_vec(j);
            T_LCL(i,j) = newt_raph(T,f);
        end
        p_LCL(i,:) = p*(T_LCL(i,:)/T).^(c_p/R);
    end
    
    z_LCL_2 = (T_0/G)*(1-(p_LCL/p).^(R*G/(g*1000)));    
    T_Bolton = 1./( 1./(T_vec.' - 55) - log(f_vec)/2840) + 55;
    
    %% Calculating height of LCL
    z_LCL = -(T_LCL - T_vec.') / G;
    z_LCL_Bolton = -(T_Bolton - T_vec.') / G;    
    Z = 120*(T_vec.'-T_d);
    
    absErr = abs(Z - z_LCL)*1000;
    
    %% Plotting    
    v = [1000,40,20,0,-20,-40];
    %v2 = [0,1,10,25,50,100,200];
    close all;
    [C,h] = contourf(f_vec*100,T_vec-273.15,T_LCL-273.15,v);
    clabel(C,h,v);
    
    figure;
    [C,h] = contour(f_vec*100,T_vec-273.15,absErr);
    clabel(C,h);
    xlabel("Relative Humidity [%]")
    ylabel("Temperature [^oC]")
    title("Difference [m]")
    set(gcf,'position',[50 100 600 350])
    
    figure;
    [C,h] = contourf(f_vec*100,T_vec-273.15,T_Bolton-273.15,v);
    clabel(C,h,v);

    
    %% Nested functions
    function F = F(x,T_o,f_0)
        F = log(f_0) + 1/R * (c_p * log(x/T_o) + ep*L * (1/x - 1/T_o));
    end

    function F_prim = F_prim(x)
        F_prim = 1/R * (c_p/x - ep*L/(x^2));
    end

    function X = newt_raph(T,f)
        x_old = T;
        err = 1;
        tol = abs(F(x_old,T,f))*1e-8;
        
        while err > tol
            x_new = x_old - F(x_old,T,f)/F_prim(x_old);
            x_old = x_new;
            
            err = abs(F(x_old,T,f));
        end
        X = x_old;
    end
end
