function y = FindLCL()
    
    %% Constants
    c_p = 1;            %[kJ/kgK]
    L = 2501;           %[kJ/kgK]
    R = 0.2869;         %[kJ/kgK]
    T_0 = 273.15;       %K
    e_s0 = 6.112;       %hPa
    ep = 0.622;
    p = 1000;           %hPa
    G = 10;              %K/km
    g = 9.81/1000;       %km/s^2
    %% Ground conditions
    T_vec = [-20,-10,0,10,20] + 273.15;
    T_o = linspace(-10,30,100) + 273.15;
    f_vec = linspace(1,100,100) / 100;
    
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
    

    
    y = T_LCL;
    
    %% Calculating height of LCL
    z_LCL = -(T_LCL - T_vec.') / G;
    z_LCL_Bolton = -(T_Bolton - T_vec.') / G;
    
    %% Plotting    
    
    v = [40,20,0,-20,-40]
    close all;
    [C,h] = contourf(f_vec*100,T_vec-273.15,T_LCL-273.15,v);
    set(gca,'Ydir','reverse');
    clabel(C,h,v);
    
    figure;
    [C,h] = contourf(f_vec*100,T_vec-273.15,T_Bolton-273.15,v);
    set(gca,'Ydir','reverse');
    clabel(C,h,v);
    
  
    figure; %height of lcl
    for i = 1:numel(T_vec)
        plot(f_vec*100,z_LCL(i,:),"LineWidth",1)
        labelT(i) = strcat("T = ",string(T_vec(i) - 273.15)," [^oC]");
        hold on;
    end
    xlabel("Relative Humidity [%]")
    ylabel("LCL Height [km]")
    legend(labelT)
    title("Lifting Condensation Level Height - Method 1")
    set(gcf,'position',[50 100 600 350])
    
    figure; %height of lcl_2
    for i = 1:numel(T_vec)
        plot(f_vec*100,z_LCL_2(i,:),"LineWidth",1)
        labelT(i) = strcat("T = ",string(T_vec(i) - 273.15)," [^oC]");
        hold on;
    end
    xlabel("Relative Humidity [%]")
    ylabel("LCL Height [km]")
    legend(labelT)
    title("Lifting Condensation Level Height - Method 2")
    set(gcf,'position',[50 100 600 350])
    
    
    figure; %height of lcl from Bolton formula
    for i = 1:numel(T_vec)
        plot(f_vec*100,z_LCL_Bolton(i,:),"LineWidth",1)
        hold on;
    end
    xlabel("Relative Humidity [%]")
    ylabel("LCL Height [km]")
    legend(labelT)
    title("Bolton's Lifting Condensation Level Height")
    set(gcf,'position',[650 100 600 350])
    
    figure; %Abosolute Error Plot
    for i = 1:numel(T_vec)
        absErr = z_LCL(i,:)-z_LCL_Bolton(i,:);
        plot(f_vec*100,absErr,"LineWidth",1)
        hold on;
    end
    xlabel("Relative Humidity [%]")
    ylabel("Difference Error [km]")
    legend(labelT)
    title("Difference between model 1 and Bolton's model")
    set(gcf,'position',[1350 100 600 250])
    
    figure; %Abosolute Error Plot
    for i = 1:numel(T_vec)
        absErr = z_LCL_2(i,:)-z_LCL_Bolton(i,:);
        plot(f_vec*100,absErr,"LineWidth",1)
        hold on;
    end
    xlabel("Relative Humidity [%]")
    ylabel("Difference Error [km]")
    legend(labelT)
    title("Difference between model 2 and Bolton's model")
    set(gcf,'position',[1350 100 600 250])
    
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
