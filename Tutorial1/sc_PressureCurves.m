function P = PressureCurves()
%Constants
    G1 = 0.01;
    G2 = 0.006;
    G3 = 0.004;
    T1 = 30+275.15;
    T2 = 10+275.15;
    T3 = -10+275.15;
    p_o = 1000; %hPa
    g = 9.8; %km*s^-2
    R = 8.3145/0.02897;

%Data Generation
    z = linspace(0,15,1000);
    p1_1 = p_o*(1-(G1*1000/T1).*z).^(g/(R*G1));
    p1_2 = p_o*(1-(G1*1000/T2).*z).^(g/(R*G1));
    p1_3 = p_o*(1-(G1*1000/T3).*z).^(g/(R*G1));
    p2_1 = p_o*(1-(G2*1000/T1).*z).^(g/(R*G2));
    p2_2 = p_o*(1-(G2*1000/T2).*z).^(g/(R*G2));
    p2_3 = p_o*(1-(G2*1000/T3).*z).^(g/(R*G2));
    p3_1 = p_o*(1-(G3*1000/T1).*z).^(g/(R*G3));
    p3_2 = p_o*(1-(G3*1000/T2).*z).^(g/(R*G3));
    p3_3 = p_o*(1-(G3*1000/T3).*z).^(g/(R*G3));
    
    ztr1_1 = (1-(200/p_o)^(R*G1/g))*T1/(1000*G1);
    ztr1_2 = (1-(200/p_o)^(R*G1/g))*T2/(1000*G1);
    ztr1_3 = (1-(200/p_o)^(R*G1/g))*T3/(1000*G1);
    ztr2_1 = (1-(200/p_o)^(R*G2/g))*T1/(1000*G2);
    ztr2_2 = (1-(200/p_o)^(R*G2/g))*T2/(1000*G2);
    ztr2_3 = (1-(200/p_o)^(R*G2/g))*T3/(1000*G2);
    ztr3_1 = (1-(200/p_o)^(R*G3/g))*T1/(1000*G3);
    ztr3_2 = (1-(200/p_o)^(R*G3/g))*T2/(1000*G3);
    ztr3_3 = (1-(200/p_o)^(R*G3/g))*T3/(1000*G3);
                        
   
%Plotting    
    close all;
    
    plot(p1_1,z,'-r','linewidth',1.5);
    hold on;
    plot(p2_1,z,'--b','linewidth',1.5);
    plot(p3_1,z,'-.m','linewidth',1.5);
    xlim([200 1000]);
    title('Pressure Curve at T(z=0) = 30^oC');
    xlabel('Pressure p [hPa]');
    ylabel('Altitude z [km]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"\Gamma_1 = 1[K/100m]" + newline + "z_{tr_1} = " + ztr1_1,"\Gamma_2 = 0.6[K/100m]" + newline + "z_{tr2} = " + ztr2_1,"\Gamma_3 = 0.4[K/100m]" + newline + "z_{tr3} = " + ztr3_1});
    set(gcf,'position',[50,50,300,400]);
    
    figure;
    plot(p1_2,z,'-r','linewidth',1.5);
    hold on;
    plot(p2_2,z,'--b','linewidth',1.5);
    plot(p3_2,z,'-.m','linewidth',1.5);
    xlim([200 1000]);
    title('Pressure Curve at T(z=0) = 10^oC');
    xlabel('Pressure p [hPa]');
    ylabel('Altitude z [km]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"\Gamma_1 = 1[K/100m]" + newline + "z_{tr_1} = " + ztr1_2,"\Gamma_2 = 0.6[K/100m]" + newline + "z_{tr2} = " + ztr2_2,"\Gamma_3 = 0.4[K/100m]" + newline + "z_{tr3} = " + ztr3_2});
    set(gcf,'position',[450,50,300,400]);
    
    figure;
    plot(p1_3,z,'-r','linewidth',1.5);
    hold on;
    plot(p2_3,z,'--b','linewidth',1.5);
    plot(p3_3,z,'-.m','linewidth',1.5);
    xlim([200 1000]);
    title('Pressure Curve at T(z=0) = -10^oC');
    xlabel('Pressure p [hPa]');
    ylabel('Altitude z [km]');
    set(gca,'XMinorTick','on','YMinorTick','on');    
    legend({"\Gamma_1 = 1[K/100m]" + newline + "z_{tr_1} = " + ztr1_3,"\Gamma_2 = 0.6[K/100m]" + newline + "z_{tr2} = " + ztr2_3,"\Gamma_3 = 0.4[K/100m]" + newline + "z_{tr3} = " + ztr3_3});
    set(gcf,'position',[900,50,300,400]);
    
    P.p1_1 = p1_1;
    P.p1_2 = p1_2;
    P.p1_3 = p1_3;
    P.p2_1 = p2_1;
    P.p2_2 = p2_2;
    P.p2_3 = p2_3;
    P.p3_1 = p3_1;
    P.p3_2 = p3_2;
    P.p3_3 = p3_3;
    
    
end
    
    
    
