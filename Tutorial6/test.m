function test = test()
    a = -45;
    [x,y] = meshgrid(linspace(-4,4,30));
    close all;
    z = exp(-x.^2/15-y.^2);
    contour(x,y,z)
    xlim([-5 5])
    ylim([-5 5])
    x = x*cosd(a) - y*sind(a);
    y = y*cosd(a) + x*sind(a);
    figure
    contour(x,y,z)
    xlim([-5 5])
    ylim([-5 5])

end