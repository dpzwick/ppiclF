close all; clear all; clc

fig = imread('ling_2016_figure_14.png');
imshow(fig);
hold on

scale_x = 2e-3;
scale_y = scale_x/(1.66*343);

x_origin_value = 0.0*scale_x;
x_max_value = 40*scale_x;
y_origin_value = 0.0*scale_y;
y_max_value = 100*scale_y;

fprintf("Select reference origin value at (%f, %f)\n", x_origin_value/scale_x ...
                                                   , y_origin_value/scale_y);
[x_origin,y_origin] = ginput(1);
fprintf("Select reference x value at %f\n", x_max_value/scale_x);
[x_max,dum] = ginput(1);
fprintf("Select reference y value at %f\n", y_max_value/scale_y);
[dum,y_max] = ginput(1);

my_y   = [0    30   60   90   120   150   180   200  ]*10^(-6);
my_upf = [5.01 5.04 5.14 5.28 5.48  5.75  6.10  6.38 ]*10^(-3) - 5e-3;
my_dpf = [7.00 7.07 7.69 8.88 10.30 11.90 13.60 15.00]*10^(-3) - 5e-3;

xx = @(x) (x - x_origin_value)*(x_max - x_origin)/(x_max_value-x_origin_value) + x_origin;
yy = @(y) (y - y_origin_value)*(y_max - y_origin)/(y_max_value-y_origin_value) + y_origin;

plot(xx(my_upf),yy(my_y),'mx-','linewidth',2);
plot(xx(my_dpf),yy(my_y),'mx-','linewidth',2);
