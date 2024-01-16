function y_px = adc2ypx(adc2,adc4,adc2i,adc4i,el_ang,mvc)

canvash = 1024;
border = 100;
forcescale = canvash - 2*border;
zerof_y = canvash - border;

y_px = zerof_y - (forcescale/mvc)*(-cosd(270 - el_ang)*(adc2 - adc2i) + sind(270 - el_ang)*(adc4-adc4i));