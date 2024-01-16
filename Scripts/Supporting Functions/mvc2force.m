function F = mvc2force(mvc,el_ang)
% MVC_volts = -cosd(270 - el_ang)*(adc2 - adc2i) + sind(270 - el_ang)*(adc4-adc4i));

% TODO: Figure out how to actually convert this

MVC_F = -cosd(270 - el_ang)*(d_adc2)*A_Fx + sind(270 - el_ang)*(d_adc4)*A_Fz;