function x_px = adc2xpx(adc2,adc4,adc2i,adc4i,el_ang,mvc,x_px_pre1,x_px_pre2)
% This function is only applicable during the choice phase, where
% horizontal position of the cursor is controlled by the participant.

canvas_centerx = 1280/2;
border = 100;
% Rleftside = canvas_centerx + border;
% Lrightside = canvas_centerx - border;
% Rleftside - Lrightside == 2*border

% set cob(Rleftside) [expr { $cob(canvascenterx) + $cob(border) }]
% set cob(Lrightside) [expr { $cob(canvascenterx) - $cob(border) }]


x_px = -sind(270 - el_ang)*(adc2 - adc2i) + cosd(270 - el_ang)*(adc4-adc4i);
x_px = canvas_centerx - x_px*2*border/(0.1*mvc);
x_px = (x_px + x_px_pre1 + x_px_pre2)/3;


%             		set ob(xcursor) [expr { -1*($ob(sin_el)*([rshm adcvolts 2] - $ob(adc2pre)) + $ob(cos_el)*([rshm adcvolts 4] - $ob(adc4pre))) }]
% 		            set ob(x_px) [expr { $cob(canvascenterx) + $ob(xcursor)*($cob(Rleftside) - $cob(Lrightside))/(0.1*$mob(mvc)) } ]		

% Smoothing...
% 		            set ob(x_px) [expr { ($ob(x_px) + $ob(x_px_pre1) + $ob(x_px_pre2))/3 }]
% 		            set ob(x_px_pre2) $ob(x_px_pre1)
% 		            set ob(x_px_pre1) $ob(x_px)