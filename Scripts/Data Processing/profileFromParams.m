function [y_px] = profileFromParams(x_px,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,topduration,toponset)

% Screen Resolution
canvasx = 1280;
canvasy = 1024;
border = 100;

% Scale Factors
timescale = (canvasx - 2*border);
forcescale = (canvasy - 2*border);

% 
zerof_y = canvasy - border;
timeend_x = canvasx - border;
timestart_x = border;

% Upper left is [0,0]

if x_px > timeend_x
    y_px = zerof_y;
%Ramp Down To End
elseif x_px > timeend_x - timescale*rampdowntime2end/totalduration
    y_px = (baselevel*forcescale)/(timescale*rampdowntime2end/totalduration)*(x_px - (timeend_x - timescale*rampdowntime2end/totalduration)) + (zerof_y - baselevel*forcescale);
%Right Base
elseif x_px > timestart_x + timescale*(toponset + rampuptime2top + topduration + rampdowntime2base)/totalduration
    y_px = zerof_y - baselevel*forcescale;
% Ramp Down To Base
elseif x_px > timestart_x + timescale*(toponset + topduration + rampuptime2top)/totalduration
    y_px = -1*(baselevel - toplevel)*forcescale/(timescale*rampdowntime2end/totalduration)*(x_px - (timescale*(toponset + rampuptime2top + topduration)/totalduration + timestart_x)) + (zerof_y - toplevel*forcescale);
% Top Phase    
elseif x_px > timestart_x + timescale*(toponset + rampuptime2top)/totalduration
    y_px = zerof_y - toplevel*forcescale;   
% Ramp Up to Top
elseif x_px > timestart_x + timescale*toponset/totalduration
    y_px = -1*(toplevel - baselevel)*forcescale/(rampuptime2base*timescale/totalduration)*(x_px - (toponset*timescale/totalduration + timestart_x)) + (zerof_y - baselevel*forcescale);
% Left Base
elseif x_px > timestart_x + timescale*rampuptime2base/totalduration
    y_px = zerof_y - baselevel*forcescale;
% Ramp Up Phase
elseif x_px > timestart_x
    y_px = -1*(baselevel*forcescale)/(timescale*rampuptime2base/totalduration)*(x_px - timestart_x) + zerof_y;
elseif x_px <= timestart_x
    y_px = zerof_y;
end