function [dfti_exact,df_dt_vsTime] = dforce_dt_timeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durat,onset);

Xtemp = [0;
         0;
         rampuptime2base;
         rampuptime2base;
         onset;
         onset;
         onset + rampuptime2top;
         onset + rampuptime2top;
         onset + rampuptime2top + durat;
         onset + rampuptime2top + durat;
         onset + rampuptime2top + durat + rampdowntime2base;
         onset + rampuptime2top + durat + rampdowntime2base;
         totalduration - rampdowntime2end;
         totalduration - rampdowntime2end;
         totalduration;
         totalduration];
Ytemp = [0;
         (baselevel - 0)/rampuptime2base;
         (baselevel - 0)/rampuptime2base;
         0;
         0;
         (toplevel - baselevel)/rampuptime2top;
         (toplevel - baselevel)/rampuptime2top;
         0;
         0;
         (baselevel - toplevel)/rampdowntime2base;
         (baselevel - toplevel)/rampdowntime2base;
         0;
         0;
         (0 - baselevel)/rampdowntime2end;
         (0 - baselevel)/rampdowntime2end;
         0];


% dfti_exact = ((Ytemp(2) - Ytemp(1))*(Xtemp(2) - Xtemp(1)) + ...   % Triangle/Ramp-Up
%             ((Ytemp(6) - Ytemp(5))*(Xtemp(6) - Xtemp(5))) + ...   % Trapezoid/Ramp-Up
%             ((Ytemp(10) - Ytemp(9))*(Xtemp(10) - Xtemp(9))) + ... % Trapezoid/Ramp-Down
%             ((Ytemp(14) - Ytemp(13))*(Xtemp(14) - Xtemp(13))));   % Triangle to End

% dfti exact should be zero
dfti_exact = (Ytemp(2) - Ytemp(1))*(Xtemp(3) - Xtemp(2)) + ...  % Triangle/Ramp-Up
             (Ytemp(6) - Ytemp(5))*(Xtemp(7) - Xtemp(6)) + ...  % Trapezoid/Ramp-Up
             (Ytemp(10) - Ytemp(9))*(Xtemp(11) - Xtemp(10)) + ... % Trapezoid/Ramp-Down
             (Ytemp(14) - Ytemp(13))*(Xtemp(15) - Xtemp(14)); % Triangle to End 


if length(unique(Xtemp)) == length(Xtemp)
    tempVec = [interp1([Xtemp(3) Xtemp(2)],[Ytemp(3) Ytemp(2)],Xtemp(2):dt:Xtemp(3)),...
               interp1([Xtemp(5) Xtemp(4)],[Ytemp(5) Ytemp(4)],Xtemp(4):dt:Xtemp(5)),...
               interp1([Xtemp(7) Xtemp(6)],[Ytemp(7) Ytemp(6)],Xtemp(6):dt:Xtemp(7)),...
               interp1([Xtemp(9) Xtemp(8)],[Ytemp(9) Ytemp(8)],Xtemp(8):dt:Xtemp(9)),...
               interp1([Xtemp(11) Xtemp(10)],[Ytemp(11) Ytemp(10)],Xtemp(10):dt:Xtemp(11)),...
               interp1([Xtemp(13) Xtemp(12)],[Ytemp(13) Ytemp(12)],Xtemp(12):dt:Xtemp(13)),...
               interp1([Xtemp(15) Xtemp(14)],[Ytemp(15) Ytemp(14)],Xtemp(14):dt:Xtemp(15)),...
               ];

   [~,indTemp,~] = unique([Xtemp(2):dt:Xtemp(3),...
                Xtemp(4):dt:Xtemp(5),Xtemp(6):dt:Xtemp(7),...
                Xtemp(8):dt:Xtemp(9),...
                Xtemp(10):dt:Xtemp(11),Xtemp(12):dt:Xtemp(13),...
                Xtemp(14):dt:Xtemp(15)]);
else
    tempVec = [interp1([Xtemp(3) Xtemp(2)],[Ytemp(3) Ytemp(2)],Xtemp(2):dt:Xtemp(3)),...
               interp1([Xtemp(5) Xtemp(4)],[Ytemp(5) Ytemp(4)],Xtemp(4):dt:Xtemp(5)),...
               interp1([Xtemp(7) Xtemp(6)],[Ytemp(7) Ytemp(6)],Xtemp(6):dt:Xtemp(7)),...
               interp1([Xtemp(9) Xtemp(8)],[Ytemp(9) Ytemp(8)],Xtemp(8):dt:Xtemp(9)),...
               interp1([Xtemp(11) Xtemp(10)],[Ytemp(11) Ytemp(10)],Xtemp(10):dt:Xtemp(11)),...
               interp1([Xtemp(13) Xtemp(12)],[Ytemp(13) Ytemp(12)],Xtemp(12):dt:Xtemp(13)),...
               interp1([Xtemp(15) Xtemp(14)],[Ytemp(15) Ytemp(14)],Xtemp(14):dt:Xtemp(15)),...
               ];
       
   [~,indTemp,~] = unique([Xtemp(2):dt:Xtemp(3),...
                Xtemp(4):dt:Xtemp(5),Xtemp(6):dt:Xtemp(7),...
                Xtemp(8):dt:Xtemp(9),...
                Xtemp(10):dt:Xtemp(11),Xtemp(12):dt:Xtemp(13),...
                Xtemp(14):dt:Xtemp(15)]);
end
       

df_dt_vsTime = tempVec(indTemp);