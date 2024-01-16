function [fti_exact,f_vsTime] = forcetimeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durat,onset);




Xtemp = [0;
         rampuptime2base;
         onset;
         onset + rampuptime2top;
         onset + rampuptime2top + durat;
         onset + rampuptime2top + durat + rampdowntime2base;
         totalduration - rampdowntime2end;
         totalduration];
Ytemp = [0;
         baselevel;
         baselevel;
         toplevel;
         toplevel;
         baselevel;
         baselevel;
         0];


fti_exact = (0.5*(Ytemp(2) - Ytemp(1))*(Xtemp(2) - Xtemp(1)) + ...  % Triangle/Ramp-Up
                Ytemp(3)*(Xtemp(3) - Xtemp(2)) + ...                  % Rectangle to Onset
            (0.5*(Ytemp(3) + Ytemp(4))*(Xtemp(4) - Xtemp(3))) + ...  % Trapezoid/Ramp-Up
                (Ytemp(4)*(Xtemp(5) - Xtemp(4))) + ...                  % Rectangle of Top Level
            (0.5*(Ytemp(5) + Ytemp(6))*(Xtemp(6) - Xtemp(5))) + ... % Trapezoid/Ramp-Down
                (Ytemp(6)*(Xtemp(7) - Xtemp(6))) + ...                 % Rectangle to Rampdown
            (0.5*(Ytemp(8) - Ytemp(7)*(Xtemp(8) - Xtemp(7)))));... % Triangle to End
      
if length(unique(Xtemp)) == length(Xtemp)
    tempVec = [interp1([Xtemp(2) Xtemp(1)],[Ytemp(2) Ytemp(1)],Xtemp(1):dt:Xtemp(2)),...
               interp1([Xtemp(3) Xtemp(2)],[Ytemp(3) Ytemp(2)],Xtemp(2):dt:Xtemp(3)),...
               interp1([Xtemp(4) Xtemp(3)],[Ytemp(4) Ytemp(3)],Xtemp(3):dt:Xtemp(4)),...
               interp1([Xtemp(5) Xtemp(4)],[Ytemp(5) Ytemp(4)],Xtemp(4):dt:Xtemp(5)),...
               interp1([Xtemp(6) Xtemp(5)],[Ytemp(6) Ytemp(5)],Xtemp(5):dt:Xtemp(6)),...
               interp1([Xtemp(7) Xtemp(6)],[Ytemp(7) Ytemp(6)],Xtemp(6):dt:Xtemp(7)),...
               interp1([Xtemp(8) Xtemp(7)],[Ytemp(8) Ytemp(7)],Xtemp(7):dt:Xtemp(8))];

   [~,indTemp,~] = unique([Xtemp(1):dt:Xtemp(2),Xtemp(2):dt:Xtemp(3),Xtemp(3):dt:Xtemp(4),...
                Xtemp(4):dt:Xtemp(5),Xtemp(5):dt:Xtemp(6),Xtemp(6):dt:Xtemp(7),...
                Xtemp(7):dt:Xtemp(8)]);
else
    tempVec = [interp1([Xtemp(2) Xtemp(1)],[Ytemp(2) Ytemp(1)],Xtemp(1):dt:Xtemp(2)),...
           interp1([Xtemp(3) Xtemp(2)],[Ytemp(3) Ytemp(2)],Xtemp(2):dt:Xtemp(3)),...
           interp1([Xtemp(5) Xtemp(4)],[Ytemp(5) Ytemp(4)],Xtemp(4):dt:Xtemp(5)),...
           interp1([Xtemp(7) Xtemp(6)],[Ytemp(7) Ytemp(6)],Xtemp(6):dt:Xtemp(7)),...
           interp1([Xtemp(8) Xtemp(7)],[Ytemp(8) Ytemp(7)],Xtemp(7):dt:Xtemp(8))];
       
    [~,indTemp,~] = unique([Xtemp(1):dt:Xtemp(2),Xtemp(2):dt:Xtemp(3),...
            Xtemp(4):dt:Xtemp(5),Xtemp(6):dt:Xtemp(7),...
            Xtemp(7):dt:Xtemp(8)]);
end
       

f_vsTime = tempVec(indTemp);