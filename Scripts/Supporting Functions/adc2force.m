function FandM = adc2force(adc1,adc2,adc3,adc4,adc5,adc6)
% Convert Robot Data to Force (N)
% TODO: What does ft_scale do? 
% TODO: Confirm Newtons?
% TODO: Investigate if bias needs to be incorporated
% Currently does not match ft_world.x and z calcs

% Combine ADCs, ensure same lengths
adc = [adc1';
       adc2';
       adc3';
       adc4';
       adc5';
       adc6'];

% Constant params for conversion
ft_cal = [[-0.41763,-0.22235,1.26797,-30.96059,-0.99539,31.98751];
          [-1.25792,35.35964,0.09228,-18.03949,1.06615,-18.27562];
          [17.48501,-0.39301,17.63390,-0.35360,17.35616,-0.69492];
          [-0.20629,-0.74400,-30.83439,0.98932,30.32068,-0.70251];
          [35.90813,-0.64760,-17.95408,-0.43155,-17.26641,1.34019];
          [0.92209,-16.87856,0.43325,-17.06802,0.86833,-17.60612]];
ft_scale = [9.05964537732276,...
            9.05964537732276,...
            2.81161408261741,...
            4.29187196412003,...
            4.29187196412003,...
            4.09280554648923]';


FandM = (ft_cal*adc).*ft_scale;
