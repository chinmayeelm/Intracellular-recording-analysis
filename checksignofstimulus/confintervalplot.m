function confintervalplot(time,Npatterns, meanMovement,sd, CIval)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% https://in.mathworks.com/matlabcentral/answers/414039-plot-confidence-interval-of-a-signal?s_tid=answers_rc1-2_p2_MLT

                                   
ySEM = sd/sqrt(Npatterns);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI = tinv([(1-CIval)/2 (1+CIval)/2], Npatterns-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI = bsxfun(@times, ySEM, CI(:)); 

curve1 = meanMovement + yCI(1,:);
curve2 = meanMovement + yCI(2,:);
x2 = [time fliplr(time)];
inBetween = [curve1 fliplr(curve2)];
figure;
fill(x2, inBetween, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
hold on;
plot(time, meanMovement, 'k', 'LineWidth', 0.5);
xlabel('time (ms)');
ylabel('Antennal Movement (deg)');
title('STA with 95% confidance interval');

end

