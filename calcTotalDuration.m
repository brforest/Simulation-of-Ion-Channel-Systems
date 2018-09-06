function [totalDuration] = calcTotalDuration(duration,nSteps)
totalDuration = 0;
for i = 1:nSteps
   totalDuration = totalDuration + duration(i);
end

end