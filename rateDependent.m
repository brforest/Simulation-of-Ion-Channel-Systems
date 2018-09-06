function [ligDependentRate, chargeRateAmp, chargeRateSign] = rateDependent(nStates, nInputRates, inputRates)
to = model.rateCoeff(:,2);
from = model.rateCoeff(:,1);

ligDependentRate = zeros(nStates,nStates);

for i=1:nInputRates
    ligDependentRate(to(i),from(i)) = inputRates(i,4);
end

chargeRateAmp = zeros(nStates,nStates);

for i=1:nInputRates
    chargeRateAmp(to(i),from(i)) = inputRates(i,5);
end

chargeRateSign = zeros(nStates,nStates);

for i=1:nInputRates
    chargeRateSign(to(i),from(i)) = inputRates(i,6);
end


end