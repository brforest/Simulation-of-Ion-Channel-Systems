function [ligDependentRate, chargeRateAmp, chargeRateSign] = rateDependent(model)
inputRates = model.rateCoeffList;
nInputRates = model.nParams;
nStates = model.nStates;

ligDependentRate = zeros(nStates,nStates);
to = inputRates(:,2);
from = inputRates(:,1);

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