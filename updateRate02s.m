function [rate] = updateRates(nStates,rateCoeff,ligDependentRate,ligAmp,voltage,chargeRateAmp,chargeRateSign,constant,stepNumber)

for num = 1:nStates
        for num1 = 1:nStates
            if (ligDependentRate(num, num1) > 0)
                    rate(num,num1) = rateCoeff(num,num1) * ligAmp(stepNumber);
            else
                rate(num,num1) = rateCoeff(num,num1) ;
            end
        end
end


for num = 1:nStates
        for num1 = 1:nStates
            rate(num,num1) = rate(num,num1) * exp(voltage * constant * chargeRateAmp(num,num1) * chargeRateSign(num,num1));
        end
end

for num = 1:nStates
    rate(num,num) = 0.0;
end


for i=1:nStates
        rate(i,i)= -sum(rate(:,i));
end

end