function [rate] = updateRates03(model, ligand, voltage,kTconstant)

for num = 1:model.nStates
        for num1 = 1:model.nStates
            if (model.ligDependentRate(num, num1) > 0)
                    rate(num,num1) = model.rateCoeff(num,num1) * ligand;
            else
                rate(num,num1) = model.rateCoeff(num,num1) ;
            end
        end
end

for num = 1:model.nStates
        for num1 = 1:model.nStates
            rate(num,num1) = rate(num,num1) * exp(voltage * model.chargeRateAmp(num,num1) * model.chargeRateSign(num,num1)*kTconstant);
        end
end

for num = 1:model.nStates
    rate(num,num) = 0.0;
end


for i=1:model.nStates
        rate(i,i)= -sum(rate(:,i));
end

end