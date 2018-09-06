function [rate] = fillRates(nInputRates, to, from, inputRates)

for i=1:nInputRates
    rate(to(i),from(i)) = inputRates(i);
end

end