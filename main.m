function main()
format long
[model(1).rateCoeffList,model(1).stateConductance,model(1).nStates,model(1).nParams] = incomingModel('Nav_model.txt');
[model(2).rateCoeffList,model(2).stateConductance,model(2).nStates,model(2).nParams] = incomingModel('Kv_model.txt');
incomingProtocol = incoming('protocol.txt');
incomingParameters = incoming('parameters.txt');

%conductance = 20-50; %in mS/cm^2 %For the squid giant axon 0.01-0.6 mS/cm^2
%constant = 1/25; %in 1/mV       %F/RT or e-/kT at room temperature

protocol.dT = incomingProtocol(1,1); %in ms        %delta time
protocol.pulseDuration = incomingProtocol(2,:); %in ms    %reads in the duration of each individual step
protocol.ligand = incomingProtocol(3,:); %the amplitude that the ligand affects rate(s)
protocol.externalCurrent = incomingProtocol(4,:); %the external current applied to the system
protocol.dExternalCurrent= incomingProtocol(5,:); %the interval of change in external current applied during different simulations
protocol.nSim = incomingProtocol(6,1); %the number of simulations performed (in this case, each with a different current in one of the steps)
[~,protocol.nSteps] = size(protocol.pulseDuration); %calculates the number of steps through the number of columns in the file
protocol.nPointsPerStep = round(protocol.pulseDuration/protocol.dT);
protocol.nPoints = sum(protocol.nPointsPerStep);

parameters.Gleak = incomingParameters(1,1);
parameters.GKir = incomingParameters(2,1);
parameters.Capacitance = incomingParameters(3,1);
parameters.iniPotential = incomingParameters(4,1);
parameters.dtCap = protocol.dT/parameters.Capacitance;

for i=1:length(model)
    model(i).rateCoeff = fillRates(model(i).nParams, model(i).rateCoeffList(:,2), model(i).rateCoeffList(:,1), model(i).rateCoeffList(:,3)) * protocol.dT; %calls a function to fill the rate matrix from the data files
    [model(i).ligDependentRate, model(i).chargeRateAmp, model(i).chargeRateSign] = rateDependent02(model(i));
    model(i).currentRates = zeros(model(i).nStates,model(i).nStates);
    model(i).previousRates = zeros(model(i).nStates,model(i).nStates);
end


[current, voltage, channel, Ileak, time] = calcResultRK01(model,protocol,parameters); %calls the function to calculate the current, voltage, population fractions, and time
%[current, voltage, channel, time] = calcResult04(model,protocol,parameters); %calls the function to calculate the current, voltage, population fractions, and time

figure(1) %plot out the voltage, current, and population fraction with respect to time
subplot(2,1,1)
plot (time, voltage);
title('Voltage');

subplot(2,1,2)
plot (time, current);
title('Current');

figure(2);
for iii=1:protocol.nSim
    toBePlotted1(:,iii) = channel(1).states(:,5,iii);
end
plot(time, toBePlotted1);
title('fraction of open Nav');
outgoing(protocol.nPoints, protocol.nSim, time, toBePlotted1, 'open_Nav.txt')

figure(3);
for iii=1:protocol.nSim
    toBePlotted2(:,iii) = (channel(1).states(:,6,iii)+channel(1).states(:,7,iii)+channel(1).states(:,8,iii));
end
plot(time, toBePlotted2);
title('fraction of inactivated Nav');
outgoing(protocol.nPoints, protocol.nSim, time, toBePlotted2, 'inactivated_Nav.txt')

figure(4);
for iii=1:protocol.nSim
    toBePlotted2(:,iii) = channel(2).states(:,6,iii);
end
plot(time, toBePlotted2);
title('fraction of open Kv');
outgoing(protocol.nPoints, protocol.nSim, time, toBePlotted2, 'open_Kv.txt')

figure(5);
for iii=1:protocol.nSim
    toBePlotted2(:,iii) = channel(1).current(:,iii);
end
plot(time, toBePlotted2);
title('Na current');
outgoing(protocol.nPoints, protocol.nSim, time, toBePlotted2, 'Na_current.txt')

figure(6);
for iii=1:protocol.nSim
    toBePlotted2(:,iii) = channel(2).current(:,iii);
end
plot(time, toBePlotted2);
title('K current');
outgoing(protocol.nPoints, protocol.nSim, time, toBePlotted2, 'K_current.txt')

figure(7);
for iii=1:protocol.nSim
    toBePlotted2(:,iii) = Ileak(:,iii);
end
plot(time, toBePlotted2);
title('Leak current');
outgoing(protocol.nPoints, protocol.nSim, time, toBePlotted2, 'Leak_current.txt')

outgoing(protocol.nPoints, protocol.nSim, time, current, 'current.txt') %output the data to .txt files
outgoing(protocol.nPoints, protocol.nSim, time, voltage, 'voltage.txt')
end
