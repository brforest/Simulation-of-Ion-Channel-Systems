function [totalCurrent, voltage, channel, time] = calcResult(model, protocol, capacitance, constant)
totalCurrent = zeros(protocol.numPoints,protocol.numSim); %creates a matrix for the current information
%channel = zeros(length(model));
reversal = [50 -90];
for i=1:length(model)
    channel(i).current = zeros(protocol.numPoints,protocol.numSim); %creates a matrix for the current information
    channel(i).statesMatrix = zeros(protocol.numPoints,model(i).numStates,protocol.numSim); %creates a matrix for the states information
    channel(i).numParams = model(i).numParams; %calculates the number of rates by the number of rows in the model.txt file
    channel(i).nStates = model(i).numStates;
    %channel(i).from = int32(model(i).rateCoeff(:,1)); %the from (i.e. state that the rate begins with) value (e.g. alpha1 would have a from of 1, beta1 would have a from of 2)
    %channel(i).to = int32(model(i).rateCoeff(:,2)); %the to (i.e. state that the rate goes to) value (e.g. alpha1 would have a to of 2, beta1 would have a to of 1)
    channel(i).state = zeros(protocol.numPoints,model(i).numStates); %creates a matrix that stores the population fraction of each state
    channel(i).rate = zeros(model(i).numStates,model(i).numStates); %creates a rate matrix, whose columns will be summed after each step and simulation
    %channel(i).state = fillStates(model(i).numStates, initialStateValues); %calls a function to fill in the initial state values
    channel(i).state = zeros(model(i).numStates);
    channel(i).rateCoeff = fillRates(model(i).numParams, model(i).rateCoeff(:,2), model(i).rateCoeff(:,1), model(i).rateCoeff(:,3)) * protocol.dT; %calls a function to fill the rate matrix from the data files
    [channel(i).ligDependentRate, channel(i).chargeRateAmp, channel(i).chargeRateSign] = rateDependent(model(i).numStates, model(i).numParams, model(i).rateCoeff(:,2), model(i).rateCoeff(:,1), model(i).rateCoeff);
    channel(i).reversalPotential = reversal(i);
end
        channel(1).state(1,1) = 1.00;
        channel(2).state(1,1) = 0.900;
        channel(2).state(1,6) = 0.100;

voltage = zeros(protocol.numPoints,protocol.numSim); %creates a matrix for the voltage information
time = zeros(protocol.numPoints,1);
for i = 2:protocol.numPoints
   time(i,1) = time((i-1),1) + protocol.dT; 
end

%jj = simulation index, sss= step index, iii = number of points (nPoints),
%j = number of states index, i = second number of states index (used for
%defining delta state (dState), see below

    for jj=1:protocol.numSim            %loops through the number of simulations defined in the input text
        voltage(1,jj) = protocol.iniPotential;
        startAt = 0;                    %sets value to start at for each step in the simulation - since each step resets state, rate, etc. matrices, must start at the correct time
        for k=1:length(model)
            channel(k).dState = zeros(1,model(k).numStates);    %resets the delta State matrix
            channel(k).newState = channel(k).state;               %resets the newState matrix, which is created so that we can reset the original value for each simulation
            channel(k).newRates = zeros(model(k).numStates,model(k).numStates);
        end      
                
        for sss = 1:protocol.numSteps                 %loops through the number of steps in the simulation
            for iii = 2:protocol.nPointsPerStep(sss)  %loops through until the end of the step is reached (nPoints is reached)
               for nnn=1:length(model)     
                  for j = 1:channel(nnn).nStates              %goes through all of the states for each point
                     channel(nnn).dState(1,j) = 0.0;                  %resets the dState matrix
                     channel(nnn).newRates = updateRates(channel(nnn).nStates,channel(nnn).rateCoeff,channel(nnn).ligDependentRate,protocol.ligAmp,voltage(iii+startAt-1,jj),channel(nnn).chargeRateAmp,channel(nnn).chargeRateSign,constant,sss); %calls a function to take the rate matrix, apply any charge, current, or ligand dependencies, and perform the sum operation in order to make the matrix useful for the calculation
                     for i = 1:model(nnn).numStates           %loops though the number of states, calculates the change in each of the state to calculate each dState value
                        channel(nnn).dState(1,j) = channel(nnn).dState(1,j)+channel(nnn).newState((startAt+iii-1),i)*channel(nnn).newRates(j,i); %adds to the dState value for each state
                     end
                     channel(nnn).newState((iii+startAt),j) = channel(nnn).newState((startAt+iii-1),j) + channel(nnn).dState(1,j); %adds dState to the value of the previous value in the matrix
                     channel(nnn).current((iii+startAt),jj) = channel(nnn).current((iii+startAt),jj) + (channel(nnn).newState((iii+startAt),j) * model(nnn).stateConductance(1,j) * (voltage((iii+startAt-1),jj)-channel(nnn).reversalPotential)); %calculates the current, summing together the concentration of each state times its conductance, then time the current voltage minus the resting voltage
                  end
                  totalCurrent((iii+startAt),jj) = totalCurrent((iii+startAt),jj)+ channel(nnn).current((iii+startAt),jj);
                  channel(nnn).statesMatrix(iii+startAt,:,jj) = channel(nnn).newState((iii+startAt),:);
               end
               totalCurrent((iii+startAt),jj) = totalCurrent((iii+startAt),jj) + (protocol.externalCurrent(sss) + protocol.currentInterval(sss)*(jj-1)); %adds the external current to current of the system
               totalCurrent((iii+startAt),jj) = totalCurrent((iii+startAt),jj) + protocol.Gleak*voltage((iii+startAt-1),jj);
               dV = (-1/capacitance) * totalCurrent((iii+startAt),jj) * protocol.dT; %calculates delta voltage based on the present current
               voltage((iii+startAt),jj) = voltage((iii+startAt-1),jj) + dV; %calculates the voltage by adding delta voltage to the previous voltage value
            end
            startAt = startAt + protocol.nPointsPerStep(sss)-1; %adds the number of points for this step, nPoints, to the startAt variable so that the next step will begin at the correct point in time
        end
    end
end
