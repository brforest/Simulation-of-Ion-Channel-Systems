function [totalCurrent, voltage, channel, time] = calcResult03(model, protocol, parameters)
%channel = zeros(length(model));
reversal = [50 -90];
nModels = length(model)

for i=1:nModels
    channel(i).current = zeros(protocol.nPoints,protocol.nSim); %creates a matrix for the current information
    channel(i).statesMatrix = zeros(protocol.nPoints,model(i).nStates,protocol.nSim); %creates a matrix for the states information
    channel(i).nParams = model(i).nParams; %calculates the number of rates by the number of rows in the model.txt file
    channel(i).states = zeros(protocol.nPoints,model(i).nStates,protocol.nSim); %creates a matrix that stores the population fraction of each state
    channel(i).dStates = zeros(1,model(i).nStates);
    channel(i).reversalPotential = reversal(1,i);
end

channel(1).states(1,1,:) = 0.9900;
channel(1).states(1,2,:) = 0.0100;
channel(2).states(1,1,:) = 0.9700;
channel(2).states(1,4,:) = 0.0150;
channel(2).states(1,6,:) = 0.0150;
        
kTconstant = 1/25;
totalCurrent = zeros(protocol.nPoints,protocol.nSim); %creates a matrix for the current information
voltage = zeros(protocol.nPoints,protocol.nSim); %creates a matrix for the voltage information
time = zeros(protocol.nPoints,1);
for i = 2:protocol.nPoints
   time(i,1) = time((i-1),1) + protocol.dT; 
end

for n=1:protocol.nSim
    voltage(1,n) = parameters.iniPotential;
    startAt = 0;
    for ss=1:protocol.nSteps
        if(ss>1)
            iiiIni = 1;
        else
            iiiIni = 2;
        end
        for iii=iiiIni:protocol.nPointsPerStep(ss)
            for mmm=1:nModels
                dCurrent = 0;
                model(mmm).currentRates = updateRates03(model(mmm),protocol.ligand,voltage(iii+startAt-1,n),kTconstant);
                %currentRates = model(mmm).currentRates
                for x=1:model(mmm).nStates
                    channel(mmm).dStates(1,x)=0.0;
                    for y=1:model(mmm).nStates
                        channel(mmm).dStates(1,x)=channel(mmm).dStates(1,x)+model(mmm).currentRates(x,y)*channel(mmm).states(iii+startAt-1,y,n); %calculating dStates
                    end %for y=1:model(mmm).nStates
                end %for x=1:model(mmm).nStates
                %channel(mmm).dStates
                for y=1:model(mmm).nStates
                    channel(mmm).states(iii+startAt,y,n) = channel(mmm).states(iii+startAt-1,y,n)+channel(mmm).dStates(1,y); %updates states
                    dCurrent = dCurrent + model(mmm).stateConductance(1,y)*channel(mmm).states(iii+startAt,y,n)*(voltage(iii+startAt-1)-channel(mmm).reversalPotential);
                end
                channel(mmm).current(iii+startAt,n) = channel(mmm).current(iii+startAt,n) + dCurrent;
                totalCurrent(iii+startAt,n)=totalCurrent(iii+startAt,n) + channel(mmm).current(iii+startAt,n);
            end %for mmmm=1:nModels
            Ileak = voltage(iii+startAt-1,n)*parameters.Gleak + (voltage(iii+startAt-1,n)-reversal(2))*parameters.GKir;
            totalCurrent(iii+startAt,n)=totalCurrent(iii+startAt,n) + Ileak + protocol.externalCurrent(ss) + protocol.externalCurrent(ss)*(n-1);
            voltage(iii+startAt,n) = voltage(iii+startAt-1,n) - parameters.dtCap*totalCurrent(iii+startAt,n);
        end %for iii=2:protocol.nPointsPerStep(ss)
        startAt = startAt+protocol.nPointsPerStep(ss);
    end %for ss=1:protocol.nSteps
end %n=1:protocol.nSim

end
