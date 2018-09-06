function [totalCurrent, voltage, channel, Ileak, time] = calcResultRK01(model, protocol, parameters)
%channel = zeros(length(model));
reversal = [50 -90];
nModels = length(model)

for i=1:nModels
    channel(i).current = zeros(protocol.nPoints,protocol.nSim); %creates a matrix for the current information
    channel(i).statesMatrix = zeros(protocol.nPoints,model(i).nStates,protocol.nSim); %creates a matrix for the states information
    channel(i).nParams = model(i).nParams; %calculates the number of rates by the number of rows in the model.txt file
    channel(i).states = zeros(protocol.nPoints,model(i).nStates,protocol.nSim); %creates a matrix that stores the population fraction of each state
    channel(i).dStates = zeros(1,model(i).nStates);
    channel(i).k1 = zeros(1,model(i).nStates);
    channel(i).k2 = zeros(1,model(i).nStates);
    channel(i).k3 = zeros(1,model(i).nStates);
    channel(i).k4 = zeros(1,model(i).nStates);
    channel(i).reversalPotential = reversal(1,i);
end

ini1 = [0.777217177123013   0.165551476053305   0.011667292789928   0.000276239246875   0.000004629681723   0.000017543335374   0.001046420274508   0.044219221495274];
ini2 = [0.848547847789901   0.142254198488903   0.008943038417585   0.000249874639101   0.000002617976096   0.000002422688414];

for ii=1:protocol.nSim
    for jj=1:model(1).nStates
        channel(1).states(1,jj,ii) = ini1(1,jj);
    end
    for jj=1:model(2).nStates
        channel(2).states(1,jj,ii) = ini2(1,jj);
    end
end
        
kTconstant = 1/25;
totalCurrent = zeros(protocol.nPoints,protocol.nSim); %creates a matrix for the current information
voltage = zeros(protocol.nPoints,protocol.nSim); %creates a matrix for the voltage information
Ileak = zeros(protocol.nPoints,protocol.nSim);
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
                model(mmm).currentRates = updateRates03(model(mmm),protocol.ligand,voltage(iii+startAt-1,n),kTconstant);
                %currentRates = model(mmm).currentRates
                for x=1:model(mmm).nStates
                    channel(mmm).k1(1,x)=0.0;
                    for y=1:model(mmm).nStates
                        channel(mmm).k1(1,x)=channel(mmm).k1(1,x)+model(mmm).currentRates(x,y)*channel(mmm).states(iii+startAt-1,y,n); %calculating k1
                    end %for y=1:model(mmm).nStates
                end
                for x=1:model(mmm).nStates
                    channel(mmm).k2(1,x)=0.0;
                    for y=1:model(mmm).nStates
                        channel(mmm).k2(1,x)=channel(mmm).k2(1,x)+model(mmm).currentRates(x,y)*(channel(mmm).states(iii+startAt-1,y,n)+0.5.*channel(mmm).k1(1,y)); %calculating dStates
                    end %for y=1:model(mmm).nStates
                end
                for x=1:model(mmm).nStates
                    channel(mmm).k3(1,x)=0.0;
                    for y=1:model(mmm).nStates
                        channel(mmm).k3(1,x)=channel(mmm).k3(1,x)+model(mmm).currentRates(x,y)*(channel(mmm).states(iii+startAt-1,y,n)+0.5.*channel(mmm).k2(1,y)); %calculating dStates
                    end %for y=1:model(mmm).nStates
                end
                for x=1:model(mmm).nStates
                    channel(mmm).k4(1,x)=0.0;
                    for y=1:model(mmm).nStates
                        channel(mmm).k4(1,x)=channel(mmm).k4(1,x)+model(mmm).currentRates(x,y)*(channel(mmm).states(iii+startAt-1,y,n)+channel(mmm).k3(1,y)); %calculating dStates
                    end %for y=1:model(mmm).nStates
                end
                channel(mmm).dStates=(channel(mmm).k1+2.*channel(mmm).k2+2.*channel(mmm).k3+channel(mmm).k4)./6; %calculating dStates
                totalStates = 0.0;
                for y=1:model(mmm).nStates
                    channel(mmm).states(iii+startAt,y,n) = channel(mmm).states(iii+startAt-1,y,n)+channel(mmm).dStates(1,y); %updates states
                    if channel(mmm).states(iii+startAt,y,n)<0.0
                        channel(mmm).states(iii+startAt,y,n)=0.0;
                    end
                    totalStates = totalStates + channel(mmm).states(iii+startAt,y,n);                    
                end
                channel(mmm).states(iii+startAt,:,n) = channel(mmm).states(iii+startAt,:,n)/totalStates;                
                for y=1:model(mmm).nStates
                    channel(mmm).current(iii+startAt,n) = channel(mmm).current(iii+startAt,n)+ model(mmm).stateConductance(1,y)*channel(mmm).states(iii+startAt,y,n)*(voltage(iii+startAt-1,n)-channel(mmm).reversalPotential);
                end
                totalCurrent(iii+startAt,n)=totalCurrent(iii+startAt,n) + channel(mmm).current(iii+startAt,n);
            end %for mmmm=1:nModels
            Ileak(iii+startAt-1,n) = voltage(iii+startAt-1,n)*parameters.Gleak + (voltage(iii+startAt-1,n)-reversal(2))*parameters.GKir;
            totalCurrent(iii+startAt,n)=totalCurrent(iii+startAt,n) + Ileak(iii+startAt-1,n) + protocol.externalCurrent(ss) + protocol.dExternalCurrent(ss)*(n-1);
            voltage(iii+startAt,n) = voltage(iii+startAt-1,n) - parameters.dtCap*totalCurrent(iii+startAt,n);
        end %for iii=2:protocol.nPointsPerStep(ss)
        startAt = startAt+protocol.nPointsPerStep(ss);
    end %for ss=1:protocol.nSteps
end %n=1:protocol.nSim

channel(1).states(protocol.nPoints,:,1)
channel(2).states(protocol.nPoints,:,1)
end
