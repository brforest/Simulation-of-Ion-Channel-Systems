function [model,conductance,nStates,nParams] = incomingModel(filename)
incomingData1 = dlmread(filename, '\t', 1,1);
nParams = incomingData1(1,1);
nStates = incomingData1(2,1);

incomingData2 = dlmread(filename, '\t', 3,1);
if nStates>5
    model = zeros(nParams,nStates);
else
    model = zeros(nParams,6);
end

for ii=1:nParams
    model(ii,:)=incomingData2(ii,:);
end

incomingData3 = dlmread(filename,'\t',3+nParams,1);
conductance(1,:)=incomingData3(1,:);
end
