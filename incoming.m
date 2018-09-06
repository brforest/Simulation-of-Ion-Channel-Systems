function [incomingData] = incoming(filename)
incomingData = dlmread(filename, '\t', 1,1);
end
