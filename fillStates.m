function [state] = fillStates(nStates, parameter)
state(1,1) = 1.0;
for i = 2:nStates
%    state(1,i) = parameter(i);
    state(1,i) = 0.0;
end

end