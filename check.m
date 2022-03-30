% check:
%		is the function that verifies the termination condition for the 
%		compact genetic algorithm global optimization method.
%
% Created by: Nathaniel Despres
%
% Inputs:	
%		p -
%			is the probability vector for the cGA
%
% Outputs: 
%		terminate - 
%			is the termination logic value
%

function terminate = check(p)

terminate = 1;

for i = 1:length(p)
	if p(i) > 0 && p(i) < 1
		terminate = 0;
		break
	end
end
end