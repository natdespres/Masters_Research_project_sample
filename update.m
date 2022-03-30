% update:
%			is the function that updates the probability vector for the 
%			compact genetic algorithm global optimization method.
%
% Created by: Nathaniel Despres
% 
% Inputs:
% 		p -
% 			is the probability vector for the cGA
% 
% 		winner -
% 			is a vector containing the encoded bit values for the 
% 			winning individual.
% 
% 		looser -
% 			is a vector containing the encoded bit values for the  
% 			looser individual.
% 
% 		popSize - 
% 			is the size of the cGA population
%			
% Outputs: 
%		p -
%			is the updated probability vector for the cGA
%

function p = update(p, winner, looser, popSize)

for i = 1:length(p)
	if winner(i) ~= looser(i)
		if winner(i) == 1
			p(i) = p(i) + (1/popSize);
		else
			p(i) = p(i) - (1/popSize);
		end
	end
end

for i = 1:length(p)
	if p(i) < 0.5 && p(i) < (1/(1000*popSize))
		p(i) = 0;
	end
	if p(i) > 0.5 && (1 - p(i)) < (1/(1000*popSize))
		p(i) = 1;
	end
end

end