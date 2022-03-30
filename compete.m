% compete:
%			is the function that performs the competion part of the compact
%			genetic algorithm global optimization method.
%
% Created by: Nathaniel Despres
%
% Inputs:	
% 		a -
% 			is a vector containing the encoded bit values for the first 
% 			individual.
% 
% 		b -
% 			is a vector containing the encoded bit values for the second 
% 			individual.
% 
% 		fA -
% 			is the function value for the "A" individual
% 
% 		fB -
% 			is the function value for the "B" individual
%			
% Outputs: 
%		winner - 
%			is a vector of encoded bits for the competion winner
%
%		looser - 
%			is a vector of encoded bits for the competion looser
%


function [winner, looser] = compete(a, b, fA, fB)

if fA >= fB
	winner = a;
	looser = b;
else
	winner = b;
	looser = a;
end
end
