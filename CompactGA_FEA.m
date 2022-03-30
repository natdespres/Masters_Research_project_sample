% CompactGA_FEA
%			Performs the compact genetic algorithm (cGA) global
%			optimization method.
%
% Created by: Nathaniel Despres
%
% Inputs:
%		maxMin -
% 			is the parameter that selects whether to maximize or
% 			minimize the evaluation function.
%
%		initialProbability -
%			is the initial probability value used to initialise all
%			elements of the probability vector.
%
%		popSize - 
%			is the size of the cGA population.
%
%		maxNumIteration -
%			is the maximum number of iterations (or generations).
%
%		selectionSize -
%			is the selection size value for the cGA.
%
% 		nodesBD - [x, y]
% 			is the array of node locations for the boundary nodes where
% 			the x-positions are in the first column and the y-positions
% 			are in the second column.
% 
% 		nodesInt - [x, y]
% 			is the array of node locations for the internal nodes where 
% 			the x-positions are in the first column and the y-positions
% 			are in the second column.
% 
% 		mechProperties - [youngsMod; crossArea; secondMomentArea]
% 			is the array that contains the Young's modulus, cross
% 			sectional area and the second moment of area of the lattice
% 			structure.
% 
% 		boundaryConditions - [Node_number, DOF, Specified_disp(m)]
% 			is the array that contains the boundary conditions.
% 
% 		loadConditions - [Node_number, DOF, Load(N)]
% 			is the array that contains the load conditions.
%
% 		geometricProperties - [strutDiameter; fullLatticeVolume]
% 			is the array that contains the strut diameter and the volume
% 			of the lattice structure that has all nodes turned on.
% 
% 		performanceProperties - 
% 			[exponentRD; stiffStressRatio; exponentEL]
% 			is the array that contains parameters used to modify the
% 			lattice structure performance equation.
%
% Outputs:
%		nodesStar - [x, y]
%			is the array of node locations for the best structure where the 
%			x-positions are in the first column and the y-positions are
%			in the second column.
%
%		elementsStar - [nodeNumA, nodeNumB]
%			is the array that defines the elements of the best structure by
%			specifying which nodes are interconnected.
%
%		forceStar - [f1x, f1y, m1, f2x, f2y, m2, ...]
%			is a vector containing the reaction forces at each node.
%
%		displacementStar - [U1x, U1y, theta1, U2x, U2y, theta2, ...]
%			is a vector containing the displacement of each node.
%
%		performanceMax -
%			is the performance value for the best individual.
%
%		extraFEAinfo -
%			[maxDisplacement, maxStress, stiffness, relativeDensity,...
%					longElementLength]
%			is the array that contains extra information about the
%			performance and characteristics of the best structure.
%
%		numFuncEval -
%			is the number of function evaluations performed until the
%			algorithm terminates
%

function [nodesStar, elementsStar, forceStar, displacementStar, ...
		performanceMax, extraFEAinfo, numFuncEval]...
				= CompactGA_FEA(maxMin, initialProbability, popSize,...
						maxNumIteration, selectionSize, nodesBD, ...
						nodesInt, mechProperties, boundaryConditions,...
						loadConditions, geometricProperties, ...
						performanceProperties)

probVector = initialProbability*ones(1, size(nodesInt,1));

terminate = 0;

i = 1;

while ~terminate && i < maxNumIteration + 1
	
	x = cell(1,selectionSize);
	f = zeros(1,selectionSize);
	
	for j = 1:selectionSize
			[x{j}, f(j)] = generateFEA(probVector, nodesBD, nodesInt, ...
				mechProperties, boundaryConditions, loadConditions, ...
				geometricProperties, performanceProperties);
	end
	
	for m = 1:(selectionSize-1)
		for n = (m + 1):selectionSize
			[winner,looser] = compete(x{m},x{n},maxMin*f(m),maxMin*f(n));
			probVector = update(probVector, winner, looser, popSize);
		end
	end
	
	terminate = check(probVector);
	
	i = i + 1;
end

solutionIteration = i - 1;
numFuncEval = solutionIteration * selectionSize;

[nodesStar, elementsStar] = decodingFEA(winner, nodesBD, nodesInt);

[forceStar, displacementStar, performanceMax, extraFEAinfo]... 
	= FEA_PostProcess(nodesStar, elementsStar, mechProperties, ...
			boundaryConditions, loadConditions, geometricProperties, ...
			performanceProperties);
end


