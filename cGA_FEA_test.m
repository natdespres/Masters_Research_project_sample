% cGA_FEA_test:
%			Performs the compact genetic algorithm (cGA) global
%			optimization method many times in order to obtain averages and
%			standard deviation.
%
% Created by: Nathaniel Despres
%
% Inputs:
% 		data_foldername -
% 			is a string of the foldername for where the folder that
% 			contains the test data and the folder that contains the
% 			best lattice structure figures are found.
% 
% 		test_foldername -
% 			is a string of the foldername for where all the data of
% 			the tests are saved.
% 
% 		struct_foldername -
% 			is a string of the foldername for where all the figures of
% 			the best lattice structure are saved.
% 
% 		numRuns2Average -
% 			is the number of runs that the GA will be performed using
% 			the same parameters.
% 
% 		maxMin -
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

function cGA_FEA_test(data_foldername, test_foldername, ...
					struct_foldername, numRuns2Average, maxMin, ...
					initialProbability, popSize, maxNumIteration, ...
					selectionSize, nodesBD, nodesInt, mechProperties, ...
					boundaryConditions, loadConditions, ...
					geometricProperties, performanceProperties)

rng('shuffle')

nodesStarData = {};
elementsStarData = {};
forceStarData = {};
displacementStarData = {};
performanceMaxData = zeros(1,numRuns2Average);
maxDisplacementData = zeros(1,numRuns2Average);
maxStressData = zeros(1,numRuns2Average);
stiffnessData = zeros(1,numRuns2Average);
relativeDensityData = zeros(1,numRuns2Average);
longElementLengthData = zeros(1,numRuns2Average);
numFuncEvalData = zeros(1,numRuns2Average);

parfor (i = 1:numRuns2Average,10)
	
	[nodesStarData{1,i}, elementsStarData{1,i}, forceStarData{1,i}, ...
	displacementStarData{1,i}, performanceMaxData(i), extraFEAinfo, ...
	numFuncEvalData(i)]...
		= CompactGA_FEA(maxMin, initialProbability, popSize, ...
				maxNumIteration, selectionSize, nodesBD, nodesInt, ...
				mechProperties, boundaryConditions, loadConditions, ...
				geometricProperties, performanceProperties);
	
	maxDisplacementData(i) = extraFEAinfo(1);
	maxStressData(i) = extraFEAinfo(2);
	stiffnessData(i) = extraFEAinfo(3);
	relativeDensityData(i) = extraFEAinfo(4);
	longElementLengthData(i) = extraFEAinfo(5);
end

averagePerformanceMax = mean(performanceMaxData);
stdPerformanceMax = std(performanceMaxData);

averageMaxDisplacement = mean(maxDisplacementData);
stdMaxDisplacement = std(maxDisplacementData);

averageMaxStress = mean(abs(maxStressData));
stdMaxStress = std(abs(maxStressData));

averageStiffness = mean(abs(stiffnessData));
stdStiffness = std(abs(stiffnessData));

averageRelativeDensity = mean(abs(relativeDensityData));
stdRelativeDensity = std(abs(relativeDensityData));

averageLongElementLength = mean(longElementLengthData);
stdLongElementLength = std(longElementLengthData);

averageNumFuncEval = mean(numFuncEvalData);
stdNumFuncEval = std(numFuncEvalData);

[best,bestIndex] = max(performanceMaxData);

StructurePlotting(nodesStarData{1,bestIndex}, ...
					elementsStarData{1,bestIndex}, ...
					displacementStarData{1,bestIndex})

title(['cGA graph ss',num2str(selectionSize), ...
				' ip',num2str(initialProbability), ...
				' ps',num2str(popSize), ...
				' eRD',num2str(performanceProperties(1)), ...
				' eEL',num2str(performanceProperties(3)), ...
				' LP',num2str(best,'%.2e'), ...
				' fe',num2str(numFuncEvalData(bestIndex),'%i')])

struct_fig_filename = ['cGAstruct_ss',num2str(selectionSize), ...
					'_ip',strrep(num2str(initialProbability),'.',''), ...
					'_ps',num2str(popSize), ...
					'_eRD',num2str(performanceProperties(1)), ...
					'_eEL',num2str(performanceProperties(3)), ...
					'_LP',strrep(num2str(best/100,'%.2e'),'.',''), ...
					'_fe',num2str(int64(averageNumFuncEval),'%i')];

saveas(gcf,[data_foldername,'/',test_foldername,'/',...
			struct_foldername,'/',struct_fig_filename])
		
saveas(gcf,[data_foldername,'/',test_foldername,'/',...
			struct_foldername,'/',struct_fig_filename,'.pdf'])

test_filename = strcat('cGAdata_ss',num2str(selectionSize), ...
					'_ip',strrep(num2str(initialProbability),'.',''), ...
					'_ps',num2str(popSize), ...
					'_eRD',num2str(performanceProperties(1)), ...
					'_eEL',num2str(performanceProperties(3)), ...
					'_fe',num2str(int64(averageNumFuncEval),'%i'));

save(fullfile(data_foldername, test_foldername, test_filename))

