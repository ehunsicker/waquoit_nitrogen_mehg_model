
addpath(genpath('modules'))
InputFile = 'InputParameters.xlsx';
OutputFile = 'OutputStochastic.mat';
nIterations = 1;
start_year = 1;
end_year = 25;

temperature_source = 'GoM temp example.xlsx'; % default values defined in InputFile
chla_source = 'default';
MeHg_source = 'default';
DOC_source = 'default';

ProbabilisticModelStochastic(InputFile, OutputFile,...
                            nIterations,start_year,end_year,...
			    temperature_source, chla_source, ...
				    MeHg_source, DOC_source)
example_plotscript