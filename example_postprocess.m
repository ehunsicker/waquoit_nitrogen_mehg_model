addpath(genpath('modules'))
addpath(genpath('processing'))
OutputFile = 'OutputStochastic.mat';
load(OutputFile,'OObj')
[numSpc, numCo, numt, blah] = size(OObj.BigY);
[numszcohorts, blah, blah] = size(OObj.szY);
cohortsperyear = OObj.cohortsperyear;
szcohortsperyear = OObj.szcohortsperyear;
start_year = OObj.t0/365.25;
end_year = OObj.tend/365.25; 
cohortyears = numCo/cohortsperyear;
Cohorts = zeros(numSpc,cohortyears,numt,2);
Ymean = OObj.Ymean;
Ymin = OObj.Ymin;
Ymax = OObj.Ymax;

% Get info for tuna of a given age:
years = start_year:end_year;
ageout = zeros(numel(years),2);

spc = 12; % bluefin tuna
agemin = 13.5;
agemax = 14.4;
for t=1:numel(years)
     %get average tuna weight and MeHg for this age 
     % ageout(:,1) = mass, ageout(:,2) = MeHg
     ageout(t,:) = get_agewindow(spc,agemin,agemax,years(t),OObj);
end

% Get info over lifespans of given cohorts:
for c=1:cohortyears-1
    for spc=1:12
        Cohorts(spc,c,:,:) = get_cohort(spc,years(c),OObj);
    end
end

SpeciesNames = OObj.SpeciesNames;
timespan = OObj.get_timespan();