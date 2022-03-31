function ProbabilisticModelStochastic(InputFile, OutputFile, ...
                                    nIterations,startyear,endyear, ...
				    temperature_source, chla_source, ...
				    MeHg_source, DOC_source)

%This script runs a stochastic version of the model by selecting parameters
%from predefined distributions. Inputs can be changed in the excel file
%listed in the variable InputFile below.

%Load input parameters from external excel worksheets

%Model Parameters

IObj = InputObj(InputFile,temperature_source,...
                  chla_source, DOC_source, MeHg_source,...
                    startyear,endyear);

numSpecies = IObj.get_numSpecies(); % Number of species being modeled
M0 = IObj.get_M0(); %Initial Species Mass concentrations
MeHg0 = IObj.get_MeHg0; %Initial Species Hg concentration
FinalAge = IObj.get_FinalAge(); %Final species age for simulation (years)

Y = cell(numSpecies,1);
T = cell(numSpecies,1);

%% Progress bar to be displayed
tic
h = waitbar(0,'Please wait...');
      
masterdt = 365.25/48;
tstart = startyear*365.25;
tend = (endyear+1)*365.25; %finish the last year
numYears = endyear-startyear;
timespan = tstart:masterdt:tend;
cohortsperyear = 1;
BigY = NaN(numSpecies,numYears*cohortsperyear,numel(timespan),2);
Yi = 0.;
n=4;
szpreps = 12;
szY = NaN((numYears-1)*szpreps,numel(timespan),2);
Ymin = cell(numSpecies,numYears*cohortsperyear);
Ymax = cell(numSpecies,numYears*cohortsperyear);
Ymean = cell(numSpecies,numYears*cohortsperyear);
Ystd = cell(numSpecies,numYears*cohortsperyear);

% Small Zooplankton
for cohort = 0:(numYears-1)*szpreps
    dt = FinalAge(n)*365.25/200;
    thistend = min(tstart+(cohort/szpreps+FinalAge(n))*365.25,tend);
    thistimespan = (tstart+cohort*365.25/szpreps):dt:thistend;
    seed_base = 1+5000;  
    y0 = [M0(n) MeHg0(n)]; %Initial Conditions
       %Solve ODE for SZP
    ode = @(t,y) odeBioaccumulationStochastic(t, y, n, BigY, szY,...
            IObj,...
            tstart, masterdt, seed_base);
    [t,y] = ode_solve(ode,thistimespan,y0);
       %First row of y variable is mass of species and second row is
       %Concentration of mercury                
    if strcmp(IObj.SpeciesNames{n,1}, 'Small Zooplankton')
        y(:,1) = y(:,1)* 2 * 5 * 1e-6; % ug-C * ug-dry wt/ug-C * ug-wet wt/ug-dry wt * (1e-6) = g;
    end
    
    [m,ti1] = min(abs(timespan-thistimespan(1)));
    [m,ti2] = min(abs(timespan-thistimespan(end)));
    szY(cohort+1,ti1:ti2,1) = interp1(thistimespan,y(:,1),timespan(ti1:ti2));
    szY(cohort+1,ti1:ti2,2) = interp1(thistimespan,y(:,2),timespan(ti1:ti2));
    end

% Whole food web
for i = 1:nIterations 
    IObj.init_external_params()       
    %Loop through each species
    for n = 5:numSpecies
        for cohort = 0:(numYears-1)*cohortsperyear
            dt = FinalAge(n)*365.25/300;
            thistend = min(tstart+(cohort/cohortsperyear+FinalAge(n))*365.25,tend);
            thistimespan = tstart+cohort*365.25/cohortsperyear:dt:thistend;
            Yi=0;        
            seed_base = i*n+1+6000;  
            y0 = [M0(n) MeHg0(n)]; %Initial Conditions
               %Solve ODE, one species at a time
            ode = @(t,y) odeBioaccumulationStochastic(t, y, n, BigY, szY,...
                    IObj,...
                    tstart, masterdt, seed_base);
            [t,y] = ode_solve(ode,thistimespan,y0);
               %First row of y variable is mass of species and second row is
               %Concentration of mercury                
            if strcmp(IObj.SpeciesNames{n,1}, 'Small Zooplankton')
                y(:,1) = y(:,1)* 2 * 5 * 1e-6; % ug-C * ug-dry wt/ug-C * ug-wet wt/ug-dry wt * (1e-6) = g;
            end
            if strcmp(IObj.SpeciesNames{n,1}, 'Heterotrophic Protist')
                   y(:,1) = y(:,1)/N_cells; % ug-C * ug-dry wt/ug-C * ug-wet wt/ug-dry wt * (1e-6) = g;
            end

            Yi = y;                  

            [m,ti1] = min(abs(timespan-thistimespan(1)));
            [m,ti2] = min(abs(timespan-thistimespan(end)));
 
            BigY(n,cohort+1,ti1:ti2,1) = interp1(thistimespan,y(:,1),timespan(ti1:ti2));
            BigY(n,cohort+1,ti1:ti2,2) = interp1(thistimespan,y(:,2),timespan(ti1:ti2));

            if i == 1
               Ymin{n,cohort+1}(:,1) = Yi(:,1);
               Ymax{n,cohort+1}(:,1) = Yi(:,1);
               Ymean{n,cohort+1}(:,1) = Yi(:,1);
               Ystd{n,cohort+1}(:,1) = Yi(:,1)*0.;
               Ymin{n,cohort+1}(:,2) = Yi(:,2);
               Ymax{n,cohort+1}(:,2) = Yi(:,2);
               Ymean{n,cohort+1}(:,2) = Yi(:,2);
               Ystd{n,cohort+1}(:,2) = Yi(:,2)*0.;
            elseif i == 2
               Ymin{n,cohort+1}(:,1) = min(Ymin{n,cohort+1}(:,1),Yi(:,1));
               Ymax{n,cohort+1}(:,1) = max(Ymax{n,cohort+1}(:,1),Yi(:,1));
               Ydel = Yi(:,1)-Ymean{n,cohort+1}(:,1);
               Ymean{n,cohort+1}(:,1) = Ydel/i+Ymean{n,cohort+1}(:,1);
               Ystd{n,cohort+1}(:,1) = Yi(:,1) - Ymean{n,cohort+1}(:,1);
               Ymin{n,cohort+1}(:,2) = min(Ymin{n,cohort+1}(:,2),Yi(:,2));
               Ymax{n,cohort+1}(:,2) = max(Ymax{n,cohort+1}(:,2),Yi(:,2));
               Ydel = Yi(:,2)-Ymean{n,cohort+1}(:,2);
               Ymean{n,cohort+1}(:,2) = Ydel/i+Ymean{n,cohort+1}(:,2);
               Ystd{n,cohort+1}(:,2) = Yi(:,2) - Ymean{n,cohort+1}(:,2);
            else
	       Ymin{n,cohort+1}(:,1) = min(Ymin{n,cohort+1}(:,1),Yi(:,1));
               Ymax{n,cohort+1}(:,1) = max(Ymax{n,cohort+1}(:,1),Yi(:,1));
               [Ymean{n,cohort+1}(:,1),Yvar] = welford(Yi(:,1), Ymean{n,cohort+1}(:,1), Ystd{n,cohort+1}(:,1).^2,i);
               Ystd{n,cohort+1}(:,1) = Yvar.^0.5;
               Ymin{n,cohort+1}(:,2) = min(Ymin{n,cohort+1}(:,2),Yi(:,2));
               Ymax{n,cohort+1}(:,2) = max(Ymax{n,cohort+1}(:,2),Yi(:,2));
               [Ymean{n,cohort+1}(:,2),Yvar] = welford(Yi(:,2), Ymean{n,cohort+1}(:,2), Ystd{n,cohort+1}(:,2).^2,i);
               Ystd{n,cohort+1}(:,2) = Yvar.^0.5;
       	    end

	end
				           
        waitbar((n-3) / (numSpecies-3),h,IObj.SpeciesNames{min(n+1,numSpecies),1})

    end     
end

% Save outputs to container 
OObj = OutputObj;
OObj.BigY = BigY;
OObj.szY = szY;
OObj.szcohortsperyear = szpreps;
OObj.t0 = timespan(1);
OObj.tend = timespan(end);
OObj.dt = masterdt;
OObj.SpeciesNames = IObj.SpeciesNames(:,1);
OObj.cohortsperyear = cohortsperyear;
OObj.Ymean = Ymean
OObj.Ymin = Ymin
OObj.Ymax = Ymax
OObj.Ystd = Ystd
save(OutputFile,'OObj')
close(h)
toc

end
