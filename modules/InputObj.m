classdef InputObj
    %INPUTOBJ 
    %   Object for handling input options and parameters
    
    properties
        ModelParams
        WaterParams
        temperature_source
        chla_source
        DOC_source
        MeHg_source
        tstart
        tend
        SpeciesParams
        SpeciesNames
        FeedPreferences
        FeedMins
        FeedMaxs
        ChlaTable
        ChlaSeasons
        isTseries
        isChlaseries
        isDOCseries
        isMeHgseries
        t_T
        t_Chla
        t_DOC
        t_MeHg
        Tinput
        DOCinput
        Chlainput
        MeHginput
   
    end
    
    methods
        function obj = InputObj(InputFile,temperature_source,...
                                chla_source, DOC_source, MeHg_source,...
                                tstart,tend)
                          
            obj.ModelParams = xlsread(InputFile,'Parameters','B4:B7'); %Model Parameters
            obj.WaterParams = xlsread(InputFile,'Parameters','B10:B24'); %Water Specific Parameters
            [obj.SpeciesParams, obj.SpeciesNames] = xlsread(InputFile,'Parameters','B29:AE40'); %Species 
            obj.SpeciesNames = obj.SpeciesNames(:,~all(cellfun('isempty',obj.SpeciesNames))); %Remove empty text columns
            obj.FeedPreferences = xlsread(InputFile,'FeedPreferences','B3:M14'); %Feeding Preferences
            obj.FeedMins = xlsread(InputFile,'FeedPreferencesMin','B3:M14'); % Feeding minima for stochastic
            obj.FeedMaxs = xlsread(InputFile,'FeedPreferencesMax','B3:M14'); % Feeding maxima room for stochastic
            obj.ChlaTable = xlsread(InputFile,'Chla','B2:M4'); % Chla Table 
            obj.ChlaSeasons = xlsread(InputFile,'Chla','A2:A4'); % and associated 'seasons'
            obj.temperature_source = temperature_source;
            obj.chla_source = chla_source;
            obj.DOC_source = DOC_source;
            obj.MeHg_source = MeHg_source;
            obj.tstart = tstart;
            obj.tend = tend;
            obj = obj.init_external_params();
        end
            
        function obj = init_temperature(obj)
            switch obj.temperature_source
                case 'default'
                    obj.isTseries = 'No';
                    obj.t_T = 'None';
                    obj.Tinput = 'None';
                otherwise
                    temperature_file = obj.temperature_source; 
                    [obj.isTseries, obj.Tinput, obj.t_T] = get_timeseries(temperature_file);
            end
        end
        
        function obj = init_chla(obj)
            switch obj.chla_source
                case 'default'
                    obj.isChlaseries = 'No';
                    obj.t_Chla = 'None';
                    obj.Chlainput = 'None';
                otherwise
                    Chla_file = obj.chla_source;
                    [obj.isChlaseries, obj.Chlainput, obj.t_Chla] = get_timeseries(Chla_file);
            end
        end

        function obj = init_DOC(obj)
            switch obj.DOC_source
                case 'default'
                    obj.isDOCseries = 'No';
                    obj.t_DOC = 'None';
                    obj.DOCinput = 'None';
                otherwise
                    DOC_file = obj.DOC_source;
                    [obj.isDOCseries, obj.DOCinput, obj.t_DOC] = get_timeseries(DOC_file);
            end
            
        end
        
        function obj = init_hg(obj)
             switch obj.MeHg_source
                case 'default'
                    obj.isMeHgseries = 'No';
                    obj.t_MeHg = 'None';
                    obj.MeHginput = 'None';
                otherwise
                    MeHg_file = obj.MeHg_source;
                    [obj.isMeHgseries, obj.MeHginput, obj.t_MeHg] = get_timeseries(MeHg_file);
            end   
        end
        
        function obj = init_external_params(obj)
            obj = obj.init_temperature();
            obj = obj.init_DOC();
            obj = obj.init_chla();
            obj = obj.init_hg();
        end
        
        function warming = get_temperature(obj,t)
            if strcmpi(obj.isTseries, 'No')
                warming = obj.ModelParams(2);      % Addition to base water temperature
            else
                warming = get_X_from_timeseries(t,obj.Tinput,obj.t_T);
            end

        end
        
        function chla = get_chla(obj,t)
            DOY = rem(t,365); % day of year for figuring out the month... no leap years so this will drift
            if strcmpi(obj.isChlaseries, 'No')
                season  = obj.ModelParams(1);      % 0 = Base Chla, 1 = High Chla, 2 = Low Chla 
                monthlyChla = obj.ChlaTable(season+1,:);
                month = floor(DOY/30.45)+1; % months as in 1/12th of the year, not calendar months
                chla = monthlyChla(month);
            else
                chla = get_X_from_timeseries(t,obj.Chlainput,obj.t_Chla);
            end

        end
        
        function DOC = get_DOC(obj,t)
            if strcmpi(obj.isDOCseries,'No')
                DOC = obj.WaterParams(5); % Dissolved Organic Carbon
            else
                DOC = get_X_from_timeseries(t,obj.DOCinput,obj.t_DOC);
            end
        end
        
        function MeHgWaterConc = get_MeHg(obj,t)
            if strcmpi(obj.isMeHgseries,'No')
                MeHgWaterConc = HgMethylation(obj.WaterParams(1),obj.WaterParams(6)); % Concentration of MeHg in water
            else
                MeHgWaterConc = get_X_from_timeseries(t,obj.MeHginput,obj.t_MeHg);
            end
        end
        
        function numSpecies = get_numSpecies(obj) 
            numSpecies = size(obj.SpeciesParams,1); % Number of species being modeled
        end
        
        function M0 = get_M0(obj)
            M0 = obj.SpeciesParams(:,22); %Initial Species Mass concentrations
        end
        
        function MeHg0 = get_MeHg0(obj)
            MeHg0 = obj.SpeciesParams(:,23); %Initial Species Hg concentration
        end
        
        function FinalAge = get_FinalAge(obj)
            FinalAge = obj.SpeciesParams(:,24); %Final species age for simulation (years)
        end
    end
end
