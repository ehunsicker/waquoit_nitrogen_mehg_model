function dy = odeBioaccumulationStochastic(t,y,...
    n,BigY, szY, ...
    IObj, ...
    tstart, masterdt, seed)
% ODEBIOACCUMULATIONSTOCHASTIC
%   Set up the equation for growth and MeHg accumulation
%   for each species.

global lastlzpprey
SpeciesParams = IObj.SpeciesParams;
SpeciesNames = IObj.SpeciesNames;
ChlaTable = IObj.ChlaTable;
ChlaSeasons = IObj.ChlaSeasons;
FeedPreferences = IObj.FeedPreferences(n,:);
FeedMins = IObj.FeedMins(n,:);
FeedMaxs = IObj.FeedMaxs(n,:);


sizeBY = size(BigY);
ti=min(sizeBY(3),int16((t-tstart)/masterdt)+1);
rng(seed*ti)
%% Define variables from input arguments 

% MODEL PARAMETERS

warming = IObj.get_temperature(t);

chla = IObj.get_chla(t);

% WATER PARAMETERS
MeHgWaterConc = IObj.get_MeHg(t);
Dow = IObj.WaterParams(2); % Octanol-Water partitioning coeff.
sat = IObj.WaterParams(3); %
MeHg_sed = IObj.WaterParams(4); % Concentration of MeHg in Sediment

DOC = IObj.get_DOC(t);    
ScEff = IObj.WaterParams(7); % Scavenging efficiency
aE = IObj.WaterParams(13); 
bE = IObj.WaterParams(14); 
cE = IObj.WaterParams(15); 
neta = (1.87+155/Dow)^-1;
% FISH PARAMETERS

Species = IObj.SpeciesNames{n,1};
TemperatureType = IObj.SpeciesNames{n,2};

%Consumption Parameters
ac = IObj.SpeciesParams(n,1); 
bc = IObj.SpeciesParams(n,2); 
Qc = IObj.SpeciesParams(n,3); 
P = IObj.SpeciesParams(n,4); 

%Respiration Parameters
ar = IObj.SpeciesParams(n,5); 
br = IObj.SpeciesParams(n,6); 
Qr = IObj.SpeciesParams(n,7); 
S = IObj.SpeciesParams(n,8); 
A = IObj.SpeciesParams(n,9); %Base of juvenile activity level
MigWeight = IObj.SpeciesParams(n,25); 
dr = IObj.SpeciesParams(n,26); 
aa = IObj.SpeciesParams(n,27); 
ba = IObj.SpeciesParams(n,28); 

%Egestion
alpha_ege = IObj.SpeciesParams(n,10); 
beta_ege = IObj.SpeciesParams(n,11); 
gamma_ege = IObj.SpeciesParams(n,12); 

%Excretion
alpha_excr = IObj.SpeciesParams(n,13);
beta_excr = IObj.SpeciesParams(n,14);
gamma_excr = IObj.SpeciesParams(n,15);

%Other
TempPref = IObj.SpeciesParams(n,16);
EnerCont = IObj.SpeciesParams(:,18);
aL = IObj.SpeciesParams(:,19);
bL = IObj.SpeciesParams(:,20);
f_assim = IObj.SpeciesParams(n,21);
del_f_assim = IObj.SpeciesParams(n,29);

%Unit Conversions
MeHgWaterConc  = MeHgWaterConc * 1e-12 * 200.59 * 1e9;  % mol/L g/L ng/L

MeHg_sed       = MeHg_sed * 1e-12 * 200.59 * 1e9;  % ng/g

%Temperature Parameters
Toc = (TempPref + 0.53)./1.05;
Tmr = 0.66 * TempPref + 16.43;
Tmc = Tmr - 3; %Harvey et al (2009), Hansen et al (1997)
Tor = Tmc;
    
Chla =chla;

[MeHgMicroP MeHgNanoP MeHgPicoP] = getPhytoplanktonHg(DOC,MeHgWaterConc);

switch TemperatureType
    
    case 'Warm'
         Temp = 6 + 6*(1-cos(2*pi*(t-30)/365)) + warming;
         
    case 'Temperate'       
         Temp = 6 + 2*(1-cos(2*pi*(t-30)/365)) + warming;

    case 'Cold'
         Temp = 6 + 1*(1-cos(2*pi*(t-30)/365)) + warming;
    otherwise
        disp(TemperatureType)
end   
    
%%
%%%%%%%%%%%%%%%%%%%%
%   Bioenergetics  %
%%%%%%%%%%%%%%%%%%%%

% Bioegnergetics equations used to calculate growth (mass / time)
%Different calculations for growth based on species
switch Species
   
    case 'Heterotrophic Protist' % 
        
        g = 0.1604*log10(Chla)+0.3711; % /d growth rate according to Berge et al. 2008 

        dy(1,1) = g * y(1); %biomass g/d
        
        Cell_mass = 6e-10*log10(Chla)+2e-9; % mass of each cell in g
       
        N_cells = dy(1,1)/Cell_mass;    
    
    case 'Small Zooplankton' % equations are from Hirst and Bunker 2003
        
        log10_g = 0.0186*Temp - 0.288 * log10(y(1)) + 0.417 * log10(Chla) - 1.209;

        g = 10^(log10_g);

        dy(1,1) = g * y(1); %mass in ug C/d
           
    case 'Polychaeta' 
        
        Tmg = 32; 
        Tog = 15; 
        Qg = 2.6;
        
        Vc   = (Tmg - Temp)/(Tmg - Tog);
        Wc   = log(Qg) * (Tmg - Tog);
        Yc   = log(Qg) * (Tmg - Tog + 2);
        Xc   = ((Wc^2) * (1 + (1+(40/Yc))^0.5)^2)/400;
        fT   = (Vc^Xc)* exp(Xc * (1-Vc));
        
        dy(1,1) = fT * 0.0107 * log(Chla) * y(1); %from relationship in Vedel&Riisgard 1993
            
    otherwise % All other species
        
        %Consumption
        Vc   = (Tmc - Temp)/(Tmc - Toc);
        Wc   = log(Qc) * (Tmc - Toc);
        Yc   = log(Qc) * (Tmc - Toc + 2);
        Xc   = ((Wc^2) * (1 + (1+(40/Yc))^0.5)^2)/400;
        rc   = (Vc^Xc)* exp(Xc * (1-Vc));
        
        Cons  = ac*(y(1)^bc) * P * rc; %Rate of Consumption

        %Set consumption to zero if it's too cold or if negative
        if  Temp < 1 || Cons <0
                Cons = 0;
                A=1;
        end

        %Egestion
        Ege = Cons * alpha_ege * Temp^(beta_ege) * exp(gamma_ege * P);
        
        %Respiration
        Vr   = (Tmr - Temp)/(Tmr - Tor);
        Wr   = log(Qr) * (Tmr - Tor);
        Yr   = log(Qr) * (Tmr - Tor + 2);
        Xr   = ((Wr^2) * (1 + (1+(40/Yr))^0.5)^2)/400;
        rr   = (Vr^Xr)* exp(Xr * (1-Vr));
        
        
        if MigWeight>0
            if y(1)>MigWeight  % weight at first migration
      
            U=aa*(y(1)^ba); %swimming speed cm/s from Rudstam et al 1988
            A = A+exp(dr*U); % Here A is the active metabolic rate and exp () is scales routine metabolism to active metabolism Schindler2002

            end
        end
        
        Resp  = ar*(y(1)^br) * rr * A + S* (Cons-Ege) ; %Rate of Respiration
        % added -egestion, because S is the energy cost of digestion,
        % see Megrey 2007   

        %Set respiration to 0 if calculated as negative
        if Resp < 0
           Resp = 0;
        end


        %Excretion
        Excr = (Cons-Ege) * alpha_excr * Temp^(beta_excr) * exp(gamma_excr * P);
        %added -egestion because egested fish doesnt get excreted see Megrey2007
        %Mass balance (g/time)
        dy(1,1) = (Cons - Resp - Ege - Excr) * y(1);   
        
end

%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MEHG BIOACCUMULATION   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Small Zooplankton, Mussels, and Polychaeta (worms) are
% treated differently in their feeding preferences 
%
% Polychaeta - Nelson et al (1995)
%

if strcmp(Species, 'Small Zooplankton') || ...
        strcmp(Species, 'Heterotrophic Protist') || ...
        strcmp(Species, 'Quahog') || ...
        strcmp(Species, 'Polychaeta')
            
    % Calculation of depth-averaged suspended matter concentration
    log10_Zeu   =  1.524-0.46*log10(Chla)-0.00051*(log10(Chla))^2+0.0282*(log10(Chla))^3;   % Morel et al.(2007)
    Zeu         =  10^(log10_Zeu);                      % euphotic depth
    Chla_ave    =  58.5*(Chla^0.546)/(1.5*Zeu);         % depth_averaged Chla               % Uitz et al.(2006) mixed waters from Table 4
    theta       =  80;                                  % C:Chla ratio, ranges between 30-130(Gallegos et al.,2011)
    SSconc      =  2.6 * Chla_ave * theta/1000;         % g/m3; (Gallegos et al.,2011) - dry wt.
    SSconc      =  SSconc/1000;                         % g/L  - dry wt.
    SSconc      =  SSconc * 5;                          % g

    %Phytoplankton feeding preferences, Updated by AS from Uitz et al
    %(2006) Plotted Mixed water data from Table 6
    % Chla restrictions are from Chrisholm
    FeedPreferences(1) = 0.1413*(Chla^-0.642); % Pico
        if Chla < 0.5
             FeedPreferences(1) = 1;
        end
    FeedPreferences(2) = 0.3291*(Chla^-0.673); % Nano
        if Chla < 0.5
             FeedPreferences(2) = 0;
        end
    FeedPreferences(3) = 1 - FeedPreferences(1) - FeedPreferences(2); %Micro
         if FeedPreferences(3) < 0
             FeedPreferences(3) = 0;
         end
    if sum(FeedPreferences(1:3))>1
        FeedPreferences
    end
    
end  

%Different calculations for meHg accumulation based on species
switch Species
    
    case 'Heterotrophic Protist' % 14-17 um in diameter eats mosltly nanoplankton of same size or slightly larger
        FeedPreferences(1) = 0;
        FeedPreferences(2) = 0.9;
        FeedPreferences(3) = 0.1;
     
        ClearanceRate = 5e-8 * Chla^-0.438 * N_cells *24; %L/cell/h *cell*h/d=L/d
              
        % MEHG BIOACCUMULATION   

        MeHgWaterUptake = (ClearanceRate) * neta * MeHgWaterConc/ y(1);  % L/d * unitless * ng/L * (1/g-biomass)

        MeHgAssim =  (ClearanceRate) * SSconc * ScEff * (MeHgNanoP*FeedPreferences(2)...
            + MeHgPicoP*FeedPreferences(1) + MeHgMicroP*FeedPreferences(3)) * f_assim/y(1);  % ng-MeHg/g-biomass/d 
      

        dy(2,1) =  MeHgWaterUptake + MeHgAssim - y(2) * (1/y(1)) * dy(1,1); % ng-MeHg/g-biomass/d
        
        
        dy(2,1) =dy(2,1)/N_cells; % ng-MeHg/g-cell/d
        
        MeHgCiliate = y(2);
       
  

    case 'Small Zooplankton'
        
       DryWt    = y(1) * 2 * 1e-3; % ug-C * (dry. wt in ug/ mg-C) * (1e-3) = dry wt. in mg          
       
       b   = 1.777 * exp(0.234 * Temp);
       n2   = 0.681 * exp(0.0199 * Temp);
       
       ClearanceRate  = b * DryWt^n2 * 24;     % ml/h * h/d = mL/d
           
       WetWt = y(1) * 2 * 5 * 1e-6; % ug-C * ug-dry wt/ug-C * ug-wet wt/ug-dry wt * (1e-6) = g
       
        % MEHG BIOACCUMULATION   
       f_nano=1;
        
       MeHgWaterUptake = (ClearanceRate/1000) * neta * MeHgWaterConc/ WetWt;  % L/d * unitless * ng/L * (1/g-ZS)

       MeHgAssim =  (ClearanceRate / 1000) * SSconc * ScEff * ((f_nano*MeHgNanoP*FeedPreferences(2) ...
            + MeHgPicoP*FeedPreferences(1) + MeHgMicroP*FeedPreferences(3)) * f_assim/(WetWt));  % ng-MeHg/g-ZS/d 


       MeHgElim = aE * (((WetWt))^(bE)) * exp(cE*Temp) * y(2);
      

       dy(2,1) =  MeHgWaterUptake + MeHgAssim - MeHgElim - y(2) * (1/y(1)) * dy(1,1); % ng-MeHg/g-ZS/d 
        
    case 'Polychaeta'
        FeedPreferences(1) = 0;
        FeedPreferences_sed = 0.2;
        FeedPreferences(2) = (FeedPreferences(2)/(FeedPreferences(2)+FeedPreferences(3)))*(1-FeedPreferences_sed); % nano  in the fraction phytoplankton diet
        FeedPreferences(3) = (FeedPreferences(3)/(FeedPreferences(2)+FeedPreferences(3)))*(1-FeedPreferences_sed); % micro in the fraction phytoplankton diet
        f_nano=1;
        
        F = 34.65 * y(1)^1.1292;  % water filtration rate in L/d - Vedel and Riisgard 1993.
       
        MeHgAssim  =  (F * SSconc * ScEff * (f_assim * (MeHgMicroP*FeedPreferences(3)+f_nano*MeHgNanoP*FeedPreferences(2)) + f_assim*MeHg_sed*FeedPreferences_sed ) / y(1));
        neta = 0.08;
        MeHgWaterUptake = F * neta *  MeHgWaterConc / y(1);  %ng-MeHg/g-fish/d

        MeHgElim = aE * (y(1)^(bE)) * exp(cE*Temp) * y(2); 

        dy(2,1) = MeHgWaterUptake + MeHgAssim - MeHgElim - y(2) * (1/y(1)) * dy(1,1);  % ng-MeHg/g-fish/d 

    case 'Quahog'
        FeedPreferences(1) = 0;
        FeedPreferences_sed = 0.2; 
        FeedPreferences(2) = (FeedPreferences(2)/(FeedPreferences(2)+FeedPreferences(3)))*(1-FeedPreferences_sed); 
        FeedPreferences(3) = (FeedPreferences(3)/(FeedPreferences(2)+FeedPreferences(3)))*(1-FeedPreferences_sed);
        
        f_nano=1;
               
        W_dry = y(1)/5;
        W_m = W_dry * 5;
        F=3*W_dry^0.84;%L/h from Riisgard 2001
        F=F*24;%L/d

        neta = 0.44;
            
        MeHgWaterUptake = F * neta * MeHgWaterConc / W_m;  %ng-MeHg/g-fish/d
        
        % Dietary uptake
        SSconc_Sed = SSconc;
        MeHgAssim = ((F * SSconc * ScEff *f_assim*(FeedPreferences(1) * MeHgPicoP + FeedPreferences(2) * MeHgNanoP*f_nano + FeedPreferences(3)*MeHgMicroP))...
                    + F * SSconc_Sed * ScEff * (FeedPreferences_sed*MeHg_sed) * f_assim/W_m);  
        
        MeHgElim = aE * ((W_m)^(bE)) * exp(cE*Temp) * y(2);

        dy(2,1) = MeHgWaterUptake + MeHgAssim - MeHgElim - y(2) * (1/W_m) * dy(1,1);
    
    
    otherwise %All other species
        
        ReqEnergy = Cons * EnerCont(n); % Cons: g prey/g pred *time EnerCont: KJ/g pred ReqEnergy: KJ/g prey 
      
       
        Iprey = find(FeedPreferences); %Index of prey organisms
        [FeedPreferences] = StochasticFeedPreferences(FeedMaxs,FeedMins,FeedPreferences,seed);
        
        
        %Construct matrix of feed preference / Energy constraints
        %
        AA = eye(length(Iprey)) - repmat(FeedPreferences(Iprey),length(Iprey),1)';
        AA = [AA; EnerCont(Iprey)'];

        B = [zeros(1,size(AA,2)) ReqEnergy]';
       
        X = AA\B;
        MeHgPrey = zeros(1,length(Iprey)); %Initialize MeHg in each prey
                
        %Estimate the length (size) of organism based on mass
        if strcmp(Species, 'Large Zooplankton')
            Length_y = ((y(1)*1e3)/aL(n))^(1/bL(n));   % in mm; reation is W(mg) = aL(mm)^b.  Krill(1) is in grams, so 1e3 is required.
            MeHgPrey(1:3) = [MeHgPicoP MeHgNanoP MeHgMicroP];
        else
            Length_y = (y(1)/aL(n))^(1/bL(n));  %Length (Size) of organism in mm
        end
        
        %Determining preferred size/weight of each prey
        
        Size_ratio=0.5593*log(Length_y)+4.2006; % calculates size ratios based on pred length 
        if strcmp(Species, 'Swordfish')||... %ratio found in Young et al 2006 (ranging from 2-5)
            strcmp(Species, 'Lobster')||...
            strcmp(Species, 'Rock Crab')
            Size_ratio=4;
        end
        if strcmp(Species, 'Bluefish')||...
            strcmp(Species, 'Striped Bass')
            Size_ratio=2.0635*log(Length_y)-7.3435; % Scharf et al 2000
        end
        if strcmp(Species, 'Squid')
            Size_ratio=1;
        end
        
       %Loop backwards from calculated prey sizes to find largest to
       %eat determined by size_ratio
       if strcmp(Species, 'Large Zooplankton')

            count=1;
            massPreyArray = squeeze(szY(:,ti,1));
            MeHgPreyArray = squeeze(szY(:,ti,2));
            aLPrey = aL(4);
            bLPrey = bL(4);
	    Size_ratio=10;
            for i = length(massPreyArray):-1:2

                Length_Prey = ((massPreyArray(i)*1e3)/aLPrey)^(1/bLPrey); 
                if massPreyArray(i) > (4.9e-5/1e3)  % this particular line ensures that the copepod is at least > 0.125 mm (copepod nauplii). otherwise, krill eats microP
                    if Length_Prey < Length_y/Size_ratio % 
                        MeHgPrey(4) = MeHgPreyArray(i);
                        count=2;
                        lastlzpprey = MeHgPreyArray(i);
                        break
                    end                        
                end
            end
                
            if count == 1
                X(4) = X(4) * EnerCont(4)/EnerCont(2);
                if numel(lastlzpprey) > 0
                MeHgPrey(4) = lastlzpprey;
                else
                MeHgPrey(4) = MeHgPrey(3);
                end
            end 
       else
           
           for j = 1:length(Iprey)              
                if Iprey(j) == 1 %If prey is phytoplankton
                    MeHgPrey(j) = MeHgPicoP;
                elseif Iprey(j) == 2
                    MeHgPrey(j) = MeHgNanoP;
                elseif Iprey(j) == 3
                    MeHgPrey(j) = MeHgMicroP;
                else
                    if Iprey(j) == 4 % szp
                        massPreyArray = squeeze(szY(:,ti,1));
                        MeHgPreyArray = squeeze(szY(:,ti,2));
                    else
                        massPreyArray = squeeze(BigY(Iprey(j),:,ti,1));
                        MeHgPreyArray = squeeze(BigY(Iprey(j),:,ti,2));
                    end
                    aLPrey = aL(Iprey(j));
                    bLPrey = bL(Iprey(j));

                    lengthPreyArray = (massPreyArray/aLPrey).^(1/bLPrey);
                    lsr = Length_y/Size_ratio;
                    alivearray = lengthPreyArray > 0;
                    smallarray = lengthPreyArray<lsr;
                    yesarray = alivearray & smallarray;
                    if sum(yesarray)>0
                        [maxval,maxi] = max(lengthPreyArray.*yesarray);
                        MeHgPrey(j) = MeHgPreyArray(maxi);
                    else
                        X(j) = X(j) * EnerCont(Iprey(j))/EnerCont(4);
                        if strcmp(Species, 'Echinodermata')
                            MeHgPrey(j)=0;
                        else
                            MeHgPrey(j) = max(szY(:,ti,2));
                           
                        end
                    end
                end
            end

        end
%%
        MeHgEaten = (MeHgPrey*X);  % Intake of MeHg in ng-MeHg/g-fish/d

        MeHgAssim  = f_assim*MeHgEaten;      % ng-MeHg/g-fish/d

        Cox = (-0.24*Temp + 14.04)*sat;
        Gv  = 1400*((y(1)/1000)^0.65)/Cox;  % weight is in kg0 for this equation - Arnot and Gobas (2004)

        MeHgWaterUptake = Gv * MeHgWaterConc * neta / y(1);  % L/d * ng/L * unitless * (1/g-ZL)

        MeHgElim = aE * (y(1)^(bE)) * exp(cE*Temp) * y(2); 

        dy(2,1) = MeHgWaterUptake + MeHgAssim - MeHgElim - y(2) * (1/y(1)) * dy(1,1);  % ng-MeHg/g-fish/d
        
   
  
end


