function [MeHgMicroP MeHgNanoP MeHgPicoP] = getPhytoplanktonHg(DOC,MeHgWaterConc)
%GETPHYTOPLANKTONHG
 % phytoplankton MeHg concentration in micro, nano, and pico plankton
 
% Cell Radius
Radius       = [1 11 110]; %[pico nano micro] um

%Volume in um3
VolumeOfCell=  4./3.*pi.*Radius.^3 ;  % um^3

%uptate rate amole/um3/h/nM
SAV=(3./Radius)-(3./Radius)*20/100;
Uptake =((0.1176*SAV)).*exp(-0.01.*DOC); %DOC in umol/L 


%MeHg in cell

%MeHgWaterConc should be in amole
MeHgWaterConc=MeHgWaterConc*1e6/200.59; %converting to fM
MeHgInCell =Uptake*(MeHgWaterConc./1000./1000).*4.*VolumeOfCell;

DensityOfCell = 1.00E-12; %ro, g um-3 marine phytoplankton density (Mason et al. 1996)
MassOfCell = DensityOfCell.*VolumeOfCell; %g

%Concentraiton in Cell pmol / g
FinalConcentration = MeHgInCell./MassOfCell./(1000*1000);

MeHgMicroP = FinalConcentration(3)*200.59/1000; %converts pmol/g to ng/g 
MeHgNanoP = FinalConcentration(2)*200.59/1000;  
MeHgPicoP = FinalConcentration(1)*200.59/1000;  

end