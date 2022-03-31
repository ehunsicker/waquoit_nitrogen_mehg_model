function [ Windowed ] = get_agewindow( species, agemin, agemax, tinterest, OutputObj )
%GET_AGEWINDOW 
%  Gives the average M and MeHg for species between agemin and
%  agemax at time tinterest
% 
    tidays = tinterest*365.25;
    BigY = OutputObj.BigY;
    t0 = OutputObj.t0;
    dt =  OutputObj.dt;
    cyoung = find_cohort_from_age(agemin, t0, tidays);
    cold = find_cohort_from_age(agemax, t0, tidays);
    ti = int16((tidays-t0)/dt);
    if cold <= 0
        Windowed = [nan,nan];
    else
        Windows = squeeze(BigY(species,cold:cyoung,ti,:));
        if size(Windows)>1
            Windowed = nanmean(Windows,1);
        else
            Windowed = Windows;
        end
    end
end

