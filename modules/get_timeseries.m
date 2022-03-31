function [isXseries, Xinput, t_X] = get_timeseries(X_filename)
% GET_TIMESERIES
% retrieve timeseries from file named X_filename
if exist(X_filename,'file') == 2
    isXseries = 'Yes';
    X_series = xlsread(X_filename);
    Xinput = X_series(:,2);
    t_X = X_series(:,1).*365.25;
else
    isXseries = 'No';
    Xinput = 'None';
    t_X = 'None';
end

end