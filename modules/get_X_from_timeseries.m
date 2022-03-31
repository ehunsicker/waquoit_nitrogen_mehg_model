function X_at_time = get_X_from_timeseries(t,Xinput,tinput)
% GET_X_FROM_TIMESERIES
% get the value of X corresponding to time t according to the time
% series Tinput
X_at_time = -999e20;
for i = 1:numel(tinput)-1
    if (tinput(i) < t) && (t < tinput(i+1))
        X_at_time = Xinput(i);
    end
end
if X_at_time < -998e20
    X_at_time = Xinput(end); % at end of timeseries, reuse last point
end
