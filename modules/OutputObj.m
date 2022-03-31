classdef OutputObj
    %OUTPUTOBJ 
    %  Holder for model output and metadata 
    
    properties
        BigY
        t0
        tend
        dt
        SpeciesNames
        cohortsperyear
        szY
        szcohortsperyear
	Ymean
	Ymin
	Ymax
	Ystd
    end
    
    methods
        function timespan = get_timespan(obj)
            s = size(obj.BigY);
            tsteps = s(3);
            timespan = zeros(tsteps,1);
            for i=0:tsteps-1
                timespan(i+1) = obj.t0 + i*obj.dt;
            end
        end
    end
    
end

