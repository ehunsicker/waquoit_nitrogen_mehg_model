function [ CohortArray ] = get_cohort( species, birth_year, OutputObj )
%GET_COHORT
%  Gets desired cohort info from a BigY-type array
%  Give a cohort birth year and species, get M,MeHg across its lifespan
    BigY = OutputObj.BigY;
    t0 = OutputObj.t0;
    cpy = OutputObj.cohortsperyear;
    ci = int16((birth_year*365.25-t0)*cpy/365.25)+1;
    Cohorts = squeeze(BigY(species,ci:ci+cpy-1,:,:));
    if cpy>1
        CohortArray = nanmean(Cohorts,1);
    else
        CohortArray = Cohorts(:,:);
    end
end

