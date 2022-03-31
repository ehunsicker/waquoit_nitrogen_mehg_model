function [ c ] = find_cohort_from_age( age, t0, t )
%FIND_COHORT_FROM_AGE 
% find which cohort is age=age at time=t
    cohorts_since_t0 = int16((t-t0)/365.25);
    integer_age = int16(age);
    c = cohorts_since_t0-integer_age;
end

