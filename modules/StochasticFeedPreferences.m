function Feed = StochasticFeedPreferences(FeedMaxs, FeedMins, FeedPrefs, seed)
% STOCHASTICFEEDPREFERENCES
%  
%  Create a sample of FeedPreferences given Max and Min values for each 
%  species, with the constraint that sum(Feed) == 1
rng(seed)
Feed = zeros(1,numel(FeedPrefs));
Iprey = find(FeedPrefs); %Index of prey organisms
[sortprey,sortind] = sort(FeedPrefs(Iprey));
sortpreyr = fliplr(sortprey);
sortindr = fliplr(sortind);
numprey = numel(sortpreyr) ;
sum=0;
for ff=1:numprey
    r = rand;
    pert = (1-r)*FeedMaxs(Iprey(sortindr(ff)))+r*FeedMins(Iprey(sortindr(ff)));
    Feed(Iprey(sortindr(ff))) = pert;
    sum = sum + pert;
end
leftovers = 1-sum; % how far from 1 did we end up?
if leftovers > 0 % if we don't have enough Feed, give to the smallest choice
    for ff=1:numprey
        feedfromlimit = FeedMaxs(Iprey(sortind(ff)))-Feed(Iprey(sortind(ff)));
        if leftovers > feedfromlimit
            Feed(Iprey(sortind(ff))) = FeedMaxs(Iprey(sortind(ff)));
            leftovers = leftovers - feedfromlimit;
        else
            Feed(Iprey(sortind(ff))) = Feed(Iprey(sortind(ff))) + leftovers;
            leftovers = 0;
        end
    end
elseif leftovers < 0 % if we have too much Feed, steal from the smallest choice
    for ff=1:numprey
        feedfromlimit = FeedMins(Iprey(sortind(ff)))-Feed(Iprey(sortind(ff)));
        if leftovers < feedfromlimit
            Feed(Iprey(sortind(ff))) = FeedMins(Iprey(sortind(ff)));
            leftovers = leftovers - feedfromlimit;
        else
            Feed(Iprey(sortind(ff))) = Feed(Iprey(sortind(ff))) + leftovers;
            leftovers = 0;
        end
    end
end