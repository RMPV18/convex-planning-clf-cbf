function regions = sortRegions(regions, ref, criteria)
    for k = 1:numel(regions) cnv(k) = vecnorm(mean(regions{k}, 1)-ref');end
    [~, ind_c] = sort(cnv);
    if  nargin < 3
        regions = regions(ind_c);
    elseif strcmp(criteria,'closervertex') 
        for k = 1:numel(regions) [mnv(k), ~] = min(vecnorm(regions{k}-ref', 2, 2)); end
        [~, ind_m] = sort(mnv(ind_c));
        regions = regions(ind_c(ind_m));
    end
end