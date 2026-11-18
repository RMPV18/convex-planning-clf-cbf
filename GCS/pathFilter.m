function [new_path, new_wp, rm_i, n] =  pathFilter(obs, path, wp)
n = [];
for i = 2:size(wp,2)-1
    n = [n, false];
    for j=1:numel(obs)
        if obs{j}.contains(wp(:,i),1e-3)
            n(end) = true;
            continue;
        end
    end
end
n = [true boolean(n) true];
rm_i = find(~n);
new_wp = wp(:, n);
new_path = path(n);
end