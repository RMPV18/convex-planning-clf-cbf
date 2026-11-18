function clean_path = clear_path(path)
    % path{1} and path{end} are ignored, matching the original indexing
    p = path(2:end-1);
    e = cell2mat(cellfun(@(s) sscanf(s, '%d-%d')', p(:), 'UniformOutput', false));
    m = size(e, 1);
    assert(m >= 2, 'need at least 2 edges between path{1} and path{end}');
    nodes = zeros(1, m+1);
    nodes(1) = setdiff(e(1,:), e(2,:));
    for k = 1:m-1
        nodes(k+1) = intersect(e(k,:), e(k+1,:));
    end
    nodes(m+1) = setdiff(e(m,:), e(m-1,:));
    clean_path = cellstr([string(nodes)]);
end