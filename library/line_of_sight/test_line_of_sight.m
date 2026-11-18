function test_line_of_sight()
%TEST_LINE_OF_SIGHT  Verification of line_of_sight / prepare_los_model.
%   1) Deterministic 2-D case with a non-convex (L-shaped) obstacle.
%   2) Deterministic 3-D cube cases (touching, crossing, endpoint inside).
%   3) Random 3-D and 5-D: closed-form vs LP agreement.
%   4) Independent dense-sampling check of the closed-form decision.
%   Throws error() on any failure.

    rng(1);

    % ---- 1) L-shaped non-convex obstacle = union of two convex boxes ----
    boxA = [0 0; 4 0; 4 1; 0 1];
    boxB = [0 1; 1 1; 1 4; 0 4];
    Lshape = {boxA, boxB};
    obs = {Lshape};
    % Inside the concavity of the L: free
    chk(line_of_sight([2 2], [2 3.5], obs), true,  '2D free inside concavity');
    % Cut across the inner corner region: free (concavity is empty)
    chk(line_of_sight([3 3], [0.5 0.5] + [0 4], obs), true, '2D free above L');
    % Crossing arm A
    chk(line_of_sight([2 -1], [2 2], obs), false, '2D blocked by arm A');
    % Crossing arm B
    chk(line_of_sight([-1 3], [2 3], obs), false, '2D blocked by arm B');
    % Convex hull of the L would block this; the true L does not
    chk(line_of_sight([3 3], [1.5 1.5] , obs), true, '2D hull would block, L does not');
    % Through arm A along its mid-line
    chk(line_of_sight([-1 0.5], [5 0.5], obs), false, '2D passes through arm A interior');
    chk(line_of_sight([1 1.0], [4 4], {boxA}), true, '2D touches corner only');

    % ---- 2) 3-D cube ----
    [X, Y, Z] = ndgrid([0 1], [0 1], [0 1]);
    cube = [X(:) Y(:) Z(:)];
    chk(line_of_sight([-1 .5 .5], [2 .5 .5], {cube}), false, '3D through cube');
    chk(line_of_sight([-1 .5 .5], [-1 .5 2], {cube}), true,  '3D beside cube');
    chk(line_of_sight([.5 .5 .5], [3 3 3], {cube}),    false, '3D endpoint inside');
    chk(line_of_sight([-1 0 0], [2 0 0], {cube}),      true,  '3D slides along edge');
    chk(line_of_sight([-1 -1 .5], [2 2 .5], {cube}),   false, '3D diagonal through');
    chk(line_of_sight([0 0 0], [0 0 0], {cube}),       true,  '3D degenerate segment on vertex');

    % ---- 3) + 4) random agreement tests ----
    run_random(3, 25, 1500, 1500);
    run_random(5, 12, 300, 600);

    % ---- flat obstacle (empty interior) fallback ----
    flat = [0 0 0; 1 0 0; 1 1 0; 0 1 0];
    [l, info] = line_of_sight([.5 .5 -1], [.5 .5 1], {flat}, struct('method', 'closed'));
    chk(l, false, '3D flat polygon crossed');
    chk(info.degenerateFallback, true, 'degenerate fallback flagged');

    fprintf('test_line_of_sight: all checks passed.\n');
end

function run_random(n, ncomp, nseg, nsamp)
    comps = cell(1, ncomp);
    for i = 1:ncomp
        c = 10 * rand(1, n);
        r = 0.4 + 0.8 * rand;
        comps{i} = c + r * (2 * rand(2 * n + 4, n) - 1);
    end
    model = prepare_los_model({comps}, n);

    nDis = 0; nBlocked = 0; nMiss = 0; nBad = 0;
    t = linspace(0, 1, nsamp);
    for s = 1:nseg
        p = 10 * rand(n, 1);
        q = 10 * rand(n, 1);
        d = q - p;
        [lc, ~] = line_of_sight(p, q, model, struct('method', 'closed'));
        ll      = line_of_sight(p, q, model, struct('method', 'lp'));
        if lc ~= ll, nDis = nDis + 1; end
        nBlocked = nBlocked + ~lc;

        % Independent sampling: any sample strictly inside a hull => must be blocked
        P = p + d * t;
        inside = false;
        for k = 1:numel(model)
            if any(all(model(k).A * P < model(k).b - 1e-6, 1))
                inside = true; break;
            end
        end
        if inside && lc, nBad = nBad + 1; end        % closed says free, sample inside: error
        if ~inside && ~lc, nMiss = nMiss + 1; end    % thin crossing missed by sampling: expected
    end
    fprintf('n=%d: segments=%d blocked=%d | closed/LP disagreements=%d | sampling contradictions=%d | thin crossings missed by sampling=%d\n', ...
            n, nseg, nBlocked, nDis, nBad, nMiss);
    if nBad > 0
        error('test_line_of_sight:samplingContradiction', 'Closed-form declared free but a sample lies inside (n=%d).', n);
    end
    if nDis > 0
        error('test_line_of_sight:methodDisagreement', 'Closed-form and LP disagree (n=%d).', n);
    end
end

function chk(got, expected, name)
    if got ~= expected
        error('test_line_of_sight:fail', 'FAILED: %s (got %d, expected %d)', name, got, expected);
    end
end
