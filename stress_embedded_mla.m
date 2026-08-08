% ============================================================
%  stress_embedded_mla.m
%
%  PER-FAMILY STRESS TEST — EMBEDDED (PLANARIZED) MLA
%
%  [Purpose] The freeform-MLA study showed "practical saturation":
%    efficiency and every angular band rise together, no angular steering,
%    selectivity pinned near the Lambertian partition
%    S = sin^2(th2) - sin^2(th1) = [0.117 0.296 0.337 0.220].
%  This script tests whether that saturation generalizes to the EMBEDDED
%  MLA family: the freeform lens array is covered/planarized by a
%  low-index resin overcoat with a flat top surface (cf. Qu et al.,
%  ACS Photonics 2018). The reduced index contrast at the lens interface
%  (lens/resin instead of lens/air) is the ONE physics change; source,
%  substrate, cavity, detectors and band readout are identical to
%  pareto_front_freeform.m (the verified reference implementation).
%
%  [Protocol — lightweight per-family budget]
%    1) N_RANDOM = 100 random valid designs at search fidelity.
%    2) ONE single-objective EQE_total optimization:
%       surrogateopt (MaxFunctionEvaluations=60, MinSurrogatePoints=25)
%       + patternsearch polish (15 evals), then high-precision
%       re-evaluation of the winner (N_FINAL_REP=3).
%    3) Outputs:
%       stress_embedded_result.mat — EVAL_LOG in the standard column format
%         [ x(1:13) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w ]
%         (saved incrementally, crash-safe)
%       stress_embedded_check.png — 3-panel signature check:
%         (a) EQE_40_60 vs EQE_total scatter + linear fit  (linear collapse?)
%         (b) best-design selectivity S_j vs Lambertian partition (bars)
%         (c) corr R(EQE_total, S_j) per band vs freeform prior
%             [+0.6 +0.7 +0.05 -0.7]                        (drift signature?)
%    4) All evaluations from this script are tagged EVAL_PHASE=5, EVAL_W=-1.
%
%  Base: pareto_front_freeform.m (verified LightTools plumbing / geometry /
%        stack / GEOM_TOL round-trip / multi-fidelity constants).
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current EVAL_LOG EVAL_PHASE EVAL_W ...
       GEOM_TOL GEOM_MISMATCH_LOG REQUIRE_MONOTONIC_X N_PLANAR T_OVER

% [Geometry verification tolerance] LightTools control-point round-trip
% mismatch allowance (same as reference).
GEOM_TOL = 1e-4;
GEOM_MISMATCH_LOG = [];   % columns: [mismatch, max_length, rescale_triggered]

% [Yield] require monotonic x2..x6 (sorted) — non-monotonic profiles make
% the spline overshoot x>1, triggering rescale -> reset mismatch -> NaN.
REQUIRE_MONOTONIC_X = true;

%% ===== THE ONE PHYSICS CHANGE: planarization (embedded MLA) =====
N_PLANAR = 1.41;    % resin refractive index (low-index overcoat)
T_OVER   = 0.005;   % [mm] resin thickness ABOVE the lens apex (flat top)

RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
restart_interval = 20;
lt = ltloc.GetLTAPI(ID_swept);
ltx= getltpointer(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== Multi-fidelity (same constants as reference) =====
% [Note] wavelength indices come from wv_list = 1:n:wavelength_num — never
% division-based sizing, so any (window, n) combination is safe.
WAVE_N_SEARCH = 10;      % wavelength step, search fidelity
WAVE_N_FINAL  = 2;       % wavelength step, verification fidelity
RAY_SEARCH    = 10000;
RAY_FINAL     = 50000;
N_FINAL_REP   = 3;

%% ===== Stress-test budget =====
N_RANDOM        = 100;   % random valid designs (phase A)
OPT_EVALS       = 60;    % surrogateopt MaxFunctionEvaluations
MIN_SURR_POINTS = 25;
N_SEED_VALID    = 20;    % valid seeds for surrogateopt
POLISH_EVALS    = 15;    % patternsearch polish budget

% Standard phase/weight markers for this family (embedded MLA stress test)
EVAL_PHASE = 5;
EVAL_W     = -1;

% Lambertian reference partition and freeform-baseline drift correlations
S_LAMBERT   = [0.117, 0.296, 0.337, 0.220];   % sin^2(th2)-sin^2(th1)
R_FREEFORM  = [0.6, 0.7, 0.05, -0.7];         % prior corr(EQE_total, S_j)
band_names  = {'0-20\circ','20-40\circ','40-60\circ','60-80\circ'};

%% Optimization Variables (13-dim, identical parameterization)
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];
nvar = numel(lb);

EVAL_LOG = [];   % [x(1:13) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w]

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, 'MeshTolerance', 1e-3, ...
    'Cache', 'on', 'Display', 'off');

%% =====================================================================
%  PHASE A — random valid designs (search fidelity)
%% =====================================================================
fprintf('\n########## EMBEDDED MLA STRESS: random designs (N=%d) ##########\n', N_RANDOM);
ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;
eval_count = 0;
Prand = genValidPoints(N_RANDOM, lb, ub);
for i = 1:N_RANDOM
    [et, eb] = simulate_both(Prand(i,:));
    if mod(i,10)==0
        fprintf('  random %3d/%d : EQE_total=%.4f  EQE_4060=%.4f\n', i, N_RANDOM, et, eb);
        save('stress_embedded_result.mat','EVAL_LOG','varNames','lb','ub', ...
             'N_PLANAR','T_OVER');   % crash-safe incremental save
    end
end
save('stress_embedded_result.mat','EVAL_LOG','varNames','lb','ub','N_PLANAR','T_OVER');

% --- geometry-rejection diagnostics (main NaN cause) ---
report_geom_rejection(GEOM_MISMATCH_LOG, GEOM_TOL);

%% =====================================================================
%  PHASE B — ONE single-objective EQE_total optimization
%% =====================================================================
fprintf('\n########## EMBEDDED MLA STRESS: EQE_total optimization ##########\n');
RenewLightTools();
lt = ltloc.GetLTAPI(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
eval_count = 0;
ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;

seedMat = genValidPoints(N_SEED_VALID, lb, ub);
sopts = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', OPT_EVALS, ...
    'MinSurrogatePoints',     MIN_SURR_POINTS, ...
    'InitialPoints',          struct('X', seedMat), ...
    'UseParallel', false, 'PlotFcn', [], 'Display', 'iter');
[xS, ~] = surrogateopt(@eqe_objconstr, lb, ub, sopts);

save('stress_embedded_result.mat','EVAL_LOG','varNames','lb','ub','N_PLANAR','T_OVER');

% --- patternsearch polish ---
bestX = [];
if ~isempty(xS) && isValidPoints(xS(:).')
    x0 = xS(:).';
    try
        xP = patternsearch(@eqe_polish, x0, [],[],[],[], lb, ub, [], psOpts);
        xP = xP(:).';
    catch
        xP = x0;
    end
    % pick the better of {x0, xP} at search fidelity
    cands = {x0};
    if ~isequal(xP, x0), cands{end+1} = xP; end
    bestScore = -inf;
    for c = 1:numel(cands)
        [et, ~] = simulate_both(cands{c});
        if isfinite(et) && et > bestScore, bestScore = et; bestX = cands{c}; end
    end
else
    fprintf('  [Warn] surrogateopt returned no feasible solution.\n');
end
save('stress_embedded_result.mat','EVAL_LOG','varNames','lb','ub','N_PLANAR','T_OVER');

%% =====================================================================
%  PHASE C — high-precision re-evaluation of the winner
%% =====================================================================
best_tot  = NaN;  best_bins = nan(1,4);
if ~isempty(bestX)
    fprintf('\n########## High-precision re-evaluation (N=%d) ##########\n', N_FINAL_REP);
    ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
    et_r  = nan(1,N_FINAL_REP);
    bin_r = nan(N_FINAL_REP,4);
    for r = 1:N_FINAL_REP
        [et_r(r), ~] = simulate_both(bestX);
        if isfinite(et_r(r))
            bin_r(r,:) = EVAL_LOG(end, nvar+2 : nvar+5);
        end
    end
    best_tot  = mean(et_r,'omitnan');
    best_bins = mean(bin_r,1,'omitnan');
    fprintf('  >>> best embedded design: EQE_total = %.5f ± %.5f (%d rays, N=%d)\n', ...
        best_tot, std(et_r,'omitnan'), RAY_FINAL, N_FINAL_REP);
end

save('stress_embedded_result.mat','EVAL_LOG','varNames','lb','ub', ...
     'N_PLANAR','T_OVER','bestX','best_tot','best_bins','S_LAMBERT','R_FREEFORM');
fprintf('\nsaved -> stress_embedded_result.mat  (EVAL_LOG %d points)\n', size(EVAL_LOG,1));
report_geom_rejection(GEOM_MISMATCH_LOG, GEOM_TOL);

%% =====================================================================
%  ANALYSIS — three saturation signatures + verdict
%  EVAL_LOG = [x(1:nvar) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w]
%% =====================================================================
nv   = nvar;
Et   = EVAL_LOG(:,nv+1);
Bins = EVAL_LOG(:, nv+2 : nv+5);
ok   = isfinite(Et) & Et > 0.05 & all(isfinite(Bins),2);   % low-EQE = noise-dominated
fprintf('\n################ EMBEDDED MLA: signature check ################\n');
fprintf('valid designs: %d of %d evaluations\n', sum(ok), size(EVAL_LOG,1));

% --- Signature 1: near-linear collapse of band vs total (40-60 band) ---
p_fit = polyfit(Et(ok), Bins(ok,3), 1);
resid = Bins(ok,3) - polyval(p_fit, Et(ok));
R2    = 1 - sum(resid.^2) / sum((Bins(ok,3) - mean(Bins(ok,3))).^2);

% --- Signature 2: best-design selectivity vs Lambertian partition ---
if all(isfinite(best_bins)) && isfinite(best_tot) && best_tot > 0
    S_best = best_bins / best_tot;
else
    % fall back to the highest-EQE logged design
    [~, ibest] = max(Et .* ok);
    S_best = Bins(ibest,:) / Et(ibest);
end
dev_best = (S_best - S_LAMBERT) ./ S_LAMBERT;

% --- Signature 3: drift correlation R(EQE_total, S_j) per band ---
R_corr = nan(1,4);
for b = 1:4
    s = Bins(ok,b) ./ Et(ok);
    c = corrcoef(Et(ok), s);
    R_corr(b) = c(1,2);
end

%% ===== 3-panel figure =====
figure('Name','Embedded MLA stress check','Color','w','Position',[80 80 1400 420]);

subplot(1,3,1);   % (a) linear collapse, 40-60 band
scatter(Et(ok), Bins(ok,3), 14, [.2 .5 .8], 'filled', 'MarkerFaceAlpha', 0.5); hold on;
xg = linspace(min(Et(ok)), max(Et(ok)), 50);
plot(xg, polyval(p_fit, xg), 'r-', 'LineWidth', 2);
xlabel('EQE_{total}'); ylabel('EQE_{40-60}'); grid on;
title(sprintf('(a) linear collapse: R^2 = %.3f', R2));
legend({'designs', sprintf('fit: slope %.3f', p_fit(1))}, 'Location','northwest','FontSize',8);

subplot(1,3,2);   % (b) best-design selectivity vs Lambertian
bar([S_best(:), S_LAMBERT(:)]); grid on;
set(gca,'XTickLabel',band_names);
ylabel('selectivity  S_j = b_j / EQE_{total}');
legend({'best embedded design','Lambertian'}, 'Location','northwest','FontSize',8);
title('(b) selectivity vs Lambertian partition');

subplot(1,3,3);   % (c) drift correlation vs freeform prior
bar([R_corr(:), R_FREEFORM(:)]); grid on;
set(gca,'XTickLabel',band_names);  ylim([-1 1]);
ylabel('corr R (EQE_{total}, S_j)');
legend({'embedded MLA','freeform prior'}, 'Location','southwest','FontSize',8);
title('(c) drift signature per band');

saveas(gcf,'stress_embedded_check.png');
fprintf('saved -> stress_embedded_check.png\n');

%% ===== text verdict =====
fprintf('\n################ VERDICT: embedded vs freeform baseline ################\n');
fprintf('[1] Linear collapse (40-60 band):  R^2 = %.3f, slope = %.3f\n', R2, p_fit(1));
fprintf('    => %s\n', ternary(R2 > 0.9, ...
    'near-linear collapse PRESENT (matches freeform saturation)', ...
    'collapse WEAK/ABSENT — embedded family deviates from freeform baseline'));

fprintf('[2] Best-design selectivity vs Lambertian [0.117 0.296 0.337 0.220]:\n');
for b = 1:4
    fprintf('    %-8s S = %.3f  (Lambertian %.3f, deviation %+.1f%%)\n', ...
        band_names{b}, S_best(b), S_LAMBERT(b), 100*dev_best(b));
end
fprintf('    => %s\n', ternary(all(abs(dev_best) < 0.15), ...
    'selectivity pinned to Lambertian partition (<15%% dev, matches freeform)', ...
    'selectivity DEVIATES >15%% in at least one band — check for steering'));

fprintf('[3] Drift correlation R(EQE_total, S_j) vs freeform prior [+0.6 +0.7 +0.05 -0.7]:\n');
for b = 1:4
    fprintf('    %-8s R = %+.2f  (freeform %+.2f)\n', band_names{b}, R_corr(b), R_FREEFORM(b));
end
same_sign = sign(R_corr) == sign(R_FREEFORM) | abs(R_corr) < 0.2 | abs(R_FREEFORM) < 0.2;
fprintf('    => %s\n', ternary(all(same_sign), ...
    'drift signature qualitatively matches the freeform baseline', ...
    'drift signature DIFFERS from freeform — family-specific behavior'));

overall = (R2 > 0.9) && all(abs(dev_best) < 0.15) && all(same_sign);
fprintf(['\n[OVERALL] %s\n'], ternary(overall, ...
    'All three signatures match: practical saturation GENERALIZES to the embedded MLA family.', ...
    'At least one signature differs: saturation does NOT fully generalize — inspect the failing panel.'));

fprintf('\n########## done ##########\n');
fprintf('  stress_embedded_result.mat / stress_embedded_check.png\n');


%% ===== evaluation wrapper (logs all four angular bands) =====
%  EVAL_LOG columns:
%    [ x(1:13) | EQE_total | EQE_0_20 | EQE_20_40 | EQE_40_60 | EQE_60_80 | phase | w ]
function [eqe_total, eqe_band] = simulate_both(pt)
global ID_swept ltml ltloc eval_count restart_interval EVAL_LOG EVAL_PHASE EVAL_W
eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    fprintf('\n[Refresh] %d sims done. Restarting LightTools...\n', eval_count);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    pause(2);
end
bins = [NaN NaN NaN NaN];
try
    out = objFcn_both(pt);
    eqe_total = out.EQE_total;
    eqe_band  = out.EQE_40_60;
    bins = [out.EQE_0_20, out.EQE_20_40, out.EQE_40_60, out.EQE_60_80];
    if eqe_total == 0
        eqe_total = NaN; eqe_band = NaN; bins = [NaN NaN NaN NaN];
    end
catch err
    fprintf('\n[Error] eval %d LightTools crash: %s\n', eval_count, err.message);
    eqe_total = NaN;  eqe_band = NaN;
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
EVAL_LOG(end+1,:) = [pt(:).', eqe_total, bins, EVAL_PHASE, EVAL_W];
% crash-safe: append-style incremental save every 10 evaluations
if mod(size(EVAL_LOG,1), 10) == 0
    try
        save('stress_embedded_result.mat','EVAL_LOG','-append');
    catch
        save('stress_embedded_result.mat','EVAL_LOG');
    end
end
end

%% ===== single-objective EQE_total (surrogateopt, constraint-coupled) =====
function out = eqe_objconstr(x)
global REQUIRE_MONOTONIC_X
x = x(:).';
% Non-monotonic x cannot be represented by LightTools anyway — reject at the
% constraint stage so surrogateopt does not waste budget on NaN.
if ~isempty(REQUIRE_MONOTONIC_X) && REQUIRE_MONOTONIC_X && any(diff(x(1:5)) < 0)
    out.Ineq = 1;  out.Fval = 1;  return;
end
if ~isValidPoints(x)
    out.Ineq = 1;  out.Fval = 1;  return;   % infeasible: return without sim
end
[et, ~] = simulate_both(x);
if ~isfinite(et)
    out.Ineq = 1;  out.Fval = 1;
else
    out.Ineq = -1;
    out.Fval = -et;                          % maximize EQE_total
end
end

%% ===== single-objective EQE_total (patternsearch) =====
function f = eqe_polish(x)
global REQUIRE_MONOTONIC_X
x = x(:).';
if ~isempty(REQUIRE_MONOTONIC_X) && REQUIRE_MONOTONIC_X && any(diff(x(1:5)) < 0)
    f = 0; return;
end
if ~isValidPoints(x), f = 0; return; end
[et, ~] = simulate_both(x);
if ~isfinite(et), f = 0; return; end
f = -et;
end

%% ===== verdict-string helper =====
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

%% ===== geometry-rejection diagnostics =====
function report_geom_rejection(mismLog, tol)
if isempty(mismLog), return; end
mism = mismLog(:,1);  maxlen = mismLog(:,2);  resc = mismLog(:,3) > 0;
rej = mism > tol;
fprintf('\n--- geometry rejection diagnostics (NaN cause) ---\n');
fprintf('  %d evals, %d rejected (%.1f%%), tol=%.1e\n', ...
    numel(mism), sum(rej), 100*mean(rej), tol);
if any(resc) || any(~resc)
    r1 = mean(rej(resc));   n1 = sum(resc);
    r0 = mean(rej(~resc));  n0 = sum(~resc);
    fprintf('  rescale triggered (curve x>1): %d times (%.1f%%), rejection rate %.1f%%\n', ...
        n1, 100*mean(resc), 100*r1);
    fprintf('  rescale not triggered        : %d times, rejection rate %.1f%%\n', n0, 100*r0);
    if any(maxlen > 1)
        fprintf('  overshoot size: median max_length=%.3f, max=%.3f\n', ...
            median(maxlen(maxlen>1)), max(maxlen));
    end
end
if any(rej)
    r = mism(rej);
    fprintf('  mismatch on rejection: median=%.2e, p90=%.2e, max=%.2e\n', ...
        median(r), prctile(r,90), max(r));
end
if any(~rej)
    fprintf('  mismatch on pass: median=%.2e, max=%.2e\n', ...
        median(mism(~rej)), max(mism(~rej)));
end
end

%% ===== random valid seed generation =====
%  With REQUIRE_MONOTONIC_X, x2..x6 are sorted ascending on generation.
function P = genValidPoints(K, lb, ub)
global REQUIRE_MONOTONIC_X
mono = ~isempty(REQUIRE_MONOTONIC_X) && REQUIRE_MONOTONIC_X;
dim = numel(lb);  P = zeros(K, dim);
for i = 1:K
    ok = false;
    while ~ok
        p = lb + rand(1, dim) .* (ub - lb);
        if mono
            p(1:5) = sort(p(1:5));          % x2..x6 ascending
        end
        if isValidPoints(p), ok = true; P(i, :) = p; end
    end
end
end


%% ===== Objective (EQE_total + four angular bands) — EMBEDDED MLA =====
%  Identical to pareto_front_freeform.m::objFcn_both, EXCEPT for the
%  planarization block (marked below).
function output = objFcn_both(point)
global ID_LT ID_swept ltml ltloc count ray_nums_current wave_n_current N_PLANAR T_OVER
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

d_sub=1.295;  r_OLED=1;  x_pattern=15;  y_pattern=15;  Lensheight=0.01;
wavelength_start=453;  wavelength_end=753;

if isempty(wave_n_current), n = 10;    else, n = wave_n_current;    end
if isempty(ray_nums_current), ray_nums = 10000; else, ray_nums = ray_nums_current; end

List=ltml.LTDbList(lt,'lens_manager[1]','SIMULATIONS');
Key=ltml.LTListByName(lt,List,'ForwardAll');
ltml.LTDbSet(lt,Key,'MaxProgress',ray_nums);
List=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
Key=ltml.LTListByName(lt,List,'Substrate');
ltml.LTDbSet(lt,Key,'Height',d_sub);
ltml.LTDbSet(lt,Key,'Y',d_sub/2);
SRList=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
SRKey=ltml.LTListAtPos(lt,SRList,2);
ltml.LTDbSet(lt,SRKey,'Y',d_sub+Lensheight/2);
List=ltml.LTDbList(lt,'lens_manager[1]','TEXTURE_ZONE_EXTENT');
Key=ltml.LTListByName(lt,List,'zone');
ltml.LTDbSet(lt,Key,'Geometry_1',x_pattern);
ltml.LTDbSet(lt,Key,'Geometry_2',y_pattern);
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

x2 = point(1);  x3 = point(2);  x4 = point(3);  x5 = point(4);  x6 = point(5);
y2 = point(6);  y3 = point(7);  y4 = point(8);  y5 = point(9);  y6 = point(10);
dETL = point(11); dHTL = point(12); stretchZ=point(13);

xy = zeros(7,2);
xy(1,:) = [0, 1];  xy(7,:) = [1, 0];
xy(2,:) = [x2, y2];  xy(3,:) = [x3, y3];  xy(4,:) = [x4, y4];
xy(5,:) = [x5, y5];  xy(6,:) = [x6, y6];

lt = ltloc.GetLTAPI(ID_swept);
ltx= getltpointer(ID_swept);
lt2 = ltloc.GetLTAPI(ID_LT);

Curve="LENS_MANAGER[1].COMPONENTS[Components].SWEPT_SOLID[SweptEntity].SWEPT_PRIMITIVE[SweptPrimitive].SWEPT_PROFILE[SweptProfile].FITTED_CURVE[SweptSurface_1]";
ltx.SetSweptProfilePoints(Curve,xy,7);
ltx.DbSet(Curve,'StartSlopeMode',"Auto");
ltx.DbSet(Curve,'EndSlopeMode',"Auto");

List=ltml.LTDbList(lt,'LENS_MANAGER[1]','FITTED_CURVE');
Key=ltml.LTListByName(lt,List,'SweptSurface_1');
ltml.LTDbSet(lt, Key,'NumFacets',100);
x_values = zeros(101,1);
for a=1:101
    x_values(a)=ltml.LTDbGet(lt,Key,'YFacetsAt',a);
end
max_length = max(x_values);
rescaled = (max_length > 1);        % rescale trigger flag (for diagnostics)
if rescaled
    xy = xy / max_length;
end
ltx.SetSweptProfilePoints(Curve,xy,7);
ltx.DbSet(Curve,'StartSlopeMode',"Auto");
ltx.DbSet(Curve,'EndSlopeMode',"Auto");

xy_l = zeros(7,2);
for j=1:7
    xy_l(j,1) = ltml.LTDbGet(lt, Key, 'YAt', j);
    xy_l(j,2) = ltml.LTDbGet(lt, Key, 'ZAt', j);
end
% [Geometry verification] round-trip check of control points set vs read.
% Large mismatch = the intended shape is not realized -> reject (EQE=0 -> NaN).
global GEOM_TOL GEOM_MISMATCH_LOG
if isempty(GEOM_TOL), GEOM_TOL = 1e-4; end
mism = max(abs(xy(:) - xy_l(:)));

% possible COM round-trip glitch -> reset once and re-check
if mism > GEOM_TOL
    ltx.SetSweptProfilePoints(Curve,xy,7);
    ltx.DbSet(Curve,'StartSlopeMode',"Auto");
    ltx.DbSet(Curve,'EndSlopeMode',"Auto");
    for j=1:7
        xy_l(j,1) = ltml.LTDbGet(lt, Key, 'YAt', j);
        xy_l(j,2) = ltml.LTDbGet(lt, Key, 'ZAt', j);
    end
    mism = max(abs(xy(:) - xy_l(:)));
end

GEOM_MISMATCH_LOG(end+1,:) = [mism, max_length, double(rescaled)];
if mism > GEOM_TOL
    output = struct('EQE_0_20',0,'EQE_20_40',0,'EQE_40_60',0,'EQE_60_80',0,'EQE_total',0);
    return;
end

rng('shuffle')
charSet = ['a':'z' 'A':'Z' '0':'9'];
index = charSet(randi(length(charSet), 1, 10));
pathname = '"C:\Users\jhkim\Desktop\Green_CE_Calculation\swept_';
pathname_unrepaired = '"C:\Users\jhkim\Desktop\Green_CE_Calculation\unrepaired\swept_unrepaired_';
totalpath = [pathname index '.ent"'];
totalpath_unrepaired = [pathname_unrepaired index '.ent"'];

ltml.LTCmd(lt, 'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt, sprintf('SaveLibrary XYZ 0,0,0 %s ', totalpath_unrepaired));
ltml.LTCmd(lt, 'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt, 'RepairEntities');
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
ltml.LTSetOption(lt2, "ShowFileDialogBox", 0);
ltml.LTCmd(lt, 'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt, sprintf('SaveLibrary XYZ 0,0,0 %s ', totalpath));
ltml.LTCmd(lt, 'Undo');
ltml.LTCmd(lt, 'Undo');

totalpathmod = [pathname index '.1.ent"'];
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'LIBRARY_ELEMENT_UNIT_CELL');
Key = ltml.LTListByName(lt2, List, 'LibraryElement');
ltml.LTDbSet(lt2, Key, 'Filename', totalpathmod);
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'TEXTURE_PARAMETER');
Key = ltml.LTListByName(lt2, List, 'StretchZ');
ltml.LTDbSet(lt2, Key, 'Value', stretchZ);

%% !!! VERIFY IN LIGHTTOOLS: planarization block creation adapted from
%%     pareto_front_freeform.m::objFcn_both CUBE_PRIMITIVE modify pattern
%%     (Substrate Height/Y sets) — check material assignment and z-extent
%%     before trusting results !!!
%
%  EMBEDDED MLA — the one physics change. A flat resin overcoat
%  (index N_PLANAR) covers the lens array from the lens base plane
%  (y = d_sub) up to a flat top at y = d_sub + Lensheight*stretchZ + T_OVER.
%  Resulting optical interfaces:  lens material / resin at the lens surface,
%  resin / air at the flat top.
%
%  NO verified precedent exists in this repo for CREATING a solid via COM —
%  every verified script only MODIFIES pre-existing CUBE_PRIMITIVEs of the
%  .lts model (by name 'Substrate', or by list position 2 for the texture
%  host cube). The code below therefore:
%   (1) looks for a cube named 'PlanarOvercoat' (add it ONCE in the GUI to
%       the array model, recommended path), and if found only updates its
%       Height / Y each evaluation using the verified LTDbSet pattern;
%   (2) if absent, attempts one-time creation via the LightTools command
%       'Cuboid' + rename + Width/Depth set. Commands in branch (2) are
%       UNVERIFIED — validate on first run (check the 3D view: block must
%       span the full x_pattern x y_pattern patch, sit on y=d_sub, and its
%       material must be a simple index n=1.41; also confirm LightTools
%       immersion rules give lens/resin at the lens surface).
%
%  z-extent logic (vertical axis is Y in this model, units mm):
%    lens base plane : y = d_sub
%    lens apex       : y = d_sub + Lensheight*stretchZ   (unit-cell height
%                      0.01 mm scaled by texture parameter StretchZ —
%                      VERIFY that StretchZ scales geometric height 1:1)
%    flat resin top  : y = d_sub + Lensheight*stretchZ + T_OVER
%    => cube Height  = Lensheight*stretchZ + T_OVER
%       cube Y (ctr) = d_sub + (Lensheight*stretchZ + T_OVER)/2
t_planar = Lensheight*stretchZ + T_OVER;
PList = ltml.LTDbList(lt2,'lens_manager[1]','CUBE_PRIMITIVE');
PKey  = [];
try
    PKey = ltml.LTListByName(lt2, PList, 'PlanarOvercoat');
catch
    PKey = [];
end
% LTListByName may return 0 / empty for a missing name depending on API
% version — treat both as "not found".
notFound = isempty(PKey) || (isnumeric(PKey) && all(PKey(:) == 0));
if notFound
    % --- one-time creation (UNVERIFIED branch — see banner above) ---
    ltml.LTCmd(lt2, sprintf('Cuboid %g,%g,%g', x_pattern, t_planar, y_pattern));
    PList = ltml.LTDbList(lt2,'lens_manager[1]','CUBE_PRIMITIVE');
    % newly created cube is expected at the end of the list; the verified
    % model contains 2 cubes (Substrate + texture host), so position 3.
    PKey  = ltml.LTListAtPos(lt2, PList, 3);
    ltml.LTDbSet(lt2, PKey, 'Name', 'PlanarOvercoat');
    ltml.LTDbSet(lt2, PKey, 'Width',  x_pattern);   % X extent = texture patch
    ltml.LTDbSet(lt2, PKey, 'Depth',  y_pattern);   % Z extent = texture patch
    ltml.LTDbSet(lt2, PKey, 'X', 0);
    ltml.LTDbSet(lt2, PKey, 'Z', 0);
    % material: simple refractive index N_PLANAR, no absorption.
    % UNVERIFIED property path — if this errors, assign the material once in
    % the GUI (user-defined material, n = N_PLANAR) and keep only the
    % Height/Y updates below.
    try
        MList = ltml.LTDbList(lt2, PKey, 'SOLID_MATERIAL');
        MKey  = ltml.LTListNext(lt2, MList);
        ltml.LTDbSet(lt2, MKey, 'IndexOfRefraction', N_PLANAR);
    catch merr
        fprintf(['[VERIFY] planarization material assignment failed (%s). ' ...
                 'Assign n=%.3f to PlanarOvercoat manually in the GUI.\n'], ...
                 merr.message, N_PLANAR);
    end
end
% per-evaluation update (verified LTDbSet pattern, cf. Substrate Height/Y)
ltml.LTDbSet(lt2, PKey, 'Height', t_planar);
ltml.LTDbSet(lt2, PKey, 'Y', d_sub + t_planar/2);
%% !!! END VERIFY SECTION !!!

%% CPS
load('nk_JH33.mat');  load('Photopic_400_800.mat');
load('CIE_1931.mat'); load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';
wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98;  horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);

no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);  cos089=cosd(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL 150];
EML_position=4;  z0=12.5;  u_data_num=499;  max_u=3;

CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_sub_CPS=CPS_result.EQE_sub;

%% bottom reflectance
TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
Reflectance=(R_p_bottom+R_s_bottom)/2;

%% Coating
lt = ltloc.GetLTAPI(ID_LT);
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0', 'DATANAME: R_Bottom_',count, 'ABSORBING: YES', 'INDEX: 1.51', 'DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n', 'AOI',j, 0, Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);   % flush before LightTools reads the file

ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

%% wavelength loop
I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
sin089=sind(0:89);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);
% [Note] explicit index vector — never division-based sizing, safe for any
% (wavelength window, n) combination.
wv_list = 1:n:wavelength_num;
K = numel(wv_list);
I_air_1_2 = zeros(90, K);
Power_output = zeros(1, wavelength_num);
for kk = 1:K
    wv = wv_list(kk);
    fileID = fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1, 90, 0, 0, 360, 90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt, 'Lens_manager[1]','DISK_SOURCE');
    SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power', weight_factor(wv));
    SRList=ltml.LTDbList(lt, 'Lens_manager[1]','Spectral_region');
    SRKey=ltml.LTListAtPos(lt,SRList,2);
    ltml.LTDbSet(lt,SRKey,'Spectral_Definition', 'Monochromatic');
    ltml.LTDbSet(lt,SRKey,'Single_Wavelength', wv+wavelength_start-1);
    List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');
    Key=ltml.LTListAtPos(lt,List,1);
    ltml.LTDbSet(lt,Key,'LoadFileName',['C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\' sprintf('AI_temp.txt')]);

    ltml.LTBegin(lt);
    ltml.LTCmd(lt,'\V3D BeginAllSimulations');
    ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,3);
    for j=1:90
        I_air_1_JH(91-j,:)=ltml.LTDbGet(lt,Key,'CellValue_UI',1,91-j);
    end
    I_air_1_2(:,kk)=smooth(I_air_1_JH);
end

weight_factor_2  = zeros(K,1);
Power_output_2   = zeros(K,1);
EQE_sub_matrix_2 = zeros(K,1);
for k = 1:K
    idx = wv_list(k);
    weight_factor_2(k)  = weight_factor(idx);
    Power_output_2(k)   = Power_output(idx);
    EQE_sub_matrix_2(k) = CPS_result.EQE_sub_matrix(idx);
end
EQE_wv_matrix = Power_output_2 ./ weight_factor_2;
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;
EQE_total = sum(EQE_wv_matrix .* EQE_sub_matrix_2);

EQE_0_20=0; EQE_20_40=0; EQE_40_60=0; EQE_60_80=0;
sin_col = sin089(:);
for k = 1:K
    contrib_k = EQE_wv_matrix(k) * EQE_sub_matrix_2(k);
    W_theta = I_air_1_2(:,k) .* sin_col;  W_tot = sum(W_theta);
    EQE_0_20  = EQE_0_20  + contrib_k * sum(W_theta(1:20))  / W_tot;
    EQE_20_40 = EQE_20_40 + contrib_k * sum(W_theta(21:40)) / W_tot;
    EQE_40_60 = EQE_40_60 + contrib_k * sum(W_theta(41:60)) / W_tot;
    EQE_60_80 = EQE_60_80 + contrib_k * sum(W_theta(61:80)) / W_tot;
end

output = struct('EQE_0_20',EQE_0_20,'EQE_20_40',EQE_20_40, ...
    'EQE_40_60',EQE_40_60,'EQE_60_80',EQE_60_80,'EQE_total',EQE_total);

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d', count) ']" Delete= \Q']);
fclose('all');
end

%% ===== Spline constraints (identical to reference) =====
function TF = isValidPoints(X)
numRows = size(X,1);  numPts = 7;  TF = true(numRows,1);
for k = 1:numRows
    x = [0, X(k,1:5), 1];
    y = [1, X(k,6:10), 0];
    violates = false;
    for i = 1:numPts - 1
        for j = i + 2:numPts - 1
            if i == 1 && j == numPts - 1, continue; end
            if checkIntersection([x(i),y(i)],[x(i+1),y(i+1)],[x(j),y(j)],[x(j+1),y(j+1)])
                violates = true; break;
            end
        end
        if violates, break; end
    end
    if ~violates
        for i = 1:numPts - 2
            if isCollinear([x(i),y(i)],[x(i+1),y(i+1)],[x(i+2),y(i+2)])
                violates = true; break;
            end
        end
    end
    if ~violates
        minD = 0.05; maxD = 1.0;
        d = hypot(diff(x), diff(y));
        if any(d < minD | d > maxD), violates = true; end
    end
    if ~violates
        maxAng = 2 * pi / 3;
        for i = 2:numPts - 1
            v1 = [x(i),y(i)] - [x(i-1),y(i-1)];
            v2 = [x(i+1),y(i+1)] - [x(i),y(i)];
            ang = atan2(norm(cross([v1,0],[v2,0])), dot(v1,v2));
            if ang > maxAng, violates = true; break; end
        end
    end
    TF(k) = ~violates;
end
    function isCol = isCollinear(p1, p2, p3)
        area = 0.5 * det([p1 1; p2 1; p3 1]);
        isCol = abs(area) < 1e-5;
    end
    function intersects = checkIntersection(p1, p2, p3, p4)
        function o = orientation(p, q, r)
            o = (q(2)-p(2))*(r(1)-q(1)) - (q(1)-p(1))*(r(2)-q(2));
        end
        o1 = orientation(p1,p2,p3);  o2 = orientation(p1,p2,p4);
        o3 = orientation(p3,p4,p1);  o4 = orientation(p3,p4,p2);
        intersects = (o1*o2 < 0) && (o3*o4 < 0);
    end
end


function RenewLightTools()
global ID_LT ID_swept ltml ltloc lt
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
model_file_path_swept = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\SweptEntity.2.lts';
model_file_path_LT = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\Lens_size_effect_for_PSO_bump_modified_v1.1.lts';

fprintf('--- Restarting LightTools ---\n');
target_user = 'jhkim';
kill_cmd = sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user);
[~, ~] = system(kill_cmd);
pause(2);

cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_swept);
status = system(cmd); %#ok<NASGU>
try
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools restart failed. Check license / installation.');
end
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);
[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    if ~isempty(tokens)
        ID_swept = str2double(tokens{1}{1});
        fprintf('PID found for user %s: %d\n', target_user, ID_swept);
    else
        error('Process found but PID extraction failed.');
    end
else
    error('No LightTools process found for user %s.', target_user);
end
cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_LT);
status = system(cmd); %#ok<NASGU>
[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    if ~isempty(tokens)
        ID_LT = str2double(tokens{3}{1});
        fprintf('PID found for user %s: %d\n', target_user, ID_LT);
    else
        error('Process found but PID extraction failed.');
    end
else
    error('No LightTools process found for user %s.', target_user);
end
pause(5);
end
