% ============================================================
%  BO_asym_v2.m  -- 비대칭 freeform 렌즈 + detuned 마이크로캐비티 공동최적화
%
%  [목적함수] phi 40도 창(위치 자동검출) 안, 특정 theta 밴드로 몰리는 EQE(파워)의
%    "절대값" 최대화. (분율 아님. 추출량 x 방위각 집중을 동시에 보상.)
%
%  [변수] 형상 3 DOF + 4개(dETL, dHTL, stretchZ, Decenter) = 총 7개
%    - xc       : 정점 x-오프셋(비대칭 강도).      0=대칭 ~ 0.5=강한 앞쏠림
%    - s_steep  : 앞(+x) 플랭크 압축(가파름).       1=대칭 ~ 2=가파름(draft<=~65도)
%    - s_gentle : 뒤(-x) 플랭크 확장(완만).         0.5=긴 완만램프 ~ 1=대칭
%    - dETL,dHTL: OLED 마이크로캐비티 두께 (나노)
%    - stretchZ : 렌즈 전체 높이(텍스처 z-스케일). 형상은 단위높이로 생성됨
%    - Decenter : 발광원-렌즈 편심 (PLANAR_REFERENCE_SURFACE.X = 15+Decenter, 기존 배선)
%
%  [geometry] SweptEntity(회전대칭 스플라인) 제거. asym_dome 을 MATLAB 이 직접
%    .ent(X,Y,Z 격자)로 생성/저장(generate_asym_ent) 후, 배열모델 텍스처의
%    LibraryElement Filename 으로 지정. -> ID_swept(스웹 인스턴스) 불필요.
%
%  ================= @@CONFIRM : 내가 내린 판단 (틀리면 알려주세요) =================
%   (C1) 형상 = 제작가능 비대칭 raised-cosine 돔(asym_dome), 형상 DOF 3개.
%        더 자유롭게 하려면 여기 파라미터화만 교체(RBF 등). BO/objFcn 나머지 재사용.
%   (C2) 높이는 h 대신 stretchZ 로 스케일(형상 .ent 는 단위높이). h/stretchZ 중복
%        자유도 제거 목적. h 를 독립 변수로 원하면 varNames 에 추가만 하면 됨.
%   (C3) 목표 theta 밴드 = [40,60], phi 창폭 = 40도. 응용 바뀌면 TH_LO/TH_HI/PHI_W
%        세 줄만 수정.
%   (C4) 원점 템플릿 freeform_template_v2.ent 가 BASE 경로에 있어야 함(레포에 있음).
%   (C5) ID_swept/LT1(스웹 모델)은 이제 미사용. 속도 위해 단일 인스턴스 런처로
%        바꿔도 되나, 기존 RenewLightTools_3(2개) 그대로 둬도 동작함(1개는 낭비).
%   (C6) far-field 2D mesh = INTENSITY_MESH 위치 3, [90 theta x 36 phi],
%        phi 중심 -175:10:175 로 가정(기존 export 형식과 동일). 다르면 read 부 수정.
% ================================================================================
clear;

%% ===== LightTools 연결 =====
global ID_swept ID_LT ltml ltloc count ray_nums_current LT1 LT2 r_pat
global FF_TEMPLATE FF_BASE
tic;
FF_BASE = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
LT1 = [FF_BASE 'SweptEntity_asym.1.lts'];      % (C5) 이제 미사용(런처 호환 위해 유지)
LT2 = [FF_BASE 'assymetric_test.1.lts'];       % 배열 모델(실사용)
FF_TEMPLATE = [FF_BASE 'freeform_template_v2.ent'];   % (C4) 원점 템플릿

RenewLightTools_3(LT1, LT2);
toc;
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 고정 설정 =====
global RESTART_INTERVAL eval_count TH_LO TH_HI PHI_W
RESTART_INTERVAL = 20;    % 시뮬 N회마다 LightTools 재시작
eval_count       = 0;
r_pat            = 2;     % 렌즈 피치(텍스처 zone). 예전엔 sweep, 이제 단일 고정값 @@CONFIRM
TH_LO = 40;  TH_HI = 60;  PHI_W = 40;   % (C3) 목표 theta 밴드 / phi 창폭 [deg]

RAY_SEARCH   = 20000;    % 탐색용 ray (빠름)
RAY_FINAL    = 100000;   % 최종 검증용 ray (정밀)
N_FINAL_REP  = 3;        % 최종 best 반복 평가 -> mean±std
POLISH_EVALS = 15;       % patternsearch 정련 예산
EXPLORATION_RATIO = 0.7;
INIT_EVAL    = 60;       % 초기 BO 예산(시드 포함)
CONV_BLOCK   = 10;  CONV_TOL = 0.002;  CONV_PATIENCE = 2;  MAX_EVAL = 150;

%% ===== 최적화 변수 (형상 3 + 4 = 7) =====
varNames = {'xc','s_steep','s_gentle','dETL','dHTL','stretchZ','Decenter'};
lb = [0.00, 1.0, 0.50,  10,  10, 0.1, 0.0];
ub = [0.50, 2.0, 1.00, 150, 150, 3.0, 7.5];
optVars = optimizableVariable.empty(0, numel(lb));
for i = 1:numel(lb)
    optVars(i) = optimizableVariable(varNames{i}, [lb(i), ub(i)]);
end
DOF = numel(lb);
fprintf('DOF = %d (형상 3: xc,s_steep,s_gentle + dETL,dHTL,stretchZ,Decenter)\n', DOF);

%% ===== 단일 BO 실행 =====
ray_nums_current = RAY_SEARCH;
initX = array2table(lb + rand(20, DOF).*(ub - lb), 'VariableNames', varNames);

fprintf('\n######## 비대칭 BO 시작: target θ∈[%d,%d], φ창=%d° (φ중심 자동검출) ########\n', ...
    TH_LO, TH_HI, PHI_W);

results = bayesopt(@bo_objective, optVars, ...
    'MaxObjectiveEvaluations', INIT_EVAL, ...
    'IsObjectiveDeterministic', false, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'ExplorationRatio', EXPLORATION_RATIO, ...
    'Verbose', 1, 'PlotFcn', {}, 'InitialX', initX);

% 수렴 판정: 개선 멈출 때까지 CONV_BLOCK 씩 추가
noImp = 0;
while noImp < CONV_PATIENCE && results.NumObjectiveEvaluations < MAX_EVAL
    prevEst = -results.MinEstimatedObjective;
    addEval = min(CONV_BLOCK, MAX_EVAL - results.NumObjectiveEvaluations);
    results = resume(results, 'MaxObjectiveEvaluations', addEval);
    newEst  = -results.MinEstimatedObjective;
    relImp  = (newEst - prevEst) / max(abs(newEst), eps);
    if relImp < CONV_TOL, noImp = noImp + 1; else, noImp = 0; end
    fprintf('[Converge] evals=%3d | bestEst EQE_region=%.5g | relImp=%+.4f | noImp %d/%d\n', ...
        results.NumObjectiveEvaluations, newEst, relImp, noImp, CONV_PATIENCE);
end

%% ===== 국소 정련 (patternsearch) =====
x0 = table2array(results.XAtMinEstimatedObjective);
fprintf('--- Local polish (patternsearch, %d evals) ---\n', POLISH_EVALS);
psOpts = optimoptions('patternsearch', 'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.05, 'MeshTolerance', 1e-3, 'Cache', 'on', 'Display', 'iter');
try
    xPol = patternsearch(@polish_objective, x0, [],[],[],[], lb, ub, [], psOpts);
catch perr
    fprintf('[Warn] patternsearch 실패(%s). BO 결과만 사용.\n', perr.message);  xPol = x0;
end

%% ===== 최종 고정밀 검증 =====
if isequal(xPol, x0), candX = {x0}; else, candX = {x0, xPol}; end
ray_nums_current = RAY_FINAL;
candMean = -inf(1, numel(candX));  candStd = zeros(1, numel(candX));
for c = 1:numel(candX)
    e = nan(1, N_FINAL_REP);
    for rrep = 1:N_FINAL_REP, e(rrep) = simulate_metric(candX{c}); end
    candMean(c) = mean(e, 'omitnan');  candStd(c) = std(e, 'omitnan');
    fprintf('  candidate %d : EQE_region = %.5g ± %.2g (%d rays)\n', ...
        c, candMean(c), candStd(c), RAY_FINAL);
end
[bestEQE, ci] = max(candMean);
if ~isfinite(bestEQE), bestEQE = -results.MinEstimatedObjective; ci = 1; end
bestX = array2table(candX{ci}, 'VariableNames', varNames);

%% ===== 결과 리포트 + 저장 (best 형상 분해 재평가) =====
ray_nums_current = RAY_FINAL;
bd = objFcn_regionPower(table2array(bestX));
save('BO_asym_result.mat', 'bestX', 'bestEQE', 'results', 'bd', ...
    'TH_LO', 'TH_HI', 'PHI_W', 'varNames', 'lb', 'ub', 'r_pat');
fprintf('\n######## Done ########\n');
fprintf('  목적함수 EQE_region(절대) = %.5g ± %.2g\n', bestEQE, candStd(ci));
fprintf('  EQE_total = %.5g | 검출 φ중심 = %+.0f° | φ대비비(전/후) = %.2f\n', ...
    bd.EQE_total, bd.phiC, bd.contrast);
disp(bestX);


%% =====================================================================
%%  1회 평가 래퍼 (주기 재시작 + 크래시 -> NaN)
%% =====================================================================
function m = simulate_metric(pt)
global ID_LT ltml ltloc eval_count RESTART_INTERVAL LT1 LT2
eval_count = eval_count + 1;
if mod(eval_count, RESTART_INTERVAL) == 0
    fprintf('\n[Refresh] 시뮬 %d회. LightTools 재시작...\n', eval_count);
    RenewLightTools_3(LT1, LT2);
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);  pause(1);
end
try
    m = objFcn_regionPower(pt).EQE_region;    % 목적함수 = 창내 절대 EQE
    if ~isfinite(m), m = NaN; end
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    m = NaN;
    RenewLightTools_3(LT1, LT2);
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end

function obj = bo_objective(Xtbl)
obj = -simulate_metric(table2array(Xtbl));    % 최대화 -> -값 최소화
end

function f = polish_objective(x)
e = simulate_metric(x);
if isnan(e), e = 0; end
f = -e;
end


%% =====================================================================
%%  Objective: asym .ent 직접생성 + 나노 CPS + phi 자동검출 창내 절대 EQE
%% =====================================================================
function output = objFcn_regionPower(point)
global ID_LT ltml ltloc count ray_nums_current r_pat FF_TEMPLATE FF_BASE
global TH_LO TH_HI PHI_W

lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

% --- 변수 언팩 ---
xc       = point(1);  s_steep = point(2);  s_gentle = point(3);
dETL     = point(4);  dHTL    = point(5);  stretchZ = point(6);  Decenter = point(7);

% --- 고정 설정 ---
d_sub = 1.3;  r_OLED = 1;  Lensheight = 0.01;
wavelength_start = 580;  wavelength_end = 590;  n = 10;
if isempty(ray_nums_current), ray_nums = 50000; else, ray_nums = ray_nums_current; end

List=ltml.LTDbList(lt,'lens_manager[1]','SIMULATIONS');
Key=ltml.LTListByName(lt,List,'ForwardAll');
ltml.LTDbSet(lt,Key,'MaxProgress',ray_nums);
List=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
Key=ltml.LTListByName(lt,List,'Substrate');
ltml.LTDbSet(lt,Key,'Height',d_sub);  ltml.LTDbSet(lt,Key,'Y',d_sub/2);
SRList=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
SRKey=ltml.LTListAtPos(lt,SRList,2);
ltml.LTDbSet(lt,SRKey,'Y',d_sub+Lensheight/2);
List=ltml.LTDbList(lt,'lens_manager[1]','TEXTURE_ZONE_EXTENT');
Key=ltml.LTListByName(lt,List,'zone');
ltml.LTDbSet(lt,Key,'Geometry_1',r_pat);  ltml.LTDbSet(lt,Key,'Geometry_2',r_pat);
% 편심(기존 배선 유지): 레퍼런스 평면 X 이동
List=ltml.LTDbList(lt,'lens_manager[1]','PLANAR_REFERENCE_SURFACE');
Key=ltml.LTListByName(lt,List,'ReferenceSurface');
ltml.LTDbSet(lt,Key,'X',15+Decenter);
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% --- (신규) asym .ent 직접 생성 + 텍스처에 물리기 ---
rng('shuffle');
charSet=['a':'z' 'A':'Z' '0':'9'];
tag = charSet(randi(numel(charSet),1,10));
entPath = [FF_BASE 'asym_' tag '.1.ent'];
ok = generate_asym_ent(xc, s_steep, s_gentle, FF_TEMPLATE, entPath);
if ~ok || ~exist(entPath,'file')
    output = fail_output();  return;
end
List = ltml.LTDbList(lt,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
Key  = ltml.LTListByName(lt,List,'LibraryElement');
ltml.LTDbSet(lt,Key,'Filename', entPath);
List = ltml.LTDbList(lt,'LENS_MANAGER[1]','TEXTURE_PARAMETER');
Key  = ltml.LTListByName(lt,List,'StretchZ');
ltml.LTDbSet(lt,Key,'Value', stretchZ);

% --- 나노 CPS + 하단 반사율 코팅 (기존과 동일) ---
load('nk_JH33.mat');  load('Photopic_400_800.mat');  load('CIE_1931.mat');  load('R_pd.mat');
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

TMF_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
Reflectance=(abs(TMF_p.r_p).^2 + abs(TMF_s.r_s).^2)/2;

fileID=fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa',count),'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0','DATANAME: R_Bottom_',count,'ABSORBING: YES','INDEX: 1.51','DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n','AOI',j,0,Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa',count) '"']);
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');  Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');  Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d',count));

I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);

% --- 파장 루프: 시뮬 + 2D far-field(90x36) 누적 ---
nLat=90; nLong=36;
K=(wavelength_num-1)/n+1;
Power_output=zeros(1,wavelength_num);
Igrids=cell(1,K);
for wv=1:n:wavelength_num
    fileID=fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1,90,0,0,360,90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','DISK_SOURCE');  SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power',weight_factor(wv));
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','Spectral_region');  SRKey=ltml.LTListAtPos(lt,SRList,2);
    ltml.LTDbSet(lt,SRKey,'Spectral_Definition','Monochromatic');
    ltml.LTDbSet(lt,SRKey,'Single_Wavelength',wv+wavelength_start-1);
    List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');  Key=ltml.LTListAtPos(lt,List,1);
    ltml.LTDbSet(lt,Key,'LoadFileName','C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt');

    ltml.LTBegin(lt);  ltml.LTCmd(lt,'\V3D BeginAllSimulations');  ltml.LTEnd(lt);

    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');  Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');
    % 2D 방향성 mesh(위치 3) -> Igrid(nLat x nLong)  (기존 인덱싱과 동일)
    Key=ltml.LTListAtPos(lt,List,3);
    Ig=zeros(nLat,nLong);
    for rr=1:nLat
        for kk=1:nLong
            v=ltml.LTDbGet(lt,Key,'CellValue_UI',kk,rr);
            if isempty(v)||~isfinite(v), v=0; end
            Ig(rr,kk)=v;
        end
    end
    Igrids{(wv+n-1)/n}=Ig;
end

% --- 파장 가중 -> EQE_total + 2D 누적 Wacc ---
weight_factor_2=zeros(K,1);  Power_output_2=zeros(K,1);  EQE_sub_matrix_2=zeros(K,1);
for k=1:K
    idx=n*(k-1)+1;
    weight_factor_2(k)=weight_factor(idx);
    Power_output_2(k)=Power_output(idx);
    EQE_sub_matrix_2(k)=CPS_result.EQE_sub_matrix(idx);
end
EQE_wv_matrix=Power_output_2./weight_factor_2;
EQE_sub_matrix_2=EQE_sub_matrix_2/sum(EQE_sub_matrix_2)*EQE_sub_CPS;
contrib=EQE_wv_matrix.*EQE_sub_matrix_2;          % 파장별 EQE 기여
EQE_total=sum(contrib);

thC = ( (1:nLat) - 0.5 ) * (90/nLat);             % theta 중심 [deg] 0.5..89.5
phC = -180 + ( (1:nLong) - 0.5 ) * (360/nLong);   % phi 중심 [deg] -175..175
sint = sind(thC(:));
Wacc = zeros(nLat, nLong);
for k=1:K
    Wk = Igrids{k} .* sint;   sk = sum(Wk(:));
    if sk>0, Wacc = Wacc + contrib(k)*(Wk/sk); end
end
% Wacc 총합 = EQE_total, 창내 합 = 창으로 나가는 절대 EQE

% --- phi 창 자동검출 -> 창내 절대 EQE ---
[PWin, phiC, contrast] = detect_phi_window(Wacc, thC, phC, TH_LO, TH_HI, PHI_W);

output = struct('EQE_region',PWin, 'EQE_total',EQE_total, 'phiC',phiC, ...
    'contrast',contrast, 'thBand',[TH_LO TH_HI], 'phiWidth',PHI_W);
fprintf('[obj] EQE_region=%.5g | EQE_total=%.5g | φc=%+.0f° | contrast=%.2f | xc=%.2f ss=%.2f sg=%.2f dec=%.2f\n', ...
    PWin, EQE_total, phiC, contrast, xc, s_steep, s_gentle, Decenter);

% 코팅 정리
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');  Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');  Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d',count) ']" Delete= \Q']);
fclose('all');
end

function output = fail_output()
output = struct('EQE_region',0,'EQE_total',0,'phiC',NaN,'contrast',0,'thBand',[NaN NaN],'phiWidth',NaN);
end


%% =====================================================================
%%  형상 생성: 비대칭 raised-cosine 돔 -> .ent 직접 쓰기 (단위높이)
%% =====================================================================
function ok = generate_asym_ent(xc, s_steep, s_gentle, templatePath, outPath)
% asym_dome: x축을 앞(+x)/뒤(-x) 비대칭 스케일한 raised-cosine 돔. 단위높이(정점=1),
% 실제 높이는 배열 텍스처 StretchZ 가 스케일. rim window 로 테두리 0(타일링) 강제.
ok = false;
try
    Ra=1.2139; Rap=1.0; n=141; tbase=0.30;
    g=linspace(-Ra,Ra,n); [X,Y]=meshgrid(g,g); r=hypot(X,Y);
    xr=X-xc;  sx=s_gentle*ones(size(X));  sx(xr>0)=s_steep;
    Xe=xr.*sx;  rho=sqrt(Xe.^2+Y.^2)/Rap;
    H=0.5*(1+cos(pi*min(max(rho,0),1)));
    W=ones(size(r)); rw=0.85*Rap;
    z=(r-rw)/(Rap-rw); m=(r>=rw)&(r<Rap);
    W(m)=1-(3*z(m).^2-2*z(m).^3); W(r>=Rap)=0;
    H=max(H.*W,0);

    Z=H; Xv=X(:); Yv=Y(:); Zv=Z(:); N=n*n;
    tpl=fileread(templatePath);
    tok=regexp(tpl,'ORAStartData;([\s\S]*?)ORAEndData;','tokenExtents');
    s0=tok{1}(1); e0=tok{1}(2);
    buf=sprintf('0 1 %d %d 0 0 %d 0 0 0',n,n,N);
    for i=1:N, buf=[buf sprintf(' %.17g %.17g %.17g',Xv(i),Yv(i),Zv(i))]; end %#ok<AGROW>
    buf=[buf ' 0 0 4 CartesianMapper 1 0 0 0 0'];
    newtxt=[tpl(1:s0-1) char(10) buf char(10) tpl(e0+1:end)];
    newtxt=regexprep(newtxt, ...
        '(CSGLensSurfacePrimitive_1[\s\S]*?setPosition:  \{ 0\. 0\. )[-0-9.eE]+(  \} ;)', ...
        ['$1' num2str(-tbase,'%g') '$2'],'once');
    newtxt=regexprep(newtxt,'restoreSmoothResample: "Yes"','restoreSmoothResample: "No"','once');
    fid=fopen(outPath,'w'); fwrite(fid,newtxt); fclose(fid);
    ok = true;
catch me
    fprintf('[Geom] .ent 생성 실패: %s\n', me.message);
end
end


%% =====================================================================
%%  phi 창 자동검출: theta 밴드 고정, phi 창(폭 phiWidth)을 wrap 슬라이딩
%% =====================================================================
function [PWin, phiC, contrast] = detect_phi_window(Wacc, thC, phC, thLo, thHi, phiWidth)
tm  = (thC>=thLo) & (thC<=thHi);
band = sum(Wacc(tm,:), 1);          % 1 x nLong : theta 밴드 내 phi 분포
nL  = numel(phC);  half = phiWidth/2;
winP = zeros(1,nL);
for c=1:nL
    d = abs(mod(phC - phC(c) + 180, 360) - 180);
    winP(c) = sum(band(d <= half));
end
[PWin, ic] = max(winP);  phiC = phC(ic);
dOpp = abs(mod(phC - (phiC+180) + 180, 360) - 180);
POpp = sum(band(dOpp <= half));
contrast = PWin / max(POpp, eps);   % 전/후 phi 대비비
end


%% =====================================================================
%%  RenewLightTools_3 (2 인스턴스, start /min + 폴링)  -- 기존 그대로
%% =====================================================================
function RenewLightTools_3(path1, path2)
global ID_LT ID_swept ltml ltloc
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
LT1 = path1;  LT2 = path2;
fprintf('--- Restarting LightTools ---\n');
target_user = 'jhkim';
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);
system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
t0=tic;
while toc(t0) < 10
    [~, cmdout] = system(find_cmd);
    if ~contains(cmdout,'lt.exe'), break; end
    pause(0.3);
end
system(sprintf('start /min "" "%s" "%s"', lt_exe_path, LT1));
try
    ltml  = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패. 라이선스/설치 확인.');
end
t1=tic; found1=false;
while toc(t1) < 20
    [status, cmdout] = system(find_cmd);
    if status==0 && contains(cmdout,'lt.exe'), found1=true; break; end
    pause(0.3);
end
if ~found1, error('lt.exe(swept) 탐색 실패'); end
tokens = regexp(cmdout, '"(\d+)"', 'tokens');
ID_swept = str2double(tokens{1}{1});
fprintf('PID(swept)=%d\n', ID_swept);
system(sprintf('start /min "" "%s" "%s"', lt_exe_path, LT2));
t2=tic; found2=false;
while toc(t2) < 20
    [status, cmdout] = system(find_cmd);
    if status==0 && contains(cmdout,'lt.exe')
        tokens = regexp(cmdout, '"(\d+)"', 'tokens');
        if numel(tokens)>=3, found2=true; break; end
    end
    pause(0.3);
end
if ~found2, error('lt.exe(LT) 탐색 실패'); end
ID_LT = str2double(tokens{3}{1});
fprintf('PID(LT)=%d\n', ID_LT);
for whichID = [ID_swept, ID_LT]
    tR=tic; ready=false;
    while toc(tR) < 20
        try
            lt = ltloc.GetLTAPI(whichID);
            ltml.LTCmd(lt, 'Message "Check Connection"');
            ready=true; break;
        catch
            pause(0.5);
        end
    end
    if ~ready
        fprintf('[경고] PID=%d COM 준비 확인 실패.\n', whichID);
    end
end
end
