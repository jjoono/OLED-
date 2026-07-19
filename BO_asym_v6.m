% ============================================================
%  BO_asym_v6.m  -- v5 + 기판두께(d_sub) 변수화
%
%  [v5 대비 변경] d_sub 를 최적화 변수로 추가, 범위 [0.5, 3.0] mm.
%    per-lenslet 소스 각크기 atan(r_OLED/d_sub) 를 옵티마이저가 직접 조절
%    -> 패턴크기(rx,ry)와의 비율이 스티어링 기하 전체를 결정 (시너지 최대).
%    [주의] 실제 데모에서 기판두께 조절엔 한계가 있음 - v5(고정 1.3) 결과와
%    비교해 "두께 자유가 주는 상한"을 보는 용도.
%  ---- 이하 v5 헤더 ----
%
%  [v4 대비 변경]
%   - MLA 패턴 크기 변수화 + 비등방: rx_pat, ry_pat 각각 [2,25]
%       (Geometry_1=rx, Geometry_2=ry; 막대형/렌티큘러형 셀 가능)
%   - 2D 편심: decX_frac(0~0.5), decY_frac(-0.5~0.5), 셀 대비 비율
%       실제 편심 = frac * 해당축 패턴크기 -> 항상 셀 안. 임의 φ 방향 조준.
%   - ReferenceSurface 기준 X=15 는 의도된 고정값(사용자 확인). r_pat=15 와
%     같았던 것은 우연 - 연동하지 않음.
%  ---- 이하 v4 헤더 ----
%
%  [v3 대비 변경]
%   (1) 실시간 진행표시: surrogateopt OutputFcn(@live_progress) -> 매 평가 best-so-far
%       곡선 갱신 + 주기적 요약(단계/φc/contrast/draft/경과/ETA) 콘솔 출력.
%   (2) 형상 자유도 대폭 확대: raised-cosine 베이스 + Nrbf RBF 요철 + Ncut 절단면
%       (min(H, 기울어진 평면) -> 날카로운 비대칭 단면). 기본 Nrbf=10, Ncut=4
%       -> 형상 4*10 + 3*4 = 52 DOF. "돔에서 거기서 거기" 탈피.
%   (3) 소자구조 변경: ITO(l_ITO_SNU_temp) -> Ag(l_Ag_McPeak),
%       thickness = [100 dETL 25 10 dHTL dAg]. dAg 를 최적화 변수 추가(범위 0~50).
%
%  [변수] 형상(4*Nrbf+3*Ncut) + 5개(dETL,dHTL,dAg,stretchZ,Decenter)
%  [목적함수] phi 40도 창(자동검출) 안, theta[40,60] 로 몰리는 EQE(파워) 절대값.
%
%  ================= @@CONFIRM =================
%   (C1) 형상 = raised-cosine + RBF + 절단면. Nrbf/Ncut 로 DOF 조절.
%   (C2) 높이는 stretchZ 변수만(형상 단위높이).
%   (C3) theta[40,60], phi 창폭 40.   (C5) r_pat=15.   (C6) mesh 45x36.
%   (C8) 소자: TAPC 다음층이 Ag(l_Ag_McPeak), 두께 dAg. nk_JH33.mat 의 material
%        구조체에 l_Ag_McPeak 필드가 있어야 함(사용자 확인).
% ============================================================
clear;

%% ===== LightTools 연결 (단일 인스턴스) =====
global ID_LT ltml ltloc count ray_nums_current LT2
global FF_TEMPLATE FF_BASE RESTART_INTERVAL eval_count TH_LO TH_HI PHI_W NRBF NCUT
global LAST_METRICS OPT_HIST SUMMARY_EVERY MAXEVAL
tic;
FF_BASE     = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
LT2         = [FF_BASE 'assymetric_test_v2.1.lts'];
FF_TEMPLATE = [FF_BASE 'freeform_template_v2.ent'];

RenewLightTools_single(LT2);
toc;
count = 1;
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 고정 설정 =====
RESTART_INTERVAL = 20;
eval_count = 0;
TH_LO = 40;  TH_HI = 60;  PHI_W = 40;
NRBF = 10;   NCUT = 4;                 % (C1) 형상 DOF = 4*NRBF + 3*NCUT = 52
LAST_METRICS = struct('phiC',NaN,'contrast',NaN,'draft',NaN);  % live 표시용
SUMMARY_EVERY = 10;                    % N 평가마다 요약 출력

RAY_SEARCH  = 20000;   RAY_FINAL = 100000;   N_FINAL_REP = 3;
MAXEVAL     = 600;     % 60변수엔 빠듯(권장 1200~3000). 아래 RESUME 로 이어서 연장 가능.
% [횟수 늘리는 법]
%  (1) 이어달리기: 1차 run 종료 후 MAXEVAL 을 더 크게(예: 1500 - "누적 총횟수"임)
%      바꾸고 RESUME_FROM_CKPT=true 로 재실행 -> checkpoint 에서 이어서 최적화.
%      (live 카운터가 이전 횟수부터 이어짐. 목적함수가 이 스크립트 안에 있으므로
%       반드시 "같은 스크립트를 다시 실행"하는 방식으로 재개할 것.)
%  (2) 평가당 비용 절감(횟수 키우려면 사실상 필수):
%      n=10 -> 30 (파장 31->11개, ~3배 빠름) + RAY_SEARCH 20000 -> 10000 (~2배)
%      => ~6배. 스펙트럼/ray 정밀도는 최종 고정밀 검증(RAY_FINAL, n 원복)에서 보상.
RESUME_FROM_CKPT = false;   % true 면 v6_ckpt 에서 이어서 (MAXEVAL 은 누적 총횟수)

%% ===== 변수 (형상 4*NRBF + 3*NCUT + 5) =====
shapeNames = {};  lbS = [];  ubS = [];
for i = 1:NRBF   % RBF: x0,y0,amp,sig
    shapeNames = [shapeNames, {sprintf('x0_%d',i),sprintf('y0_%d',i),sprintf('amp_%d',i),sprintf('sig_%d',i)}]; %#ok<AGROW>
    lbS = [lbS, -0.8, -0.8, -0.60, 0.12];  ubS = [ubS, 0.8, 0.8, 0.70, 0.60]; %#ok<AGROW>
end
for i = 1:NCUT   % 절단면: z0(높이), m(기울기), phi0(방향)
    shapeNames = [shapeNames, {sprintf('cz_%d',i),sprintf('cm_%d',i),sprintf('cphi_%d',i)}]; %#ok<AGROW>
    lbS = [lbS, 0.20, 0.0, 0.0];  ubS = [ubS, 1.40, 1.6, 2*pi]; %#ok<AGROW>
end
varNames = [shapeNames, {'dETL','dHTL','dAg','stretchZ','decX_frac','decY_frac','rx_pat','ry_pat','d_sub'}];
lb = [lbS,  10,  10,  0,  0.1, 0.0, -0.5,  2,  2, 0.5];
ub = [ubS, 150, 150, 50,  3.0, 0.5,  0.5, 25, 25, 3.0];
NV = numel(lb);   S_DIM = 4*NRBF + 3*NCUT;
fprintf('총 변수 = %d (형상 %d | + dETL,dHTL,dAg,stretchZ,decX,decY,rx_pat,ry_pat,d_sub)\n', ...
    NV, S_DIM);

%% ===== surrogateopt 실행 (실시간 진행표시) =====
ray_nums_current = RAY_SEARCH;
optsSurr = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', MAXEVAL, ...
    'UseParallel', false, ...
    'MinSurrogatePoints', max(2*NV+1, 60), ...
    'OutputFcn', @live_progress, ...      % <-- 실시간 진행
    'PlotFcn', [], ...
    'CheckpointFile', 'BO_asym_v6_ckpt.mat', ...
    'Display', 'iter');

fprintf('\n######## surrogateopt 시작: %d DOF, target θ∈[%d,%d], φ창=%d° ########\n', ...
    NV, TH_LO, TH_HI, PHI_W);
ckptFile = 'BO_asym_v6_ckpt.mat';
if RESUME_FROM_CKPT && isfile(ckptFile)
    fprintf('[Resume] %s 에서 이어서 최적화 (MAXEVAL=%d 은 누적 총횟수)\n', ckptFile, MAXEVAL);
    [xBest, fBest, exitflag, outS] = surrogateopt(ckptFile, optsSurr); %#ok<ASGLU>
else
    [xBest, fBest, exitflag, outS] = surrogateopt(@surr_objective, lb, ub, optsSurr); %#ok<ASGLU>
end

%% ===== 최종 고정밀 검증 =====
ray_nums_current = RAY_FINAL;
e = nan(1, N_FINAL_REP);
for rrep = 1:N_FINAL_REP, e(rrep) = simulate_metric(xBest); end
bestEQE = mean(e,'omitnan');  bestStd = std(e,'omitnan');
bd = objFcn_regionPower(xBest);
save('BO_asym_v6_result.mat', 'xBest','bestEQE','bestStd','bd','varNames','lb','ub', ...
    'NRBF','NCUT','TH_LO','TH_HI','PHI_W','outS','OPT_HIST');
fprintf('\n######## Done ########\n');
fprintf('  EQE_region(절대) = %.5g ± %.2g (surrogate fBest=%.5g)\n', bestEQE, bestStd, -fBest);
fprintf('  EQE_total=%.5g | φ중심=%+.0f° | φ대비비=%.2f | maxDraft=%.1f°\n', ...
    bd.EQE_total, bd.phiC, bd.contrast, bd.maxDraft);


%% =====================================================================
%%  실시간 진행표시 OutputFcn
%% =====================================================================
function stop = live_progress(~, optimValues, state)
stop = false;
persistent t0 fcArr eqeArr bestArr bestSoFar figH ax
global LAST_METRICS SUMMARY_EVERY MAXEVAL OPT_HIST
switch state
    case 'init'
        t0=tic; fcArr=[]; eqeArr=[]; bestArr=[]; bestSoFar=-inf;
        figH=figure('Name','BO asym 실시간 진행','Color','w');
        ax=axes(figH); hold(ax,'on'); grid(ax,'on');
        xlabel(ax,'평가 횟수'); ylabel(ax,'EQE\_region (창내 절대)');
        title(ax,'파랑=매 평가, 빨강=best-so-far');
    case 'done'
        if ~isempty(bestSoFar)
            fprintf('\n===== [완료] 총 %d 평가 | best EQE_region=%.4g | %.1f분 =====\n', ...
                optimValues.funccount, bestSoFar, toc(t0)/60);
        end
    otherwise  % 'iter'
        f = optimValues.fval;
        if isfinite(f), eqe = -f; else, eqe = NaN; end
        n = optimValues.funccount;
        fcArr(end+1)=n;  eqeArr(end+1)=eqe;
        if isfinite(eqe) && eqe>bestSoFar, bestSoFar=eqe; end
        bestArr(end+1)=bestSoFar;
        if ~isempty(figH) && isvalid(figH)
            cla(ax);
            plot(ax, fcArr, eqeArr, 'b.', 'MarkerSize',8);
            plot(ax, fcArr, bestArr, 'r-', 'LineWidth',1.6);
            drawnow limitrate;
        end
        if mod(numel(fcArr), SUMMARY_EVERY)==0
            el=toc(t0); rate=numel(fcArr)/max(el,eps); eta=(MAXEVAL-n)/max(rate,eps);
            ph=''; try, ph=char(string(optimValues.currentFlag)); catch, end
            lm=LAST_METRICS;
            fprintf(['\n===== [진행] %d/%d 평가 | best EQE_region=%.4g | ' ...
                'φc=%+.0f° contrast=%.2f draft=%.0f° | 단계=%s | %.1f분 경과, ETA~%.1f분 =====\n\n'], ...
                n, MAXEVAL, bestSoFar, lm.phiC, lm.contrast, lm.draft, ph, el/60, eta/60);
        end
        OPT_HIST=struct('fc',fcArr,'eqe',eqeArr,'best',bestArr);
end
end

%% ===== surrogateopt 목적함수 (최소화) =====
function f = surr_objective(x)
m = simulate_metric(x);
if isnan(m), f = NaN; else, f = -m; end
end

%% ===== 1회 평가 래퍼 =====
function m = simulate_metric(pt)
global ID_LT ltml ltloc eval_count RESTART_INTERVAL LT2
eval_count = eval_count + 1;
if mod(eval_count, RESTART_INTERVAL) == 0
    fprintf('\n[Refresh] 시뮬 %d회. LightTools 재시작...\n', eval_count);
    RenewLightTools_single(LT2);
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);  pause(1);
end
try
    out = objFcn_regionPower(pt);
    m = out.EQE_region;
    if ~isfinite(m) || out.loadFailed, m = NaN; end
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    m = NaN;
    RenewLightTools_single(LT2);
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end


%% =====================================================================
%%  Objective
%% =====================================================================
function output = objFcn_regionPower(point)
global ID_LT ltml ltloc count ray_nums_current FF_TEMPLATE FF_BASE
global TH_LO TH_HI PHI_W NRBF NCUT LAST_METRICS

lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

S = 4*NRBF + 3*NCUT;
pShape   = point(1:S);
dETL     = point(S+1);  dHTL = point(S+2);  dAg = point(S+3);
stretchZ = point(S+4);
decX_frac= point(S+5);  decY_frac = point(S+6);
rx_pat   = point(S+7);  ry_pat    = point(S+8);
d_sub    = point(S+9);            % 기판두께 변수 (v6)
DecenterX = decX_frac * rx_pat;   % 편심 = 셀 대비 비율 -> 항상 셀 안
DecenterY = decY_frac * ry_pat;

r_OLED=1; Lensheight=0.01;   % d_sub 는 변수(위에서 언팩)
wavelength_start=450; wavelength_end=750; n=10;
if isempty(ray_nums_current), ray_nums=50000; else, ray_nums=ray_nums_current; end

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
ltml.LTDbSet(lt,Key,'Geometry_1',rx_pat);  ltml.LTDbSet(lt,Key,'Geometry_2',ry_pat);
List=ltml.LTDbList(lt,'lens_manager[1]','PLANAR_REFERENCE_SURFACE');
Key=ltml.LTListByName(lt,List,'ReferenceSurface');
ltml.LTDbSet(lt,Key,'X',15+DecenterX);   % 기준 X=15 는 의도된 고정값(사용자 확인)
ltml.LTDbSet(lt,Key,'Y',15+DecenterY);   % @@CONFIRM 모델의 ReferenceSurface Y 기본값이 15가 아니면 그 값으로 교체
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% --- 고DOF freeform .ent 생성(RBF + 절단면) ---
rng('shuffle');  charSet=['a':'z' 'A':'Z' '0':'9'];
tag = charSet(randi(numel(charSet),1,10));
entPath = [FF_BASE 'asym_' tag '.1.ent'];
[ok, maxDraft] = generate_freeform_ent(pShape, NRBF, NCUT, FF_TEMPLATE, entPath);
if ~ok || ~exist(entPath,'file'), output=fail_output(); output.loadFailed=true; return; end
List=ltml.LTDbList(lt,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
Key=ltml.LTListByName(lt,List,'LibraryElement');
ltml.LTDbSet(lt,Key,'Filename', entPath);
List=ltml.LTDbList(lt,'LENS_MANAGER[1]','TEXTURE_PARAMETER');
Key=ltml.LTListByName(lt,List,'StretchZ');
ltml.LTDbSet(lt,Key,'Value', stretchZ);

% --- 나노 CPS (구조 변경: ITO -> Ag, 마지막 두께 = dAg) ---
load('nk_JH_total.mat');  load('Photopic_400_800.mat');  load('CIE_1931.mat');  load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';
wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98; horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);
no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_Ag_McPeak 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_Ag_McPeak 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);  cos089=cosd(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL dAg];      % <-- 마지막층 = dAg (변수)
EML_position=4; z0=12.5; u_data_num=499; max_u=3;
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

% --- 파장 루프: 시뮬 + 2D far-field(45x36) 누적 ---
nLat=45; nLong=36;  K=(wavelength_num-1)/n+1;
Power_output=zeros(1,wavelength_num);  Igrids=cell(1,K);
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
loadFailed = ~any(isfinite(Power_output)) || all(Power_output(1:n:end)<=0);

weight_factor_2=zeros(K,1); Power_output_2=zeros(K,1); EQE_sub_matrix_2=zeros(K,1);
for k=1:K
    idx=n*(k-1)+1;
    weight_factor_2(k)=weight_factor(idx);
    Power_output_2(k)=Power_output(idx);
    EQE_sub_matrix_2(k)=CPS_result.EQE_sub_matrix(idx);
end
EQE_wv_matrix=Power_output_2./weight_factor_2;
EQE_sub_matrix_2=EQE_sub_matrix_2/sum(EQE_sub_matrix_2)*EQE_sub_CPS;
contrib=EQE_wv_matrix.*EQE_sub_matrix_2;
EQE_total=sum(contrib);

thC=((1:nLat)-0.5)*(90/nLat);
phC=-180+((1:nLong)-0.5)*(360/nLong);
sint=sind(thC(:));
Wacc=zeros(nLat,nLong);
for k=1:K
    Wk=Igrids{k}.*sint;  sk=sum(Wk(:));
    if sk>0, Wacc=Wacc+contrib(k)*(Wk/sk); end
end
[PWin, phiC, contrast]=detect_phi_window(Wacc,thC,phC,TH_LO,TH_HI,PHI_W);

output=struct('EQE_region',PWin,'EQE_total',EQE_total,'phiC',phiC, ...
    'contrast',contrast,'thBand',[TH_LO TH_HI],'phiWidth',PHI_W, ...
    'maxDraft',maxDraft,'loadFailed',loadFailed);
LAST_METRICS = struct('phiC',phiC,'contrast',contrast,'draft',maxDraft);   % live 표시용
fprintf('[obj] EQEreg=%.5g | EQEtot=%.5g | φc=%+.0f° | ctr=%.2f | drft=%.0f° | dAg=%.1f | cell=%.1fx%.1f dec=(%.2f,%.2f)%s\n', ...
    PWin, EQE_total, phiC, contrast, maxDraft, dAg, rx_pat, ry_pat, DecenterX, DecenterY, ternary(loadFailed,' [LOAD FAIL]',''));

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');  Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');  Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d',count) ']" Delete= \Q']);
fclose('all');
end

function output = fail_output()
output=struct('EQE_region',0,'EQE_total',0,'phiC',NaN,'contrast',0, ...
    'thBand',[NaN NaN],'phiWidth',NaN,'maxDraft',NaN,'loadFailed',false);
end
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end


%% =====================================================================
%%  형상: raised-cosine 베이스 + RBF 요철 + 절단면 -> .ent (단위높이)
%% =====================================================================
function [ok, maxDraftDeg] = generate_freeform_ent(pShape, Nrbf, Ncut, templatePath, outPath)
ok=false; maxDraftDeg=NaN;
try
    Ra=1.2139; Rap=1.0; n=141; tbase=0.30;
    g=linspace(-Ra,Ra,n); [X,Y]=meshgrid(g,g); r=hypot(X,Y);
    H = 0.5*(1+cos(pi*min(max(r/Rap,0),1)));  H(r>Rap)=0;   % 베이스(단위높이)
    for i=1:Nrbf                                            % RBF 국소 요철
        b=4*(i-1);
        x0=pShape(b+1); y0=pShape(b+2); amp=pShape(b+3); sig=pShape(b+4);
        H = H + amp.*exp(-((X-x0).^2+(Y-y0).^2)./(2*sig^2));
    end
    W=ones(size(r)); rw=0.85*Rap;                          % rim window(테두리 0)
    z=(r-rw)/(Rap-rw); m=(r>=rw)&(r<Rap);
    W(m)=1-(3*z(m).^2-2*z(m).^3); W(r>=Rap)=0;
    H=max(H.*W,0);
    base=4*Nrbf;                                           % 절단면(날카로운 단면)
    for k=1:Ncut
        z0=pShape(base+3*(k-1)+1); mm=pShape(base+3*(k-1)+2); phi0=pShape(base+3*(k-1)+3);
        H = min(H, z0 + mm*(X*cos(phi0) + Y*sin(phi0)));
    end
    H=max(H,0);
    [gx,gy]=gradient(H,g,g);  slope=atand(hypot(gx,gy));
    maxDraftDeg=max(slope(r<=Rap));
    % .ent 쓰기
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
    ok=true;
catch me
    fprintf('[Geom] .ent 생성 실패: %s\n', me.message);
end
end


%% ===== phi 창 자동검출 =====
function [PWin, phiC, contrast] = detect_phi_window(Wacc, thC, phC, thLo, thHi, phiWidth)
tm=(thC>=thLo)&(thC<=thHi);
band=sum(Wacc(tm,:),1);
nL=numel(phC); half=phiWidth/2;
winP=zeros(1,nL);
for c=1:nL
    d=abs(mod(phC-phC(c)+180,360)-180);
    winP(c)=sum(band(d<=half));
end
[PWin,ic]=max(winP);  phiC=phC(ic);
dOpp=abs(mod(phC-(phiC+180)+180,360)-180);
POpp=sum(band(dOpp<=half));
contrast=PWin/max(POpp,eps);
end


%% ===== RenewLightTools_single =====
function RenewLightTools_single(modelPath)
global ID_LT ltml ltloc
lt_exe_path='C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
fprintf('--- Restarting LightTools (single) ---\n');
target_user='jhkim';
find_cmd=sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);
system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
t0=tic;
while toc(t0)<10
    [~,cmdout]=system(find_cmd);
    if ~contains(cmdout,'lt.exe'), break; end
    pause(0.3);
end
system(sprintf('start /min "" "%s" "%s"', lt_exe_path, modelPath));
try
    ltml =actxserver('ltcom64.LTAPI2');
    ltloc=actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패. 라이선스/설치 확인.');
end
t1=tic; found=false;
while toc(t1)<20
    [status,cmdout]=system(find_cmd);
    if status==0 && contains(cmdout,'lt.exe'), found=true; break; end
    pause(0.3);
end
if ~found, error('lt.exe 탐색 실패'); end
tokens=regexp(cmdout,'"(\d+)"','tokens');
ID_LT=str2double(tokens{1}{1});
fprintf('PID(LT)=%d\n', ID_LT);
tR=tic; ready=false;
while toc(tR)<20
    try
        lt=ltloc.GetLTAPI(ID_LT);
        ltml.LTCmd(lt,'Message "Check Connection"');
        ready=true; break;
    catch
        pause(0.5);
    end
end
if ~ready, fprintf('[경고] COM 준비 확인 실패.\n'); end
end
