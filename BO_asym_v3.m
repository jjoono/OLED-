% ============================================================
%  BO_asym_v3.m  -- 고자유도(30+ DOF) 비대칭 freeform + detuned 캐비티 공동최적화
%
%  [v2 대비 대변경]
%   (1) 옵티마이저: bayesopt(GP, ~20 DOF 한계) -> surrogateopt (내장, 비싼함수
%       전용 RBF surrogate, 연속 ~100변수). 30+ DOF 에서도 동작.
%   (2) 형상: 3 DOF 비대칭 돔 -> "raised-cosine 베이스 + Nrbf 국소 RBF 요철"
%       고자유도 3D freeform. 기본 Nrbf=8 -> 형상 32 DOF (+4 = 36 변수).
%       single-valued height field + 사각격자 + rim window -> 항상 유효 solid.
%   (3) SweptEntity(LT1/ID_swept) 완전 제거 -> 단일 인스턴스(배열모델만).
%       .ent 는 MATLAB 이 직접 생성. RepairEntities 불필요.
%   (4) 로드 실패 런타임 가드: 시뮬 후 TotalPower 가 0/비유한 이면 .ent 로드
%       실패로 판정 -> 목적함수 NaN 반환 -> surrogateopt 가 그 점을 버리고 회피.
%       (swept 안전망 없이도 하드 실패가 전체 run 을 죽이지 않음)
%
%  [변수] 형상 4*Nrbf + 4개(dETL,dHTL,stretchZ,Decenter)
%    - 형상: 각 RBF 마다 [x0,y0,amp,sigma]. 임의 위치 볼록/오목 요철 -> 비대칭 자유형
%    - dETL,dHTL: 마이크로캐비티 두께 / stretchZ: 렌즈 전체 높이(형상은 단위높이)
%    - Decenter: 편심 (PLANAR_REFERENCE_SURFACE.X = 15+Decenter, 기존 배선)
%
%  [목적함수] phi 40도 창(자동검출) 안, theta[40,60] 로 몰리는 EQE(파워) 절대값.
%
%  ================= @@CONFIRM =================
%   (C1) 형상 = raised-cosine + RBF freeform, Nrbf=8 (형상 32 DOF). Nrbf 로 조절.
%   (C2) 높이는 stretchZ 변수로만. 형상 .ent 는 단위높이(정점~1). h 는 변수 아님.
%   (C3) 목표 theta[40,60], phi 창폭 40도.
%   (C4) freeform_template_v2.ent 가 BASE 에 있어야 함(복사 완료).
%   (C5) r_pat = 15 고정.
%   (C6) far-field mesh = 위치3, 90x36, phi 중심 -175:10:175.
%   (C7) 병렬평가 불가(단일 COM 인스턴스). 64thread/512GB 는 레이트레이서용.
% ============================================================
clear;

%% ===== LightTools 연결 (단일 인스턴스) =====
global ID_LT ltml ltloc count ray_nums_current LT2 r_pat
global FF_TEMPLATE FF_BASE RESTART_INTERVAL eval_count TH_LO TH_HI PHI_W NRBF
tic;
FF_BASE     = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
LT2         = [FF_BASE 'assymetric_test.1.lts'];       % 배열 모델(유일 인스턴스)
FF_TEMPLATE = [FF_BASE 'freeform_template_v2.ent'];

RenewLightTools_single(LT2);
toc;
count = 1;
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 고정 설정 =====
RESTART_INTERVAL = 20;   % 시뮬 N회마다 재시작
eval_count       = 0;
r_pat            = 15;   % (C5) 렌즈 피치 고정
TH_LO = 40;  TH_HI = 60;  PHI_W = 40;    % (C3) 목표 theta 밴드 / phi 창폭
NRBF  = 8;               % (C1) RBF 개수 -> 형상 DOF = 4*NRBF

RAY_SEARCH  = 20000;     % 탐색용 ray
RAY_FINAL   = 100000;    % 최종 검증용 ray
N_FINAL_REP = 3;
MAXEVAL     = 400;       % surrogateopt 총 평가 예산 (36 DOF 고려; 시간 보고 조절)

%% ===== 변수 (형상 4*NRBF + 4) =====
shapeNames = cell(1, 4*NRBF);
lbS = zeros(1, 4*NRBF);  ubS = zeros(1, 4*NRBF);
for i = 1:NRBF
    b = 4*(i-1);
    shapeNames{b+1}=sprintf('x0_%d',i);  lbS(b+1)=-0.8; ubS(b+1)=0.8;
    shapeNames{b+2}=sprintf('y0_%d',i);  lbS(b+2)=-0.8; ubS(b+2)=0.8;
    shapeNames{b+3}=sprintf('amp_%d',i); lbS(b+3)=-0.45;ubS(b+3)=0.55;
    shapeNames{b+4}=sprintf('sig_%d',i); lbS(b+4)=0.15; ubS(b+4)=0.55;
end
varNames = [shapeNames, {'dETL','dHTL','stretchZ','Decenter'}];
lb = [lbS,  10,  10, 0.1, 0.0];
ub = [ubS, 150, 150, 3.0, 7.5];
NV = numel(lb);
fprintf('총 변수 = %d (형상 %d = 4*%d RBF + dETL,dHTL,stretchZ,Decenter)\n', NV, 4*NRBF, NRBF);

%% ===== surrogateopt 실행 =====
ray_nums_current = RAY_SEARCH;
optsSurr = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', MAXEVAL, ...
    'UseParallel', false, ...            % (C7) 단일 COM 인스턴스 -> 직렬
    'MinSurrogatePoints', max(2*NV+1, 40), ...
    'PlotFcn', [], ...
    'CheckpointFile', 'BO_asym_v3_ckpt.mat', ...   % 크래시 시 재개용
    'Display', 'iter');

fprintf('\n######## surrogateopt 시작: %d DOF, target θ∈[%d,%d], φ창=%d° ########\n', ...
    NV, TH_LO, TH_HI, PHI_W);
objfun = @(x) surr_objective(x);
[xBest, fBest, exitflag, outS] = surrogateopt(objfun, lb, ub, optsSurr); %#ok<ASGLU>

%% ===== 최종 고정밀 검증 =====
ray_nums_current = RAY_FINAL;
e = nan(1, N_FINAL_REP);
for rrep = 1:N_FINAL_REP, e(rrep) = simulate_metric(xBest); end
bestEQE = mean(e, 'omitnan');  bestStd = std(e, 'omitnan');

bd = objFcn_regionPower(xBest);
save('BO_asym_v3_result.mat', 'xBest', 'bestEQE', 'bestStd', 'bd', ...
    'varNames', 'lb', 'ub', 'r_pat', 'NRBF', 'TH_LO', 'TH_HI', 'PHI_W', 'outS');
fprintf('\n######## Done ########\n');
fprintf('  EQE_region(절대) = %.5g ± %.2g (surrogate fBest=%.5g)\n', bestEQE, bestStd, -fBest);
fprintf('  EQE_total=%.5g | φ중심=%+.0f° | φ대비비=%.2f | maxDraft=%.1f°\n', ...
    bd.EQE_total, bd.phiC, bd.contrast, bd.maxDraft);


%% =====================================================================
%%  surrogateopt 목적함수 (최소화) : -EQE_region
%% =====================================================================
function f = surr_objective(x)
m = simulate_metric(x);
if isnan(m)
    f = NaN;          % 실패 -> surrogateopt 가 버림
else
    f = -m;           % 최대화 -> 최소화
end
end

%% ===== 1회 평가 래퍼 (주기 재시작 + 크래시 -> NaN) =====
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
    if ~isfinite(m) || out.loadFailed, m = NaN; end   % 로드 실패 가드
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    m = NaN;
    RenewLightTools_single(LT2);
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end


%% =====================================================================
%%  Objective: RBF freeform .ent 직접생성 + 나노 CPS + phi 자동검출 절대 EQE
%% =====================================================================
function output = objFcn_regionPower(point)
global ID_LT ltml ltloc count ray_nums_current r_pat FF_TEMPLATE FF_BASE
global TH_LO TH_HI PHI_W NRBF

lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

% --- 변수 언팩 ---
pShape   = point(1:4*NRBF);
dETL     = point(4*NRBF+1);  dHTL = point(4*NRBF+2);
stretchZ = point(4*NRBF+3);  Decenter = point(4*NRBF+4);

% --- 고정 설정 ---
d_sub=1.3; r_OLED=1; Lensheight=0.01;
wavelength_start=580; wavelength_end=590; n=10;
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
ltml.LTDbSet(lt,Key,'Geometry_1',r_pat);  ltml.LTDbSet(lt,Key,'Geometry_2',r_pat);
List=ltml.LTDbList(lt,'lens_manager[1]','PLANAR_REFERENCE_SURFACE');
Key=ltml.LTListByName(lt,List,'ReferenceSurface');
ltml.LTDbSet(lt,Key,'X',15+Decenter);
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% --- RBF freeform .ent 직접 생성 + 텍스처에 물리기 ---
rng('shuffle');  charSet=['a':'z' 'A':'Z' '0':'9'];
tag = charSet(randi(numel(charSet),1,10));
entPath = [FF_BASE 'asym_' tag '.1.ent'];
[ok, maxDraft] = generate_rbf_ent(pShape, NRBF, FF_TEMPLATE, entPath);
if ~ok || ~exist(entPath,'file')
    output = fail_output();  output.loadFailed = true;  return;
end
List=ltml.LTDbList(lt,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
Key=ltml.LTListByName(lt,List,'LibraryElement');
ltml.LTDbSet(lt,Key,'Filename', entPath);
List=ltml.LTDbList(lt,'LENS_MANAGER[1]','TEXTURE_PARAMETER');
Key=ltml.LTListByName(lt,List,'StretchZ');
ltml.LTDbSet(lt,Key,'Value', stretchZ);

% --- 나노 CPS + 하단 반사 코팅 ---
load('nk_JH33.mat');  load('Photopic_400_800.mat');  load('CIE_1931.mat');  load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';
wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98; horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);
no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);  cos089=cosd(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL 150];
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

% --- 파장 루프: 시뮬 + 2D far-field 누적 ---
nLat=90; nLong=36;  K=(wavelength_num-1)/n+1;
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

% --- 로드 실패 가드: 파워가 전혀 안 나오면 .ent 로드 실패로 판정 ---
loadFailed = ~any(isfinite(Power_output)) || all(Power_output(1:n:end)<=0);

% --- 파장 가중 -> EQE_total + 2D 누적 Wacc ---
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
fprintf('[obj] EQE_region=%.5g | EQE_total=%.5g | φc=%+.0f° | contrast=%.2f | draft=%.0f° | dec=%.2f%s\n', ...
    PWin, EQE_total, phiC, contrast, maxDraft, Decenter, ternary(loadFailed,' [LOAD FAIL]',''));

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
%%  형상: raised-cosine 베이스 + Nrbf 국소 RBF 요철 -> .ent (단위높이)
%% =====================================================================
function [ok, maxDraftDeg] = generate_rbf_ent(pShape, Nrbf, templatePath, outPath)
ok=false; maxDraftDeg=NaN;
try
    Ra=1.2139; Rap=1.0; n=141; tbase=0.30;
    g=linspace(-Ra,Ra,n); [X,Y]=meshgrid(g,g); r=hypot(X,Y);
    % 베이스: raised-cosine 돔(테두리 높이·기울기 0 -> 최상 몰더빌리티), 단위높이
    H = 0.5*(1+cos(pi*min(max(r/Rap,0),1)));  H(r>Rap)=0;
    % RBF 국소 요철 (임의 위치 볼록/오목) -> 고자유도 비대칭
    for i=1:Nrbf
        b=4*(i-1);
        x0=pShape(b+1); y0=pShape(b+2); amp=pShape(b+3); sig=pShape(b+4);
        H = H + amp.*exp(-((X-x0).^2+(Y-y0).^2)./(2*sig^2));
    end
    % rim window: 테두리에서 0 강제(타일링) + RBF 누설 차단
    W=ones(size(r)); rw=0.85*Rap;
    z=(r-rw)/(Rap-rw); m=(r>=rw)&(r<Rap);
    W(m)=1-(3*z(m).^2-2*z(m).^3); W(r>=Rap)=0;
    H=max(H.*W,0);
    % 제작성 모니터: 조리개 내 최대 draft(경사) [deg]  (제약 아님, 로깅용)
    [gx,gy]=gradient(H,g,g);  slope=atand(hypot(gx,gy));
    maxDraftDeg = max(slope(r<=Rap));
    % .ent 쓰기 (Z=+H, rear=-tbase, SmoothResample No)
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


%% =====================================================================
%%  phi 창 자동검출
%% =====================================================================
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


%% =====================================================================
%%  RenewLightTools_single : 배열모델 1개만 재시작 (start /min + 폴링)
%% =====================================================================
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
