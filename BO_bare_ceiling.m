% ============================================================
%  BO_bare_ceiling.m  -- 렌즈 없이 "소자구조(캐비티)만"으로 EQE_region 상한 계산
%
%  [목적] v4(렌즈+캐비티 공동최적화)의 비교 기준선.
%    렌즈 자리에 평면(flat) .ent 를 물려 "렌즈 없음"으로 만들고, 캐비티 3변수
%    (dETL,dHTL,dAg)만 최적화해 target 영역(θ[40,60] x φ 40°창)의 EQE_region
%    "절대 최대값"을 구한다. 렌즈가 이 값을 넘어야 렌즈의 기여가 입증됨.
%
%  [주의 - 물리] bare 소자는 방위각(φ) 대칭이라 특정 φ 로 못 몬다. 따라서
%    EQE_region(bare) ≈ (θ밴드 EQE) x (φ창폭/360°) 수준이 상한이다. 이 값을
%    v4 렌즈 결과와 비교하면 "φ 접힘 이득(fold gain)"이 순수하게 드러난다.
%
%  [구조/변수] Ag 스택(= v4 와 동일), thickness=[100 dETL 25 10 dHTL dAg].
%    변수 = dETL[10,150], dHTL[10,150], dAg[0,50].  렌즈 관련 변수 없음.
%  [옵티마이저] surrogateopt (3 DOF 라 가벼움). 실시간 진행표시 포함.
% ============================================================
clear;
global ID_LT ltml ltloc count ray_nums_current LT2 r_pat FF_TEMPLATE FF_BASE
global TH_LO TH_HI PHI_W RESTART_INTERVAL eval_count LAST_METRICS SUMMARY_EVERY MAXEVAL
tic;
FF_BASE     = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
LT2         = [FF_BASE 'assymetric_test.1.lts'];
FF_TEMPLATE = [FF_BASE 'freeform_template_v2.ent'];

RenewLightTools_single(LT2);
toc;
count = 1;  eval_count = 0;  RESTART_INTERVAL = 20;
r_pat = 15;  TH_LO = 40;  TH_HI = 60;  PHI_W = 40;
LAST_METRICS = struct('phiC',NaN,'contrast',NaN,'draft',0);  SUMMARY_EVERY = 10;
RAY_SEARCH = 20000;  RAY_FINAL = 100000;  N_FINAL_REP = 3;  MAXEVAL = 120;

lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 렌즈 없음: 평면 .ent 1회 생성 + 물리기 =====
flatEnt = [FF_BASE 'flat_nolens.1.ent'];
if ~make_flat_ent(FF_TEMPLATE, flatEnt)
    error('평면 .ent 생성 실패 - 템플릿 경로 확인');
end
List=ltml.LTDbList(lt,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
Key=ltml.LTListByName(lt,List,'LibraryElement');  ltml.LTDbSet(lt,Key,'Filename', flatEnt);
List=ltml.LTDbList(lt,'LENS_MANAGER[1]','TEXTURE_PARAMETER');
Key=ltml.LTListByName(lt,List,'StretchZ');  ltml.LTDbSet(lt,Key,'Value', 1);
% 편심 없음: 레퍼런스면 기본 위치(X=15, Decenter=0)
List=ltml.LTDbList(lt,'lens_manager[1]','PLANAR_REFERENCE_SURFACE');
Key=ltml.LTListByName(lt,List,'ReferenceSurface');  ltml.LTDbSet(lt,Key,'X',15);

%% ===== surrogateopt: 캐비티 3변수만 =====
varNames = {'dETL','dHTL','dAg'};
lb = [ 10,  10,  0];
ub = [150, 150, 50];
ray_nums_current = RAY_SEARCH;
opts = optimoptions('surrogateopt','MaxFunctionEvaluations',MAXEVAL, ...
    'UseParallel',false,'MinSurrogatePoints',20,'OutputFcn',@live_progress, ...
    'PlotFcn',[],'CheckpointFile','BO_bare_ckpt.mat','Display','iter');
fprintf('\n######## bare ceiling: 캐비티만(dETL,dHTL,dAg), θ∈[%d,%d] φ창=%d° ########\n', ...
    TH_LO, TH_HI, PHI_W);
[xBest,fBest,~,outS] = surrogateopt(@bare_obj, lb, ub, opts); %#ok<ASGLU>

%% ===== 최종 고정밀 검증 =====
ray_nums_current = RAY_FINAL;
e = nan(1,N_FINAL_REP);
for rrep=1:N_FINAL_REP, e(rrep)=bare_metric(xBest); end
bareEQE = mean(e,'omitnan');  bareStd = std(e,'omitnan');
bd = objFcn_bare(xBest);
save('BO_bare_ceiling_result.mat','xBest','bareEQE','bareStd','bd','varNames','lb','ub', ...
    'TH_LO','TH_HI','PHI_W','outS');
fprintf('\n######## Bare ceiling ########\n');
fprintf('  EQE_region(bare, 상한) = %.5g ± %.2g\n', bareEQE, bareStd);
fprintf('  최적 캐비티: dETL=%.1f dHTL=%.1f dAg=%.1f\n', xBest(1),xBest(2),xBest(3));
fprintf('  EQE_total=%.5g | θ밴드 EQE=%.5g | EQE_region/θ밴드 = %.3f (φ균일이면 ~%.3f)\n', ...
    bd.EQE_total, bd.EQE_theta, bd.EQE_region/max(bd.EQE_theta,eps), PHI_W/360);
fprintf('  --> v4 렌즈 EQE_region 이 %.5g 를 넘으면 그 초과분이 순수 렌즈 이득(φ fold).\n', bareEQE);


%% ===== 실시간 진행 =====
function stop = live_progress(~, optimValues, state)
stop=false; persistent t0 fc eqe best bs figH ax
global SUMMARY_EVERY MAXEVAL
switch state
    case 'init'
        t0=tic; fc=[]; eqe=[]; best=[]; bs=-inf;
        figH=figure('Name','bare ceiling 진행','Color','w'); ax=axes(figH); hold(ax,'on'); grid(ax,'on');
        xlabel(ax,'평가'); ylabel(ax,'EQE\_region(bare)'); title(ax,'파랑=매평가, 빨강=best');
    case 'done'
        fprintf('\n===== [완료] %d 평가 | bare 상한 EQE_region=%.4g | %.1f분 =====\n', ...
            optimValues.funccount, bs, toc(t0)/60);
    otherwise
        f=optimValues.fval; v=-f; if ~isfinite(f), v=NaN; end
        fc(end+1)=optimValues.funccount; eqe(end+1)=v;
        if isfinite(v)&&v>bs, bs=v; end
        best(end+1)=bs;
        if ~isempty(figH)&&isvalid(figH)
            cla(ax); plot(ax,fc,eqe,'b.','MarkerSize',8); plot(ax,fc,best,'r-','LineWidth',1.5); drawnow limitrate;
        end
        if mod(numel(fc),SUMMARY_EVERY)==0
            el=toc(t0); rate=numel(fc)/max(el,eps); eta=(MAXEVAL-optimValues.funccount)/max(rate,eps);
            fprintf('\n===== [진행] %d/%d | bare best EQE_region=%.4g | %.1f분, ETA~%.1f분 =====\n\n', ...
                optimValues.funccount, MAXEVAL, bs, el/60, eta/60);
        end
end
end

function f = bare_obj(x)
m=bare_metric(x); if isnan(m), f=NaN; else, f=-m; end
end
function m = bare_metric(pt)
global ID_LT ltml ltloc eval_count RESTART_INTERVAL LT2
eval_count=eval_count+1;
if mod(eval_count,RESTART_INTERVAL)==0
    RenewLightTools_single(LT2);
    lt=ltloc.GetLTAPI(ID_LT); ltml.LTSetOption(lt,"ShowFileDialogBox",0);
    % 재시작 후 평면 .ent 재설정
    global FF_BASE
    flatEnt=[FF_BASE 'flat_nolens.1.ent'];
    List=ltml.LTDbList(lt,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
    Key=ltml.LTListByName(lt,List,'LibraryElement'); ltml.LTDbSet(lt,Key,'Filename',flatEnt);
    List=ltml.LTDbList(lt,'lens_manager[1]','PLANAR_REFERENCE_SURFACE');
    Key=ltml.LTListByName(lt,List,'ReferenceSurface'); ltml.LTDbSet(lt,Key,'X',15);
    pause(1);
end
try
    out=objFcn_bare(pt); m=out.EQE_region;
    if ~isfinite(m)||out.loadFailed, m=NaN; end
catch err
    fprintf('\n[Error] bare eval %d 충돌: %s\n', eval_count, err.message);
    m=NaN; RenewLightTools_single(LT2);
    lt=ltloc.GetLTAPI(ID_LT); ltml.LTSetOption(lt,"ShowFileDialogBox",0);
end
end


%% ===== bare objective: 캐비티만 바꿔 시뮬 (평면 .ent 고정) =====
function output = objFcn_bare(point)
global ID_LT ltml ltloc count ray_nums_current r_pat TH_LO TH_HI PHI_W
lt=ltloc.GetLTAPI(ID_LT); ltml.LTSetOption(lt,"ShowFileDialogBox",0);
dETL=point(1); dHTL=point(2); dAg=point(3);

d_sub=1.3; r_OLED=1; Lensheight=0.01;
wavelength_start=580; wavelength_end=590; n=10;
if isempty(ray_nums_current), ray_nums=50000; else, ray_nums=ray_nums_current; end
List=ltml.LTDbList(lt,'lens_manager[1]','SIMULATIONS'); Key=ltml.LTListByName(lt,List,'ForwardAll');
ltml.LTDbSet(lt,Key,'MaxProgress',ray_nums);
List=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE'); Key=ltml.LTListByName(lt,List,'Substrate');
ltml.LTDbSet(lt,Key,'Height',d_sub); ltml.LTDbSet(lt,Key,'Y',d_sub/2);
SRList=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE'); SRKey=ltml.LTListAtPos(lt,SRList,2);
ltml.LTDbSet(lt,SRKey,'Y',d_sub+Lensheight/2);
List=ltml.LTDbList(lt,'lens_manager[1]','TEXTURE_ZONE_EXTENT'); Key=ltml.LTListByName(lt,List,'zone');
ltml.LTDbSet(lt,Key,'Geometry_1',r_pat); ltml.LTDbSet(lt,Key,'Geometry_2',r_pat);
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE'); Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% --- CPS (Ag 스택, thickness 마지막=dAg) ---
load('nk_JH33.mat'); load('Photopic_400_800.mat'); load('CIE_1931.mat'); load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).'; wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98; horizontal_dipole_ratio=0.865; bottom_air_refractive_index=ones(wavelength_num,1);
no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_Ag_McPeak 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_Ag_McPeak 1.51*ones(401,1)];
layer_num=size(no_bar,2);  sin089=sind(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:); ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL dAg];
CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,4,12.5,499,3,wavelength);
EQE_sub_CPS=CPS_result.EQE_sub;
TMF_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
Reflectance=(abs(TMF_p.r_p).^2+abs(TMF_s.r_s).^2)/2;
fileID=fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa',count),'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0','DATANAME: R_Bottom_',count,'ABSORBING: YES','INDEX: 1.51','DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89, fprintf(fileID,'%s  %d  %d  %.3f\n','AOI',j,0,Reflectance(i-wavelength_start+1,j+1)); end
end
fclose(fileID);
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa',count) '"']);
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY'); Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE'); Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d',count));

I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
P_white=I_white.*repmat(sin089,wavelength_num,1); weight_factor=sum(P_white,2);

nLat=45; nLong=36; K=(wavelength_num-1)/n+1;
Power_output=zeros(1,wavelength_num); Igrids=cell(1,K);
for wv=1:n:wavelength_num
    fileID=fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1,90,0,0,360,90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','DISK_SOURCE'); SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power',weight_factor(wv));
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','Spectral_region'); SRKey=ltml.LTListAtPos(lt,SRList,2);
    ltml.LTDbSet(lt,SRKey,'Spectral_Definition','Monochromatic'); ltml.LTDbSet(lt,SRKey,'Single_Wavelength',wv+wavelength_start-1);
    List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER'); Key=ltml.LTListAtPos(lt,List,1);
    ltml.LTDbSet(lt,Key,'LoadFileName','C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt');
    ltml.LTBegin(lt); ltml.LTCmd(lt,'\V3D BeginAllSimulations'); ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH'); Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');
    Key=ltml.LTListAtPos(lt,List,3);
    Ig=zeros(nLat,nLong);
    for rr=1:nLat, for kk=1:nLong
        v=ltml.LTDbGet(lt,Key,'CellValue_UI',kk,rr); if isempty(v)||~isfinite(v), v=0; end
        Ig(rr,kk)=v;
    end, end
    Igrids{(wv+n-1)/n}=Ig;
end
loadFailed = ~any(isfinite(Power_output)) || all(Power_output(1:n:end)<=0);

weight_factor_2=zeros(K,1); Power_output_2=zeros(K,1); EQE_sub_matrix_2=zeros(K,1);
for k=1:K
    idx=n*(k-1)+1;
    weight_factor_2(k)=weight_factor(idx); Power_output_2(k)=Power_output(idx); EQE_sub_matrix_2(k)=CPS_result.EQE_sub_matrix(idx);
end
EQE_wv=Power_output_2./weight_factor_2;
EQE_sub_matrix_2=EQE_sub_matrix_2/sum(EQE_sub_matrix_2)*EQE_sub_CPS;
contrib=EQE_wv.*EQE_sub_matrix_2; EQE_total=sum(contrib);

thC=((1:nLat)-0.5)*(90/nLat); phC=-180+((1:nLong)-0.5)*(360/nLong); sint=sind(thC(:));
Wacc=zeros(nLat,nLong);
for k=1:K
    Wk=Igrids{k}.*sint; sk=sum(Wk(:)); if sk>0, Wacc=Wacc+contrib(k)*(Wk/sk); end
end
[PWin,phiC,contrast]=detect_phi_window(Wacc,thC,phC,TH_LO,TH_HI,PHI_W);
% θ밴드 전체(모든 φ) EQE = 참고용 (φ 균일이면 PWin ≈ EQE_theta * PHI_W/360)
tm=(thC>=TH_LO)&(thC<=TH_HI); EQE_theta=sum(sum(Wacc(tm,:)));

output=struct('EQE_region',PWin,'EQE_total',EQE_total,'EQE_theta',EQE_theta, ...
    'phiC',phiC,'contrast',contrast,'loadFailed',loadFailed);
fprintf('[bare] EQE_region=%.5g | EQE_theta=%.5g | EQE_total=%.5g | dETL=%.0f dHTL=%.0f dAg=%.1f%s\n', ...
    PWin, EQE_theta, EQE_total, dETL, dHTL, dAg, ternary(loadFailed,' [FAIL]',''));

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY'); Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE'); Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d',count) ']" Delete= \Q']);
fclose('all');
end

function s=ternary(c,a,b), if c, s=a; else, s=b; end, end

%% ===== 평면(flat) .ent 생성: Z=0 (렌즈 없음 = 평평한 계면) =====
function ok = make_flat_ent(templatePath, outPath)
ok=false;
try
    Ra=1.2139; n=141; tbase=0.30;
    g=linspace(-Ra,Ra,n); [X,Y]=meshgrid(g,g);
    Z=zeros(n,n);   % 완전 평면
    Xv=X(:); Yv=Y(:); Zv=Z(:); N=n*n;
    tpl=fileread(templatePath);
    tok=regexp(tpl,'ORAStartData;([\s\S]*?)ORAEndData;','tokenExtents');
    s0=tok{1}(1); e0=tok{1}(2);
    buf=sprintf('0 1 %d %d 0 0 %d 0 0 0',n,n,N);
    for i=1:N, buf=[buf sprintf(' %.17g %.17g %.17g',Xv(i),Yv(i),Zv(i))]; end %#ok<AGROW>
    buf=[buf ' 0 0 4 CartesianMapper 1 0 0 0 0'];
    newtxt=[tpl(1:s0-1) char(10) buf char(10) tpl(e0+1:end)];
    newtxt=regexprep(newtxt,'(CSGLensSurfacePrimitive_1[\s\S]*?setPosition:  \{ 0\. 0\. )[-0-9.eE]+(  \} ;)',['$1' num2str(-tbase,'%g') '$2'],'once');
    newtxt=regexprep(newtxt,'restoreSmoothResample: "Yes"','restoreSmoothResample: "No"','once');
    fid=fopen(outPath,'w'); fwrite(fid,newtxt); fclose(fid);
    ok=true;
catch me
    fprintf('[flat] 생성 실패: %s\n', me.message);
end
end

%% ===== phi 창 자동검출 =====
function [PWin, phiC, contrast] = detect_phi_window(Wacc, thC, phC, thLo, thHi, phiWidth)
tm=(thC>=thLo)&(thC<=thHi); band=sum(Wacc(tm,:),1);
nL=numel(phC); half=phiWidth/2; winP=zeros(1,nL);
for c=1:nL, d=abs(mod(phC-phC(c)+180,360)-180); winP(c)=sum(band(d<=half)); end
[PWin,ic]=max(winP); phiC=phC(ic);
dOpp=abs(mod(phC-(phiC+180)+180,360)-180); POpp=sum(band(dOpp<=half));
contrast=PWin/max(POpp,eps);
end

%% ===== RenewLightTools_single =====
function RenewLightTools_single(modelPath)
global ID_LT ltml ltloc
lt_exe_path='C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
fprintf('--- Restarting LightTools (single) ---\n'); target_user='jhkim';
find_cmd=sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);
system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
t0=tic; while toc(t0)<10, [~,c]=system(find_cmd); if ~contains(c,'lt.exe'), break; end, pause(0.3); end
system(sprintf('start /min "" "%s" "%s"', lt_exe_path, modelPath));
try, ltml=actxserver('ltcom64.LTAPI2'); ltloc=actxserver('ltlocator.Locator');
catch, error('LightTools 재시작 실패.'); end
t1=tic; found=false;
while toc(t1)<20, [st,c]=system(find_cmd); if st==0&&contains(c,'lt.exe'), found=true; break; end, pause(0.3); end
if ~found, error('lt.exe 탐색 실패'); end
tk=regexp(c,'"(\d+)"','tokens'); ID_LT=str2double(tk{1}{1}); fprintf('PID(LT)=%d\n', ID_LT);
tR=tic; ready=false;
while toc(tR)<20
    try, lt=ltloc.GetLTAPI(ID_LT); ltml.LTCmd(lt,'Message "Check Connection"'); ready=true; break;
    catch, pause(0.5); end
end
if ~ready, fprintf('[경고] COM 준비 확인 실패.\n'); end
end
