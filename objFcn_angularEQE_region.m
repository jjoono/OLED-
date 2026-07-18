function output = objFcn_angularEQE_region(point)
% OBJFCN_ANGULAREQE_REGION  v4 objFcn 의 (theta 고정밴드 x phi 자동검출) 확장판.
%
% [v4 대비 변경점 - 딱 3곳]
%   (1) far-field mesh(위치 3)를 1D(theta) 대신 2D(theta x phi)로 읽는다
%       -> read_ff_mesh2d(). mesh 는 Symmetry="No Symmetry" 필수. @@VERIFY
%   (2) 파장가중 누적그리드 Wacc 에서 phi 창을 자동검출 -> detect_phi_window().
%       theta 밴드는 응용에 맞게 아래 REGION_* 로 고정 (기본 40~60).
%   (3) bare 기준(같은 캐비티, 렌즈 없음, capture_bare_reference.m 로 1회 저장)과
%       "같은 창"에서 비교해 이득 분해를 출력:
%         G_abs = 창내 절대파워비,  G_extract = 총추출비,  G_redist = 분율비
%         (G_abs ≈ G_extract x G_redist ; 재분배 기여를 분리해 추적)
%
% 반환 output 필드:
%   .EQE_region  : 자동검출 창의 region EQE  <- BO 목적함수로 사용 권장
%   .phiC        : 검출된 phi 중심 [deg]
%   .contrast    : phi 대비비 / .EQE_total / .G_abs .G_extract .G_redist (bare 있을때)
%
% 형상 주입부(swept entity)는 v4 그대로 두었다 - 이동/편심 수정은 사용자가 진행.

global ID_LT ID_swept ltml ltloc count r_pat ray_nums_current

% ===== region 설정 (응용에 따라 조정) =====
REGION_TH_LO   = 40;    % [deg] theta 밴드 하한 (고각 응용 예시)
REGION_TH_HI   = 60;    % [deg] theta 밴드 상한
REGION_PHIW    = 80;    % [deg] phi 창 전체폭 (중심은 자동검출)
MESH_POS       = 3;     % far-field 2D mesh 위치
MESHCFG = struct('nLong',36,'nLat',45,'longMin',-180,'longMax',180, ...
                 'latMin',0,'latMax',90);   % @@VERIFY 모델 mesh bin 수와 일치
BARE_REF_MAT   = 'bare_ref.mat';            % capture_bare_reference.m 산출물

lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

d_sub = 1.3;  r_OLED = 1;  x_pattern = r_pat;  y_pattern = r_pat;
Lensheight = 0.01;
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
ltml.LTDbSet(lt,Key,'Geometry_1',x_pattern);  ltml.LTDbSet(lt,Key,'Geometry_2',y_pattern);
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% ---- 형상 파라미터 언팩 + swept entity 주입 (v4 그대로; 사용자 수정 예정) ----
x2=point(1); x3=point(2); x4=point(3); x5=point(4); x6=point(5);
y2=point(6); y3=point(7); y4=point(8); y5=point(9); y6=point(10);
dETL=point(11); dHTL=point(12); stretchZ=point(13);

xy = zeros(7,2);
xy(1,:)=[0,1]; xy(7,:)=[1,0];
xy(2,:)=[x2,y2]; xy(3,:)=[x3,y3]; xy(4,:)=[x4,y4]; xy(5,:)=[x5,y5]; xy(6,:)=[x6,y6];

lt  = ltloc.GetLTAPI(ID_swept);
ltx = getltpointer(ID_swept);
lt2 = ltloc.GetLTAPI(ID_LT);

Curve="LENS_MANAGER[1].COMPONENTS[Components].SWEPT_SOLID[SweptEntity].SWEPT_PRIMITIVE[SweptPrimitive].SWEPT_PROFILE[SweptProfile].FITTED_CURVE[SweptSurface_1]";
ltx.SetSweptProfilePoints(Curve,xy,7);
ltx.DbSet(Curve,'StartSlopeMode',"Auto");  ltx.DbSet(Curve,'EndSlopeMode',"Auto");
List=ltml.LTDbList(lt,'LENS_MANAGER[1]','FITTED_CURVE');
Key=ltml.LTListByName(lt,List,'SweptSurface_1');
ltml.LTDbSet(lt,Key,'NumFacets',100);
x_values=zeros(101,1);
for a=1:101, x_values(a)=ltml.LTDbGet(lt,Key,'YFacetsAt',a); end
if max(x_values) > 1, xy = xy / max(x_values); end
ltx.SetSweptProfilePoints(Curve,xy,7);
ltx.DbSet(Curve,'StartSlopeMode',"Auto");  ltx.DbSet(Curve,'EndSlopeMode',"Auto");
xy_l=zeros(7,2);
for j=1:7
    xy_l(j,1)=ltml.LTDbGet(lt,Key,'YAt',j);
    xy_l(j,2)=ltml.LTDbGet(lt,Key,'ZAt',j);
end
if max(abs(xy(:)-xy_l(:))) > 1e-4
    output = fail_output();  return;
end

rng('shuffle');
charSet=['a':'z' 'A':'Z' '0':'9'];
index=charSet(randi(numel(charSet),1,10));
pathname='"C:\Users\jhkim\Desktop\Green_CE_Calculation\swept_';
pathname_unrepaired='"C:\Users\jhkim\Desktop\Green_CE_Calculation\unrepaired\swept_unrepaired_';
totalpath=[pathname index '.ent"'];
totalpath_unrepaired=[pathname_unrepaired index '.ent"'];
ltml.LTCmd(lt,'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt,sprintf('SaveLibrary XYZ 0,0,0 %s ',totalpath_unrepaired));
ltml.LTCmd(lt,'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt,'RepairEntities');
ltml.LTSetOption(lt,"ShowFileDialogBox",0);
ltml.LTSetOption(lt2,"ShowFileDialogBox",0);
ltml.LTCmd(lt,'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt,sprintf('SaveLibrary XYZ 0,0,0 %s ',totalpath));
ltml.LTCmd(lt,'Undo');  ltml.LTCmd(lt,'Undo');
totalpathmod=[pathname index '.1.ent"'];
List=ltml.LTDbList(lt2,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
Key=ltml.LTListByName(lt2,List,'LibraryElement');
ltml.LTDbSet(lt2,Key,'Filename',totalpathmod);
List=ltml.LTDbList(lt2,'LENS_MANAGER[1]','TEXTURE_PARAMETER');
Key=ltml.LTListByName(lt2,List,'StretchZ');
ltml.LTDbSet(lt2,Key,'Value',stretchZ);

% ---- 나노 CPS + 하단반사 코팅 (v4 그대로) ----
load('nk_JH33.mat');  load('Photopic_400_800.mat');  load('CIE_1931.mat');  load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';
wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98;  horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);
no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL 150];
EML_position=4;  z0=12.5;  u_data_num=499;  max_u=3;
CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_sub_CPS=CPS_result.EQE_sub;
TMF_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
Reflectance=(abs(TMF_p.r_p).^2+abs(TMF_s.r_s).^2)/2;
lt = ltloc.GetLTAPI(ID_LT);
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

% ===== 파장 루프: (1) 2D mesh 읽기로 교체 =====
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
    [Igrids{(wv+n-1)/n}, thC, phC] = read_ff_mesh2d(ltml, lt, MESH_POS, MESHCFG);
end

% ---- 파장 가중 (v4 그대로) ----
weight_factor_2=zeros(K,1);  Power_output_2=zeros(K,1);  EQE_sub_matrix_2=zeros(K,1);
for k=1:K
    idx=n*(k-1)+1;
    weight_factor_2(k)=weight_factor(idx);
    Power_output_2(k)=Power_output(idx);
    EQE_sub_matrix_2(k)=CPS_result.EQE_sub_matrix(idx);
end
EQE_wv_matrix=Power_output_2./weight_factor_2;
EQE_sub_matrix_2=EQE_sub_matrix_2/sum(EQE_sub_matrix_2)*EQE_sub_CPS;
contrib=EQE_wv_matrix.*EQE_sub_matrix_2;   % 파장별 EQE 기여
EQE_total=sum(contrib);

% ===== (2) 파장가중 누적 Wacc + phi 자동검출 =====
% 각 파장 그리드를 "그 파장의 EQE 기여" 로 정규화-가중해 합산:
%   Wacc(theta,phi) 총합 = EQE_total 이 되도록 -> 창내 합 = EQE_region 그 자체.
Wacc=zeros(MESHCFG.nLat,MESHCFG.nLong);
sint=sind(thC(:));
for k=1:K
    Wk=Igrids{k}.*sint;                    % 파워 가중
    sk=sum(Wk(:));
    if sk>0, Wacc=Wacc+contrib(k)*(Wk/sk); end
end
R=detect_phi_window(Wacc,thC,phC,REGION_TH_LO,REGION_TH_HI,REGION_PHIW);

output=struct();
output.EQE_total  = EQE_total;
output.EQE_region = R.PWin;        % 자동검출 창의 EQE  <- BO 목적함수
output.phiC       = R.phiC;
output.fracWin    = R.fracWin;
output.contrast   = R.contrast;
output.thBand     = [REGION_TH_LO REGION_TH_HI];
output.phiWidth   = REGION_PHIW;

% ===== (3) bare 기준과 "같은 창" 비교 (있을 때만) =====
% bare_ref.mat: capture_bare_reference.m 가 저장 (Wacc_bare, EQE_total_bare 등).
% 같은 (dETL,dHTL) 의 bare 와 비교해야 캐비티 효과가 아니라 렌즈 효과가 분리된다.
persistent BARE
if isempty(BARE) && exist(BARE_REF_MAT,'file'), BARE=load(BARE_REF_MAT); end
if ~isempty(BARE)
    tm=(thC>=REGION_TH_LO)&(thC<=REGION_TH_HI);
    d=abs(mod(phC-R.phiC+180,360)-180);  pm=d<=REGION_PHIW/2;
    PWin_bare=sum(sum(BARE.Wacc(tm,pm)));
    output.G_abs     = R.PWin / max(PWin_bare,eps);                 % 창내 절대이득
    output.G_extract = EQE_total / max(BARE.EQE_total,eps);         % 총추출 이득
    output.G_redist  = output.G_abs / max(output.G_extract,eps);    % phi 재분배 이득
    fprintf('[Region] phiC=%+.0f | EQE_reg=%.4f | G_abs=%.2fx = G_ext %.2fx * G_redist %.2fx | contrast=%.2f\n', ...
        R.phiC, output.EQE_region, output.G_abs, output.G_extract, output.G_redist, output.contrast);
else
    fprintf('[Region] phiC=%+.0f | EQE_reg=%.4f | contrast=%.2f (bare_ref.mat 없음: 이득분해 생략)\n', ...
        R.phiC, output.EQE_region, output.contrast);
end

% 코팅 정리 (v4 그대로)
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');  Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');  Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d',count) ']" Delete= \Q']);
fclose('all');
end

function output = fail_output()
output=struct('EQE_total',0,'EQE_region',0,'phiC',NaN,'fracWin',0, ...
    'contrast',0,'thBand',[NaN NaN],'phiWidth',NaN);
end
