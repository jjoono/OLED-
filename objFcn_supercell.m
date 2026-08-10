function output = objFcn_supercell(hyp, seed)
% OBJFCN_SUPERCELL  무작위 조립(random supercell) MLA 1회 평가.
%   stress_random_mla.m 과 test_random_supercell.m 이 **같은 코드**를 쓰도록
%   별도 파일로 분리했다. 스크립트 안에 복사본을 두면 한쪽만 고쳐질 위험이 있다.
%
%   입력: hyp = [fill, rJitter, posJitter, aspect, aspectJitter, profileMix]
%         seed = 슈퍼셀 난수 시드 (결정론성 보장)
%   전역: ID_LT, ltml, ltloc, count, ray_nums_current, wave_n_current
global ID_LT ltml ltloc count ray_nums_current wave_n_current
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

% [주의] x_pattern/y_pattern 은 다른 캠페인 스크립트와 반드시 같아야 한다.
%   patch 크기가 다르면 EQE_total 이 통째로 달라져 계열 간 비교가 무의미해진다.
PATCH_XY = 25;
d_sub=1.295;  r_OLED=1;  x_pattern=PATCH_XY;  y_pattern=PATCH_XY;  Lensheight=0.01;

% [슈퍼셀 스케일] 슈퍼셀 하나에 nCols x nCols 개의 렌즈렛이 들어가므로, 텍스처
%   배치 간격을 같은 배수로 키우지 않으면 렌즈렛 하나의 물리 크기가 1/nCols 로
%   줄어든다. 그러면 "무작위성의 효과" 가 아니라 "렌즈렛이 작아진 효과" 를 재게 된다.
%   따라서 간격 = 기준 간격 x nCols 로 설정해 렌즈렛 크기를 보존한다.
SUPERCELL_NCOLS = 8;
SPACING_X0 = 0.0866;   % [mm] 단일 렌즈렛 기준 간격 (모델 .lts 의 초기값)
SPACING_Y0 = 0.1000;
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

%% !!! VERIFY IN LIGHTTOOLS: 텍스처 배치 간격 확대 (슈퍼셀 스케일 보존) !!!
%  모델 .lts 의 ZoneTextureHexagonalPlacement 는 setXSpacing/setYSpacing 을 갖는다.
%  DB 목록 이름 선례가 레포에 없어 후보를 순회하며 설정하고, 읽어서 검증한다.
%  설정에 실패하면 렌즈렛 크기가 1/nCols 로 줄어든 채 계산되므로 즉시 멈춘다.
tgtX = SPACING_X0 * SUPERCELL_NCOLS;
tgtY = SPACING_Y0 * SUPERCELL_NCOLS;
placedOK = false;
placeLists = {'ZONE_TEXTURE_HEXAGONAL_PLACEMENT','HEXAGONAL_PLACEMENT', ...
              'ZONE_TEXTURE_PLACEMENT','TEXTURE_PLACEMENT'};
for ip = 1:numel(placeLists)
    try
        PL = ltml.LTDbList(lt,'lens_manager[1]',placeLists{ip});
        PK = ltml.LTListAtPos(lt,PL,1);
        ltml.LTDbSet(lt,PK,'XSpacing',tgtX);
        ltml.LTDbSet(lt,PK,'YSpacing',tgtY);
        gx = ltml.LTDbGet(lt,PK,'XSpacing');
        gy = ltml.LTDbGet(lt,PK,'YSpacing');
        if abs(gx-tgtX) < 1e-9 && abs(gy-tgtY) < 1e-9
            placedOK = true;  break;
        end
    catch
        % 이 목록 이름은 이 모델에 없음 -> 다음 후보
    end
end
if ~placedOK
    error(['텍스처 배치 간격(XSpacing/YSpacing) 설정 실패. 슈퍼셀을 그대로 쓰면 ' ...
           '렌즈렛이 1/%d 크기가 되어 비교가 무의미해진다. LightTools Database ' ...
           'Browser 에서 ZoneTextureHexagonalPlacement 의 실제 DB 목록 이름을 ' ...
           '확인해 placeLists 에 추가할 것.'], SUPERCELL_NCOLS);
end

%% --- 슈퍼셀 .ent 생성 + unit-cell 교체 (objFcn_both 의 swept 블록 대체) ---
BASE = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
scDir = [BASE 'supercell_ents\'];
if ~exist(scDir,'dir'), mkdir(scDir); end

params = struct();
params.fill         = hyp(1);
params.rJitter      = hyp(2);
params.posJitter    = hyp(3);
params.aspect       = hyp(4);
params.aspectJitter = hyp(5);
params.profileMix   = hyp(6);
params.templatePath = [BASE 'freeform_template_v2.ent'];   % 원점정렬 검증 템플릿
params.nCols        = SUPERCELL_NCOLS;    % 위 배치 간격 확대와 반드시 같은 값
% 파일명: objFcn_both 의 swept_<tag> 관례를 따르되 결정론성은 (seed,hyp) 가 보장
rng('shuffle');
charSet = ['a':'z' 'A':'Z' '0':'9'];
index = charSet(randi(length(charSet), 1, 10));
params.outPath = [scDir 'supercell_' index '.1.ent'];

generate_random_supercell_ent(seed, params);   % 내부에서 rng(seed) -> 결정론적

lt2 = ltloc.GetLTAPI(ID_LT);
totalpathmod = ['"' params.outPath '"'];       % objFcn_both 와 동일하게 따옴표 포함
%% !!! VERIFY IN LIGHTTOOLS: supercell entity loading — check on first run !!!
%  아래 3 call 은 objFcn_both 이 swept_XXX.1.ent 를 물리던 검증된 패턴 그대로다.
%  다만 (i) SaveLibrary 산출물이 아닌 MATLAB 직접 작성 freeform .ent (SurfacePairLens/
%  FreeformEntity) 가 같은 unit-cell 에 물리는지, (ii) RepairEntities 없이 유효
%  solid 로 로드되는지, (iii) StretchZ=1 이 이 엔티티에 항등으로 작용하는지를
%  첫 실행에서 GUI 로 육안 확인할 것. (generate_freeform_ent.m 헤더가 같은 방식
%  연결을 전제하며, ent_tests/ 의 무작위 .ent 들이 수동 로드 검증된 전례.)
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'LIBRARY_ELEMENT_UNIT_CELL');
Key = ltml.LTListByName(lt2, List, 'LibraryElement');
ltml.LTDbSet(lt2, Key, 'Filename', totalpathmod);
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'TEXTURE_PARAMETER');
Key = ltml.LTListByName(lt2, List, 'StretchZ');
ltml.LTDbSet(lt2, Key, 'Value', 1);            % 높이는 .ent 에 이미 구워짐

%% CPS  (objFcn_both 와 동일; dETL/dHTL 은 이 가족의 변수가 아니므로 기준값 고정)
dETL = 60;  dHTL = 60;    % 스택 두께 기준값 (하이퍼벡터에 포함되지 않음)
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

%% bottom reflectance  (objFcn_both 와 동일)
TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
Reflectance=(R_p_bottom+R_s_bottom)/2;

%% Coating  (objFcn_both 와 동일)
lt = ltloc.GetLTAPI(ID_LT);
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0', 'DATANAME: R_Bottom_',count, 'ABSORBING: YES', 'INDEX: 1.51', 'DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n', 'AOI',j, 0, Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);   % LightTools 가 읽기 전에 플러시

ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

%% 파장 루프  (objFcn_both 와 동일)
I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
sin089=sind(0:89);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);
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

% 임시 슈퍼셀 .ent 정리 (디스크 누적 방지; 시뮬 종료 후이므로 안전)
try, delete(params.outPath); catch, end
end