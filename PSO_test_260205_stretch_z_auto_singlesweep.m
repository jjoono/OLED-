% MoE 관련 설정 및 GP, Acquisition Function 제거됨
% 순수 PSO + LightTools 연결 구조
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count
RenewLightTools();
% 기존에 연결된 세션이 있다면 재사용, 없으면 생성 (에러 방지용)
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
lt = ltloc.GetLTAPI(ID_swept); % swept entity
ltx= getltpointer(ID_swept);  % swept entity
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);


restart_interval=20;
eval_count=0;

% Define segment length and other necessary parameters
lt = ltloc.GetLTAPI(ID_LT);  % lenssizeeffect
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
d_sub=1.295;
r_OLED=1;
x_pattern=100;
y_pattern=100;
Lensheight=0.01;
wavelength_start=580;
wavelength_end=590;
n=10; % step size for wavelength
ray_nums=50000;

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

% passing input points
% point=[0.309194810202001	0.411173391675431	0.732395969768461	0.908340737923524	1	1.06709972561775	1.09409236301265	0.797932699506210	0.430532616637382	0.0854931028115802	44.5316455206464	49.7763138173905	2.27786560296760];
% point=[0.10559	0.26625	0.54841	0.26025	0.61016	1.1452	1.4330	1.2009	0.95352	0.74498	70	80	0.9235]; %JinHKim
% point=[0.369245997866821   0.538417655727910   0.551303271071089   0.468476063366523   0.788419473908105   0.949853677040117   0.930573493068931   0.773274091459503   0.343299563377685   0.336860702702061   71.7552036643924   85.3539450945621   2.12620442284361];
point=[0.258819045102521	0.500000000000000	0.707106781186548	0.866025403784439	0.965925826289068 0.965925826289068	0.866025403784439	0.707106781186548	0.500000000000000	0.258819045102521 65 70 1];
% point(6:10)=point(13)*point(6:10);

x2 = point(1);  x3 = point(2);  x4 = point(3);  x5 = point(4);  x6 = point(5);
y2 = point(6);  y3 = point(7);  y4 = point(8);  y5 = point(9);  y6 = point(10);
dETL = point(11); dHTL = point(12);
stretchZ = point(13);
% stretchZ=1;

% Create spline control points
xy = zeros(7,2);
xy(1,:) = [0, 1];
xy(7,:) = [1, 0];
xy(2,:) = [x2, y2];
xy(3,:) = [x3, y3];
xy(4,:) = [x4, y4];
xy(5,:) = [x5, y5];
xy(6,:) = [x6, y6];

lt = ltloc.GetLTAPI(ID_swept); % swept entity
ltx= getltpointer(ID_swept);  % swept entity
lt2 = ltloc.GetLTAPI(ID_LT); % LT simulation

Curve="LENS_MANAGER[1].COMPONENTS[Components].SWEPT_SOLID[SweptEntity].SWEPT_PRIMITIVE[SweptPrimitive].SWEPT_PROFILE[SweptProfile].FITTED_CURVE[SweptSurface_1]";
ltx.SetSweptProfilePoints(Curve,xy,7); % 7*2 double
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

if max_length > 1
    xy = xy / max_length;
end

ltx.SetSweptProfilePoints(Curve,xy,7); % 7*2 double
ltx.DbSet(Curve,'StartSlopeMode',"Auto");
ltx.DbSet(Curve,'EndSlopeMode',"Auto");

xy_l = zeros(7,2); % x,y coordinates in LightTools

for j=1:7
    xy_l(j,1) = ltml.LTDbGet(lt, Key, 'YAt', j);
    xy_l(j,2) = ltml.LTDbGet(lt, Key, 'ZAt', j);
end

if ~isequal(xy, xy_l)
    output = struct();
    output.EQE_0_20 = 0;
    output.EQE_20_40 = 0;
    output.EQE_40_60 = 0;
    output.EQE_60_80 = 0;
    output.EQE_total = 0;
    return;
end

% File name and path configuration
strLength = 10;
charSet = ['a':'z' 'A':'Z' '0':'9'];
numChars = length(charSet);
randIndices = randi(numChars, 1, strLength);
index = charSet(randIndices);

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

%% Define layer (CPS)
load('nk_JH33.mat');
load('Photopic_400_800.mat');
load('CIE_1931.mat');
load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';

wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98;
horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);

no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);
cos089=cosd(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL 150];

EML_position=4; % count from left side (+air)
z0=12.5;
u_data_num=997;
max_u=3;

CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_air_CPS=CPS_result.EQE_air;
EQE_sub_CPS=CPS_result.EQE_sub;

%% bottom reflectance
TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);

R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
T_p_bottom=no_bar(:,1)./no_bar(:,layer_num)*(1./cos089).*sqrt(1-(ne_bar(:,layer_num)./ne_bar(:,1)*sin089).^2).*abs(TMF_OLED_bottom_p.t_p).^2;

R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
T_s_bottom=no_bar(:,1)./no_bar(:,layer_num)*(1./cos089).*sqrt(1-(no_bar(:,layer_num)./no_bar(:,1)*sin089).^2).*abs(TMF_OLED_bottom_s.t_s).^2;

for i=1:wavelength_num
    T_p_bottom(i,ceil(asind(ne_bar(i,1)/ne_bar(i,layer_num)))+1:end)=0;
    T_s_bottom(i,ceil(asind(no_bar(i,1)/no_bar(i,layer_num)))+1:end)=0;
end

Transmittance=(T_p_bottom+T_s_bottom)/2;
Reflectance=(R_p_bottom+R_s_bottom)/2;

%% Coating (.mat to .coa)
lt = ltloc.GetLTAPI(ID_LT); % LT simulation
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0', 'DATANAME: R_Bottom_',count, 'ABSORBING: YES', 'INDEX: 1.51', 'DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n', 'AOI',j, 0, Reflectance(i-wavelength_start+1,j+1));
    end
end

ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

%%
I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p); % s랑 p 따로 구분하지 않음 일단
sin089=sind(0:89);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2); % I_white : I_sub의 파장별 intensity 301x90행렬
I_white_ang=sum(P_white);
%     weight_factor(1,1)=weight_factor(2,1);

wavelength_num=length(wavelength);

I_air_1_2=zeros(90,(wavelength_num+n-1)/n);
Luminance=cell((wavelength_num+n-1)/n,1);
Ray_wv=zeros(1,(wavelength_num+n-1)/n);
Cell_flux= zeros((wavelength_num+n-1)/n,9);
for wv=1:n:wavelength_num
    fileID = fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1, 90, 0, 0, 360, 90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt, 'Lens_manager[1]','DISK_SOURCE');
    SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power', weight_factor(wv)); % 파장에 따른 파워를 다르게 설정, 그 안에서 각도별 파워는 grid에서 조정
    for k=1:1  % 예전에 광원 많았을때는 k=1:광원수 였었음
        SRList=ltml.LTDbList(lt, 'Lens_manager[1]','Spectral_region');
        SRKey=ltml.LTListAtPos(lt,SRList,k+1);
        ltml.LTDbSet(lt,SRKey,'Spectral_Definition', 'Monochromatic');
        ltml.LTDbSet(lt,SRKey,'Single_Wavelength', wv+wavelength_start-1);
        List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');
        Key=ltml.LTListAtPos(lt,List,k);
        pathname='C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\'; % have to change pathname
        ltml.LTDbSet(lt,Key,'LoadFileName',[pathname sprintf('AI_temp.txt')]);
    end
    %% 시뮬레이션 및 후처리
    ltml.LTBegin(lt);
    ltml.LTCmd(lt,'\V3D BeginAllSimulations');
    ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');  % [W]
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,2);
    Power_output_30(wv)=ltml.LTDbGet(lt,Key,'TotalPower');  % [W]
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,3);
    for j=1:90
        I_air_1_JH(91-j,:)=ltml.LTDbGet(lt,Key,'CellValue_UI',1,91-j);
    end
    I_air_1_2(:,(wv+n-1)/n)=smooth(I_air_1_JH);
    %     I_air_1_2(:,(wv+n-1)/n)=I_air_1_JH;
end

K = (wavelength_num-1)/n + 1;

weight_factor_2  = zeros(K,1);
Power_output_2   = zeros(K,1);
EQE_sub_matrix_2 = zeros(K,1);

for k = 1:K
    idx = n*(k-1) + 1;

    weight_factor_2(k)  = weight_factor(idx);
    Power_output_2(k)   = Power_output(idx);
    EQE_sub_matrix_2(k) = CPS_result.EQE_sub_matrix(idx);
end

EQE_wv_matrix = Power_output_2 ./ weight_factor_2;  % (Kx1)

% 3) Normalize CPS spectral EQE_sub distribution to match EQE_sub_CPS
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;  % (Kx1)

% 4) Total EQE after optics
EQE_total = sum(EQE_wv_matrix .* EQE_sub_matrix_2);

% 5) Angular EQEs using LT angular intensity distribution per sampled wavelength
EQE_0_20   = 0;
EQE_20_40  = 0;
EQE_40_60  = 0;
EQE_60_80  = 0;

sin_col = sin089(:);  % 90x1 for elementwise multiply

for k = 1:K
    % Per-wavelength contribution to total EQE
    contrib_k = EQE_wv_matrix(k) * EQE_sub_matrix_2(k);

    % Angular radiant intensity vs theta for this wavelength sample
    I_theta = I_air_1_2(:,k);  % 90x1, theta = 0..89 deg

    % Convert to proportional angular power weights (constants cancel in fractions)
    W_theta = I_theta .* sin_col;  % 90x1, proportional to dP/dtheta integrated over azimuth
    W_tot   = sum(W_theta);

    % Fractions in bins (using [a,b) convention)
    f_0_20   = sum(W_theta(1:20))   / W_tot;  % 0..19 deg
    f_20_40  = sum(W_theta(21:40))  / W_tot;  % 20..39 deg
    f_40_60  = sum(W_theta(41:60))  / W_tot;  % 40..59 deg
    f_60_80  = sum(W_theta(61:80))  / W_tot;  % 60..79 deg

    % Accumulate angular EQEs
    EQE_0_20   = EQE_0_20   + contrib_k * f_0_20;
    EQE_20_40  = EQE_20_40  + contrib_k * f_20_40;
    EQE_40_60  = EQE_40_60  + contrib_k * f_40_60;
    EQE_60_80  = EQE_60_80  + contrib_k * f_60_80;
end

output = struct();
output.EQE_0_20 = EQE_0_20;
output.EQE_20_40 = EQE_20_40;
output.EQE_40_60 = EQE_40_60;
output.EQE_60_80 = EQE_60_80;
output.EQE_total = EQE_total;

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d', count) ']" Delete= \Q']);
fclose('all');


