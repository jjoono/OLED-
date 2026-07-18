%% capture_bare_reference.m
% bare 기준(같은 detuned 캐비티, 렌즈 없음)의 파장가중 far-field 그리드를 1회
% 계산해 bare_ref.mat 로 저장한다. objFcn_angularEQE_region 이 이 파일을 읽어
% "같은 창" 이득분해(G_abs = G_extract x G_redist)를 출력한다.
%
% [사용법]
%  1) LightTools 배열 모델에서 렌즈 텍스처를 비활성 (예: zone 크기 0 또는
%     LibraryElement 를 평판으로) - 이 부분은 GUI 에서 직접. @@VERIFY
%  2) 아래 dETL/dHTL 을 "비교할 캐비티와 동일하게" 설정 후 실행.
%     (bare 도 같은 캐비티여야 렌즈 효과만 분리됨)
%  3) bare_ref.mat 생성 확인 후 렌즈 텍스처 원복.
%
% [주의] BO 중 (dETL,dHTL) 이 변수로 움직이면 bare 도 함께 변해야 엄밀하지만,
% 평가마다 bare 재시뮬은 비용 2배라 비실용적. 실무 절차:
%   - BO 는 절대 EQE_region 최대화로 돌리고 (bare 불필요),
%   - 이득분해 리포트는 "고정 기준 캐비티" bare 대비로 해석하거나,
%   - 최종 best (dETL,dHTL) 로 이 스크립트를 한 번 더 돌려 최종 수치를 확정.

global ID_LT ltml ltloc count
if isempty(count), count = 999; end

% ===== 기준 캐비티 (비교 대상과 동일하게!) =====
dETL = 50;   % @@VERIFY
dHTL = 50;   % @@VERIFY
ray_nums = 100000;   % 기준값이므로 고정밀

MESH_POS = 3;
MESHCFG = struct('nLong',36,'nLat',45,'longMin',-180,'longMax',180, ...
                 'latMin',0,'latMax',90);
d_sub=1.3; r_OLED=1; wavelength_start=580; wavelength_end=590; n=10;

lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt,"ShowFileDialogBox",0);
List=ltml.LTDbList(lt,'lens_manager[1]','SIMULATIONS');
Key=ltml.LTListByName(lt,List,'ForwardAll');
ltml.LTDbSet(lt,Key,'MaxProgress',ray_nums);
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% ---- CPS + 하단반사 (objFcn 과 동일 코드경로) ----
load('nk_JH33.mat'); load('Photopic_400_800.mat'); load('CIE_1931.mat'); load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';
wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98; horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);
no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL 150];
CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,4,12.5,499,3,wavelength);
EQE_sub_CPS=CPS_result.EQE_sub;
TMF_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
Reflectance=(abs(TMF_p.r_p).^2+abs(TMF_s.r_s).^2)/2;
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
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY'); Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE'); Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d',count));

I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);

K=(wavelength_num-1)/n+1;
Power_output=zeros(1,wavelength_num);
Igrids=cell(1,K);
for wv=1:n:wavelength_num
    fileID=fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1,90,0,0,360,90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','DISK_SOURCE'); SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power',weight_factor(wv));
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','Spectral_region'); SRKey=ltml.LTListAtPos(lt,SRList,2);
    ltml.LTDbSet(lt,SRKey,'Spectral_Definition','Monochromatic');
    ltml.LTDbSet(lt,SRKey,'Single_Wavelength',wv+wavelength_start-1);
    List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER'); Key=ltml.LTListAtPos(lt,List,1);
    ltml.LTDbSet(lt,Key,'LoadFileName','C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt');
    ltml.LTBegin(lt); ltml.LTCmd(lt,'\V3D BeginAllSimulations'); ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH'); Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');
    [Igrids{(wv+n-1)/n}, thC, phC] = read_ff_mesh2d(ltml, lt, MESH_POS, MESHCFG);
end

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

Wacc=zeros(MESHCFG.nLat,MESHCFG.nLong);
sint=sind(thC(:));
for k=1:K
    Wk=Igrids{k}.*sint;  sk=sum(Wk(:));
    if sk>0, Wacc=Wacc+contrib(k)*(Wk/sk); end
end

save('bare_ref.mat','Wacc','EQE_total','thC','phC','dETL','dHTL','ray_nums');
fprintf('bare_ref.mat 저장: EQE_total(bare)=%.4f, 캐비티(dETL,dHTL)=(%g,%g)\n', EQE_total, dETL, dHTL);
