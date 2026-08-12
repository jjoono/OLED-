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
%% Parameter Definition
lb = [10, 10];   % lower bounds
ub = [150, 150]; % upper bounds
dim = numel(lb);

%% PSO Parameters
max_iter = 5;      % 최대 반복 횟수 (Simulation 비용 고려하여 조절 필요)
pop_size = 20;      % 입자(Particle)의 개수 (보통 20~50)
w = 0.7;            % 관성 가중치 (Inertia Weight)
c1 = 1.5;           % 개인 인지 계수 (Cognitive Coefficient)
c2 = 1.5;           % 사회적 인지 계수 (Social Coefficient)

% 속도 제한 (변수 범위의 20% 정도로 설정)
v_max = 0.2 * (ub - lb);
v_min = -v_max;
restart_interval=20;
eval_count=0;
%% Initialization
particles.position = zeros(pop_size, dim);
particles.velocity = zeros(pop_size, dim);
particles.cost = zeros(pop_size, 1);
particles.best_position = zeros(pop_size, dim);
particles.best_cost = zeros(pop_size, 1);

global_best.position = zeros(1, dim);
global_best.cost = -inf; % 우리는 EQE를 Maximize 하므로 초기값은 -무한대

% 초기 입자 생성
disp('Initializing Swarm...');
for i = 1:pop_size
    valid_pos = false;
    while ~valid_pos
        % 랜덤 위치 생성
        pos = lb + rand(1, dim) .* (ub - lb);
        % 기하학적 제약조건 확인
        if isValidPoints(pos)
            valid_pos = true;
            particles.position(i, :) = pos;
        end
    end

    % 초기 속도 0
    particles.velocity(i, :) = zeros(1, dim);

    eval_count = eval_count + 1;
    if mod(eval_count, restart_interval) == 0
        fprintf('\n[Init Refresh] 초기화 시뮬레이션 %d회 수행. LightTools를 재시작합니다...\n', eval_count);

        RenewLightTools(); % LightTools 껐다 켜고 ID 갱신
        
        % 재시작 후 옵션 재설정
        lt = ltloc.GetLTAPI(ID_swept); 
        ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
        pause(2); 
    end

    % 초기 평가
    val = objFcn_angularEQE(particles.position(i, :)).EQE_40_60;

    % 평가 결과 저장
    particles.cost(i) = val;
    particles.best_cost(i) = val;
    particles.best_position(i, :) = particles.position(i, :);

    % Global Best 업데이트
    if particles.cost(i) > global_best.cost
        global_best.cost = particles.cost(i);
        global_best.position = particles.position(i, :);
    end

    fprintf('Init Particle %d/%d: EQE = %.4f\n', i, pop_size, val);
end

% 결과 기록용 변수
history_best_cost = zeros(max_iter, 1);

%% PSO Main Loop
disp('Starting PSO Loop...');
eval_count=0;
RenewLightTools();
lt = ltloc.GetLTAPI(ID_swept); % swept entity
ltx= getltpointer(ID_swept);  % swept entity
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
for it = 1:max_iter
    fprintf('=== Iteration %d Start ===\n', it);
    for i = 1:pop_size

        % 1. 속도 업데이트
        % v = w*v + c1*r1*(pBest - x) + c2*r2*(gBest - x)
        r1 = rand(1, dim);
        r2 = rand(1, dim);

        particles.velocity(i, :) = w * particles.velocity(i, :) ...
            + c1 * r1 .* (particles.best_position(i, :) - particles.position(i, :)) ...
            + c2 * r2 .* (global_best.position - particles.position(i, :));

        % 속도 제한 적용
        particles.velocity(i, :) = max(particles.velocity(i, :), v_min);
        particles.velocity(i, :) = min(particles.velocity(i, :), v_max);

        % 2. 위치 업데이트
        new_position = particles.position(i, :) + particles.velocity(i, :);

        % 경계 처리 (Bounds Handling)
        new_position = max(new_position, lb);
        new_position = min(new_position, ub);

        % 3. 기하학적 제약조건 (isValidPoints) 체크
        % 만약 이동한 위치가 유효하지 않은 형상이라면?
        % 전략: 이번 이동을 취소하거나, 랜덤하게 다시 생성하거나, 해당 입자 평가 스킵
        % 여기서는 평가를 스킵하고 매우 낮은 점수를 부여하여 도태되게 함

        if ~isValidPoints(new_position)
            cost_new = -1; % Invalid shape penalty
        else
            % 4. 목적 함수 평가 (LightTools Simulation)
            particles.position(i, :) = new_position; % 유효할 때만 위치 확정
            eval_count=eval_count+1;
            if mod(eval_count,restart_interval)==0
                fprintf('\n[Refresh] 시뮬레이션 %d회 수행. LightTools를 재시작합니다...\n', eval_count);

                RenewLightTools(); % 여기서 LightTools 껐다 켜고 ID 갱신됨
                lt = ltloc.GetLTAPI(ID_swept); % swept entity
                ltx= getltpointer(ID_swept);  % swept entity
                ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
%                 ltml.LTSetOption(lt2, "ShowFileDialogBox", 0);
                pause(2); % 잠시 대기
            end
            cost_new = objFcn_angularEQE(new_position).EQE_40_60;
        end

        % 5. pBest 업데이트
        if cost_new > particles.best_cost(i)
            particles.best_cost(i) = cost_new;
            particles.best_position(i, :) = new_position;
        end

        % 6. gBest 업데이트
        if cost_new > global_best.cost
            global_best.cost = cost_new;
            global_best.position = new_position;
        end

    end % end of particles loop

    history_best_cost(it) = global_best.cost;

    fprintf('Iteration %d/%d : Best EQE = %.5f\n', it, max_iter, global_best.cost);

    % 간단한 그래프 업데이트
    figure(1);
    plot(1:it, history_best_cost(1:it), '-o', 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Best EQE_40_60');
    title(['PSO Optimization (Iter ' num2str(it) ')']);
    grid on; drawnow;

end

%% Save Results
result_pso = struct();
result_pso.best_position = global_best.position;
result_pso.best_cost = global_best.cost;
result_pso.history = history_best_cost;
result_pso.final_particles = particles;

save('PSO_Result_Trial1.mat', 'result_pso');
disp('Optimization Finished.');


%% Objective Function (기존 코드 그대로 유지)
function output = objFcn_angularEQE(point)
global ID_LT ID_swept ltml ltloc count
% Define segment length and other necessary parameters
lt = ltloc.GetLTAPI(ID_LT);  % lenssizeeffect
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
d_sub=1.3;
r_OLED=1;
x_pattern=25;
y_pattern=25;
Lensheight=0.01;
wavelength_start=450;
wavelength_end=750;
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
x2 = point(1);  x3 = point(2);  x4 = point(3);  x5 = point(4);  x6 = point(5);
y2 = point(6);  y3 = point(7);  y4 = point(8);  y5 = point(9);  y6 = point(10);
dETL = point(11); dHTL = point(12);

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

end

%% Spline Constraints Function (기존 코드 그대로 유지)
function TF = isValidPoints(X)
% 원본 코드의 isValidPoints 함수 내용을 그대로 유지하세요.
% X: N x 12 matrix (numeric)
numRows = size(X,1);
numPts  = 7;
TF = true(numRows,1);

for k = 1:numRows
    x = [0, X(k,1:5), 1];    % x2~x6
    y = [1, X(k,6:10), 0];   % y2~y6

    violates = false;

    % (1) Intersection
    for i = 1:numPts - 1
        for j = i + 2:numPts - 1
            if i == 1 && j == numPts - 1
                continue;
            end
            if checkIntersection([x(i), y(i)], [x(i+1), y(i+1)], ...
                    [x(j), y(j)], [x(j+1), y(j+1)])
                violates = true;
                break;
            end
        end
        if violates, break; end
    end

    % (2) Collinearity
    if ~violates
        for i = 1:numPts - 2
            if isCollinear([x(i), y(i)], [x(i+1), y(i+1)], [x(i+2), y(i+2)])
                violates = true;
                break;
            end
        end
    end

    % (3) Spacing
    if ~violates
        minD = 0.05; maxD = 1.0;
        d = hypot(diff(x), diff(y));
        if any(d < minD | d > maxD)
            violates = true;
        end
    end

    % (4) Angle
    if ~violates
        maxAng = 2 * pi / 3;
        for i = 2:numPts - 1
            v1 = [x(i), y(i)] - [x(i-1), y(i-1)];
            v2 = [x(i+1), y(i+1)] - [x(i), y(i)];
            ang = atan2(norm(cross([v1,0], [v2,0])), dot(v1, v2));
            if ang > maxAng
                violates = true;
                break;
            end
        end
    end

    TF(k) = ~violates;
end

% === Helper Functions ===
    function isCol = isCollinear(p1, p2, p3)
        area = 0.5 * det([p1 1; p2 1; p3 1]);
        isCol = abs(area) < 1e-5;
    end

    function intersects = checkIntersection(p1, p2, p3, p4)
        function o = orientation(p, q, r)
            o = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2));
        end
        o1 = orientation(p1, p2, p3);
        o2 = orientation(p1, p2, p4);
        o3 = orientation(p3, p4, p1);
        o4 = orientation(p3, p4, p2);
        intersects = (o1 * o2 < 0) && (o3 * o4 < 0);
    end
end


function RenewLightTools()
global ID_LT ID_swept ltml ltloc lt
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
model_file_path_swept = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\SweptEntity.2.lts';
model_file_path_LT = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\Lens_size_effect_for_PSO_bump_modified_v1.1.lts';
% =========================================================================

fprintf('--- Restarting LightTools ---\n');

% 1. 기존 LightTools 강제 종료
target_user = 'jhkim';
kill_cmd = sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user);
[~, ~] = system(kill_cmd);
pause(2);

% 2. 시스템 명령어로 .lts 파일 직접 실행
% 명령어 구조: "실행파일경로" "모델파일경로" &
% 끝에 '&'를 붙여야 MATLAB이 멈추지 않고 다음 줄로 넘어갑니다.
cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_swept);
status = system(cmd);
% 2. LightTools 재실행 및 연결
try
    % 새 인스턴스 생성
    % (LTAPI2를 새로 부르면 LightTools가 켜집니다)
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패. 라이선스나 설치 상태를 확인하세요.');
end

find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);

[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    % 정규식으로 PID(숫자) 추출
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');

    if ~isempty(tokens)
        % 결과가 여러 개일 경우(혹시 덜 꺼진게 있다면),
        % 보통 가장 마지막에 뜬 것(리스트의 끝 혹은 PID가 큰 것)이 새 프로세스일 확률이 높으나,
        % 여기서는 깨끗이 지우고 켰으므로 첫 번째 것을 가져옵니다.
        pid_str = tokens{1}{1};
        ID_swept = str2double(pid_str);

        fprintf('PID found for user %s: %d\n', target_user, ID_swept);
    else
        error('프로세스는 찾았으나 PID 추출 실패. 정규식 확인 필요.');
    end
else
    error('사용자 %s 로 실행된 LightTools(lt.exe)를 찾을 수 없습니다.', target_user);
end
cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_LT);

status = system(cmd);
% 2. LightTools 재실행 및 연결
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);

[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    % 정규식으로 PID(숫자) 추출
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');

    if ~isempty(tokens)
        % 결과가 여러 개일 경우(혹시 덜 꺼진게 있다면),
        % 보통 가장 마지막에 뜬 것(리스트의 끝 혹은 PID가 큰 것)이 새 프로세스일 확률이 높으나,
        % 여기서는 깨끗이 지우고 켰으므로 첫 번째 것을 가져옵니다.
        pid_str = tokens{3}{1};
        ID_LT = str2double(pid_str);
        
        fprintf('PID found for user %s: %d\n', target_user, ID_LT);
    else
        error('프로세스는 찾았으나 PID 추출 실패. 정규식 확인 필요.');
    end
else
    error('사용자 %s 로 실행된 LightTools(lt.exe)를 찾을 수 없습니다.', target_user);
end
pause(5);
end
