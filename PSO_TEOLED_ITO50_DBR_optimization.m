%% =====================================================================
%  Top-emitting OLED (TEOLED) 최적화
%  - Top 전극: thin Ag  ->  ITO 50 nm 로 교체
%  - ITO 위에 DBR(TiO2/SiO2) 을 N pair 만큼 적층
%  - 최적화 변수 : dETL, dHTL, dHigh(DBR high), dLow(DBR low), N(DBR pair 수)
%  최적화기      : particleswarm (PSO)
%
%  * N(정수) 은 particleswarm 이 정수를 직접 지원하지 않으므로
%    목적함수 내부에서 round() 로 처리한다.
%  * 필요한 사용자 함수/데이터 (원본과 동일):
%      CPS_for_Isub_v3, TMF_birefringence_whole_p, TMF_birefringence_whole_s
%      nk_JH_total.mat, Photopic_400_800.mat, CIE_1931.mat,
%      hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat
% =====================================================================

%% --- 1. PSO 기본 설정 ---------------------------------------------------
%          [ dETL   dHTL   dHigh   dLow    N   ]
nvars = 5;
lb    = [ 10     10     10      10      1   ];   % 하한
ub    = [ 250    250    150     200     10  ];   % 상한  (N: 1~10 pair)

%% --- 2. 병렬 처리 설정 -------------------------------------------------
numWorkers = 10;
if isempty(gcp('nocreate'))
    parpool('local', numWorkers);
end

options = optimoptions('particleswarm', ...
    'SwarmSize',      100, ...                 % 파티클 개수
    'MaxIterations',  50, ...                  % 최대 반복 (원본 5 -> 실사용 시 상향)
    'PlotFcn',        @pswplotbestf, ...
    'UseParallel',    true, ...
    ...  % HybridFcn(@fmincon) 은 정수변수 N 때문에 제거.
    ...  % 연속변수만 추가로 다듬고 싶으면 N 을 고정한 뒤 fmincon 별도 실행 권장.
    'OutputFcn',      @realtime_pso_monitor);

%% --- 3. PSO 실행 ------------------------------------------------------
fprintf('PSO 최적화를 시작합니다...\n');
[x_optimal, fval] = particleswarm(@evaluateEQE, nvars, lb, ub, options);

%% --- 4. 결과 출력 -----------------------------------------------------
fprintf('------------------------------------\n');
fprintf('최적화가 완료되었습니다.\n');
fprintf('  > 최적의 dETL      : %.4f nm\n', x_optimal(1));
fprintf('  > 최적의 dHTL      : %.4f nm\n', x_optimal(2));
fprintf('  > 최적의 dHigh(DBR): %.4f nm\n', x_optimal(3));
fprintf('  > 최적의 dLow (DBR): %.4f nm\n', x_optimal(4));
fprintf('  > 최적의 N (pair)  : %d\n',      round(x_optimal(5)));
fprintf('  > 최대 Totalpower  : %.4f\n',   -fval);
fprintf('------------------------------------\n');


%% =====================================================================
%  목적 함수 (Objective Function)
% =====================================================================
function output = evaluateEQE(point)
load('nk_JH_total.mat')
load('Photopic_400_800.mat');
load('CIE_1931.mat')
load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat')
BSDF=BSDF_MLA(:,:,10);

%% 1. 최적화 변수 -------------------------------------------------------
dETL  = point(1);
dHTL  = point(2);
dHigh = point(3);          % DBR high-index (TiO2) 두께
dLow  = point(4);          % DBR low-index  (SiO2) 두께
N     = round(point(5));   % DBR pair 수 (정수화)
if N < 1, N = 1; end

%% 2. 변수 정의 ---------------------------------------------------------
wavelength_start=450;
wavelength_end=750;
wavelength=(wavelength_start:5:wavelength_end).';

z0=12.5;
u_data_num=995;
max_u=3;
wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd;
emission_spectrum=emission_spectrum(51:5:351);
emission_spectrum=emission_spectrum/sum(emission_spectrum);
eta_rad=0.1;                         % IQE
horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);   % **렌즈(MLA) 사용 시 유지
high=real(material.l_TiO2_SJ_RTP);   % DBR high-index n (lossless)
low =real(material.l_SiO2_SJ_RTP);   % DBR low-index  n (lossless)

%% 3. 층 구조 정의 (top 전극 ITO 50nm + DBR N pair) --------------------
% ITO 굴절률: 실제 ITO nk 데이터가 있으면 아래 nk_ITO 를 교체할 것.
%   예) nk_ITO = material.l_ITO_xxx;
nk_ITO = 1.9*ones(401,1);            % (임시) ITO 굴절률
n21    = 2.1*ones(401,1);            % 버퍼/주입층
exit_n = 1.77*ones(401,1);           % 출광(반무한) 매질

% --- head : air ~ (버퍼) ~ ITO 까지 (원본과 동일, top Ag 만 ITO 로 교체) ---
%   col :  air   Ag(back)          BPhen             B3(ETL)          TCTA:B3(EML)          TCTA               TAPC(HTL)          buffer  ITO(top)
no_head=[ones(401,1) material.l_Ag_McPeak material.l_BPhen_CS material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO n21 nk_ITO];
ne_head=[ones(401,1) material.l_Ag_McPeak material.l_BPhen_CS material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO n21 nk_ITO];
th_head=[100 0 dETL 25 10 dHTL 5 50];   % 대응 두께 (ITO=50nm 고정)

% --- DBR : [high low] 를 N pair 반복 (유전체 -> ne=no) ---
dbr_nk = repmat([high low], 1, N);      % 401 x (2N)
dbr_th = repmat([dHigh dLow], 1, N);    % 1   x (2N)

% --- 전체 stack (마지막은 반무한 출광 매질) ---
no_bar    = [no_head dbr_nk exit_n];
ne_bar    = [ne_head dbr_nk exit_n];
thickness = [th_head dbr_th];           % 반무한 매질(air, exit)은 두께 미포함

EML_position=5;                         % TCTA:B3 (변경 없음)

% 파장 슬라이싱 (원본과 동일 방식)
no_bar=no_bar(wavelength_start-449:5:wavelength_end-449,:);
ne_bar=ne_bar(wavelength_start-449:5:wavelength_end-449,:);
layer_num=size(no_bar,2);
sin089=sind(0:89);
cos089=cosd(0:89);

%% 4. CPS 계산 ----------------------------------------------------------
CPS_result=CPS_for_Isub_v3(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_CPS=CPS_result.EQE_air;
EQE_sub_CPS=CPS_result.EQE_sub;

%% 5. TMF 계산 ----------------------------------------------------------
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
I_sub=(CPS_result.I_sub_s+CPS_result.I_sub_p)/2;

%% 6. Power 계산 --------------------------------------------------------
Psub_norm=I_sub.*sin089;
temp=repmat(sum(Psub_norm,2),1,90);
Psub_norm=Psub_norm./temp; clear temp;
modifiedMatrix = Psub_norm;
nanIndices = isnan(modifiedMatrix);
modifiedMatrix(nanIndices) = 0;
Psub_norm=modifiedMatrix; clear modifiedMatrix nanIndices;
ROLED=Reflectance;
P0=CPS_result.EQE_sub_matrix;
R_bot=repmat(ROLED,1,1,90);

R_step = 20;

BSDF_R = BSDF(180:-1:91, :);
BSDF_T_total = sum(BSDF(1:90, :));
BSDF_T=BSDF(1:1:90,:);

BSDF_CE = zeros(180, 90);
for i = 1:90
    for j = 1:180
        BSDF_CE(j,i) = BSDF(j,i)*sind(i-0.5)/sind(j-0.5);
    end
end
BSDF_T_CE = BSDF_CE(1:1:90, :);

R_matrix_1 = repmat(reshape(BSDF_R', 1, 90, 90), wavelength_num, 1, 1) .* R_bot;

for j=1:R_step
    for i = 1:wavelength_num
        temp_psub = Psub_norm(i, :);
        temp_rmat = reshape(R_matrix_1(i, :, :), 90, 90);
        temp_power(i,j) = temp_psub * temp_rmat^(j-1) * BSDF_T_total';
        temp_power_total(i,:,j) = CPS_result.P_sub(i,:) * temp_rmat^(j-1);
        temp_power_total_for_CIE(i,:,j) = I_sub(i,:) * temp_rmat^(j-1);
        temp_power_30(i,:,j)=temp_psub*temp_rmat^(j-1)*BSDF_T;
    end
    temp_power_30_mla(:,j)=sum(temp_power_30(:,1:30,j),2);
    power30_mla(j)=temp_power_30_mla(:,j)'*P0;
    Power(j)=temp_power(:,j)'*P0;
    P_MLA(:,:,j)=temp_power_total(:,:,j)*BSDF_T_CE';
    I_MLA_for_CIE(:,:,j)=temp_power_total_for_CIE(:,:,j)*BSDF_T_CE';
    mean_P_MLA = mean(P_MLA(:,2:5,j),2);
end
Totalpower=sum(Power);
Power_30_total=sum(power30_mla);

CIE_angle=30;
I_MLA_for_CIE_total=sum(I_MLA_for_CIE,3);
I_MLA_norm=I_MLA_for_CIE_total(:,1);
I_MLA_norms=I_MLA_norm/max(I_MLA_norm);
FWHM_MLA=sum(I_MLA_norms>=0.5);

A=I_MLA_for_CIE_total(:,1);
B=I_MLA_for_CIE_total;
x = (1:length(A))';
I_max = max(A);
Total_Area = trapz(x, A);
Integral_Width = Total_Area / I_max;
[~, spec_peak] = max(B, [], 1);
spec_shift=max(spec_peak)-min(spec_peak);

%% 7. 최종 출력 (PSO 는 최소화 -> Totalpower 최대화를 위해 음수 반환) ----
% 스펙트럼 안정성 제약을 걸고 싶으면 아래 조건문을 활성화하세요.
% if (Integral_Width < 15) && (spec_shift < 15)
    output = -Totalpower;
% else
%     output = 1e3;     % Hard penalty (제약 위반 영역 차단)
% end

fclose('all');
end


%% =====================================================================
%  PSO 실시간 모니터링 및 중간 저장 (Output Function)
% =====================================================================
function stop = realtime_pso_monitor(optimValues, state)
persistent iteration_log nvars_local
stop = false;

if strcmp(state, 'init')
    iteration_log = [];
    nvars_local = length(optimValues.bestx);
    if nvars_local > 0
        iteration_log = zeros(0, nvars_local + 1);
    end
    fprintf('--- 실시간 PSO 모니터링 및 중간 저장 시작 ---\n');
    fprintf(' Iter  |  Current Max Power (Totalpower) | 저장 상태\n');
    fprintf('----------------------------------------------------------\n');

elseif strcmp(state, 'iter')
    iter = optimValues.iteration;
    current_best_point = optimValues.bestx;
    current_best_fval  = optimValues.bestfval;

    fprintf(' %5d |      %.6f                 |', iter, -current_best_fval);

    if isempty(nvars_local) || nvars_local == 0
        nvars_local = length(current_best_point);
        iteration_log = zeros(0, nvars_local + 1);
    end
    new_row = [current_best_point, current_best_fval];
    iteration_log(end+1, :) = new_row;

    try
        save('data_in_progress_TEOLED_ITO50nm_DBR_optimization.mat', 'iteration_log');
        fprintf(' 저장됨\n');
    catch ME
        fprintf(' 저장 실패: %s\n', ME.message);
    end

elseif strcmp(state, 'done')
    fprintf('----------------------------------------------------------\n');
    fprintf('실시간 모니터링 및 중간 저장이 완료되었습니다.\n');
    fprintf('최종 iteration별 최고 기록이 "data_in_progress_TEOLED_ITO50nm_DBR_optimization.mat" 에 저장되었습니다.\n');
end
end
