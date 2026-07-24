%% =====================================================================
%  Top-emitting OLED (TEOLED) 최적화 - Bayesian Optimization (bayesopt)
%  - Top 전극 : thin Ag -> ITO 50 nm
%  - ITO 위에 DBR(TiO2/SiO2) N pair 적층
%  - 최적화 변수 : dETL, dHTL, dHigh, dLow, N(DBR pair, 정수)
%
%  필요한 사용자 함수/데이터 (원본과 동일):
%    CPS_for_Isub_v3, TMF_birefringence_whole_p, TMF_birefringence_whole_s
%    nk_JH_total.mat, hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat
% =====================================================================

%% --- 1. 최적화 변수 정의 (N 은 정수형 -> bayesopt 네이티브 지원) -------
vars = [
    optimizableVariable('dETL',  [10  250])
    optimizableVariable('dHTL',  [10  250])
    optimizableVariable('dHigh', [10  150])
    optimizableVariable('dLow',  [10  200])
    optimizableVariable('N',     [1   10],  'Type','integer')
];

%% --- 2. 병렬 처리 ------------------------------------------------------
numWorkers = 10;
if isempty(gcp('nocreate'))
    parpool('local', numWorkers);
end

%% --- 3. Bayesian Optimization 실행 (빠른 수렴 옵션) -------------------
fprintf('Bayesian 최적화를 시작합니다...\n');
results = bayesopt(@evaluateEQE, vars, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ... % 탐색+수렴 균형, local minima 회피
    'MaxObjectiveEvaluations', 120, ...     % 총 평가 횟수
    'NumSeedPoints',           20, ...      % 초기 랜덤 샘플 (GP 초기화)
    'ExplorationRatio',        0.5, ...      % EI-plus 탐색 강도
    'GPActiveSetSize',         300, ...      % GP 학습 포인트 상한 (속도)
    'IsObjectiveDeterministic',true, ...     % 결정론적 시뮬 -> GP 적합/수렴 향상
    'UseParallel',             true, ...     % 병렬 평가로 가속
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, ...
    'OutputFcn',               @saveBOProgress, ...
    'Verbose',                 1);

%% --- 4. 결과 출력 -----------------------------------------------------
xbest = results.XAtMinObjective;    % 관측 최적점
fbest = results.MinObjective;
fprintf('------------------------------------\n');
fprintf('최적화가 완료되었습니다.\n');
fprintf('  > 최적의 dETL      : %.4f nm\n', xbest.dETL);
fprintf('  > 최적의 dHTL      : %.4f nm\n', xbest.dHTL);
fprintf('  > 최적의 dHigh(DBR): %.4f nm\n', xbest.dHigh);
fprintf('  > 최적의 dLow (DBR): %.4f nm\n', xbest.dLow);
fprintf('  > 최적의 N (pair)  : %d\n',      xbest.N);
fprintf('  > 최대 Totalpower  : %.4f\n',   -fbest);
fprintf('------------------------------------\n');


%% =====================================================================
%  목적 함수 : Totalpower 최대화 -> -Totalpower 최소화
% =====================================================================
function objective = evaluateEQE(x)
load('nk_JH_total.mat')
load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat')
BSDF = BSDF_MLA(:,:,10);

%% 변수
dETL  = x.dETL;
dHTL  = x.dHTL;
dHigh = x.dHigh;
dLow  = x.dLow;
N     = double(x.N);

%% 시뮬 파라미터
wavelength_start = 450;
wavelength_end   = 750;
wavelength = (wavelength_start:5:wavelength_end).';
wavelength_num = length(wavelength);

z0 = 12.5;  u_data_num = 995;  max_u = 3;
emission_spectrum = spectrum.l_I_Irdmppyph2tmd;
emission_spectrum = emission_spectrum(51:5:351);
emission_spectrum = emission_spectrum / sum(emission_spectrum);
eta_rad = 0.1;
horizontal_dipole_ratio = 0.865;
bottom_air_refractive_index = ones(wavelength_num,1);   % MLA(렌즈) 사용 시 유지
high = real(material.l_TiO2_SJ_RTP);   % DBR high (TiO2)
low  = real(material.l_SiO2_SJ_RTP);   % DBR low  (SiO2)

%% 층 구조 : top 전극 ITO 50nm + DBR N pair
nk_ITO = 1.9*ones(401,1);   % (임시) 실제 ITO nk 로 교체 권장: material.l_ITO_xxx
n21    = 2.1*ones(401,1);
exit_n = 1.77*ones(401,1);

no_head = [ones(401,1) material.l_Ag_McPeak material.l_BPhen_CS material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO n21 nk_ITO];
ne_head = [ones(401,1) material.l_Ag_McPeak material.l_BPhen_CS material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO n21 nk_ITO];
th_head = [100 0 dETL 25 10 dHTL 5 50];         % ITO = 50nm 고정

dbr_nk = repmat([high low], 1, N);              % 401 x 2N (유전체: ne=no)
dbr_th = repmat([dHigh dLow], 1, N);            % 1   x 2N

no_bar    = [no_head dbr_nk exit_n];
ne_bar    = [ne_head dbr_nk exit_n];
thickness = [th_head dbr_th];
EML_position = 5;

no_bar = no_bar(wavelength_start-449:5:wavelength_end-449, :);
ne_bar = ne_bar(wavelength_start-449:5:wavelength_end-449, :);
layer_num = size(no_bar,2);
sin089 = sind(0:89);  cos089 = cosd(0:89);

%% CPS
CPS_result = CPS_for_Isub_v3(no_bar,ne_bar,thickness,emission_spectrum,eta_rad, ...
    horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);

%% TMF (Totalpower 계산에 필요한 Reflectance 만)
TMF_p = TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s = TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
Reflectance = (abs(TMF_p.r_p).^2 + abs(TMF_s.r_s).^2) / 2;
I_sub = (CPS_result.I_sub_s + CPS_result.I_sub_p) / 2;

%% Power 계산
Psub_norm = I_sub .* sin089;
Psub_norm = Psub_norm ./ repmat(sum(Psub_norm,2),1,90);
Psub_norm(isnan(Psub_norm)) = 0;

P0    = CPS_result.EQE_sub_matrix;
R_bot = repmat(Reflectance,1,1,90);
R_step = 20;

BSDF_R       = BSDF(180:-1:91, :);
BSDF_T_total = sum(BSDF(1:90, :));
R_matrix_1   = repmat(reshape(BSDF_R', 1, 90, 90), wavelength_num, 1, 1) .* R_bot;

temp_power = zeros(wavelength_num, R_step);
Power      = zeros(1, R_step);
for j = 1:R_step
    for i = 1:wavelength_num
        temp_rmat = reshape(R_matrix_1(i,:,:), 90, 90);
        temp_power(i,j) = Psub_norm(i,:) * temp_rmat^(j-1) * BSDF_T_total';
    end
    Power(j) = temp_power(:,j)' * P0;
end
Totalpower = sum(Power);

%% 목적값 (최대화 -> 음수 최소화)
objective = -Totalpower;

fclose('all');
end


%% =====================================================================
%  진행상황 저장 (bayesopt Output Function)
% =====================================================================
function stop = saveBOProgress(results, state)
stop = false;
switch state
    case 'initial'
        fprintf('--- BO 모니터링/중간 저장 시작 ---\n');
    case 'iteration'
        try
            save('data_in_progress_TEOLED_ITO50nm_DBR_BO.mat', 'results');
        catch ME
            fprintf(' 저장 실패: %s\n', ME.message);
        end
    case 'done'
        save('data_final_TEOLED_ITO50nm_DBR_BO.mat', 'results');
        fprintf('--- BO 완료, 최종 결과 저장됨 ---\n');
end
end
