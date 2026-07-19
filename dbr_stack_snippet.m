% ============================================================
%  dbr_stack_snippet.m  -- objFcn 의 스택 정의를 Ag -> DBR 로 교체하는 조각
%
%  [목적] "무손실 DBR 캐비티가 효율-방향성 캡을 깨는가"를 사장님 실제 CPS 로 검증.
%  [사용] objFcn_regionPower(또는 bare_ceiling_CPS)의
%         no_bar / ne_bar / thickness 정의 3줄을 아래로 교체.
%         Ndbr(= DBR pair 수)를 sweep 하거나 변수로 추가해 비교.
%
%  [물리 근거 - TMM 확인됨] 같은 finesse(=같은 방향성)에서
%     Ag 50nm: outcoupling 24%  vs  DBR 4pair: 56%  (2.3x)
%     Ag 는 두꺼울수록 흡수로 효율 붕괴, DBR 은 무손실이라 유지.
%
%  [!] 주의 3가지
%   (1) DBR 은 비전도성 -> 실제 소자는 투명전극(얇은 Ag/ITO) + DBR 필요.
%       광학 상한 확인용이면 순수 DBR, 현실적이면 얇은 Ag 1층 + DBR 로.
%   (2) TiO2/SiO2 n,k 데이터가 material 구조체에 없으면 아래처럼 상수 n(무손실)
%       근사 사용. 정밀히 하려면 실제 분산 n,k 를 material 에 추가.
%   (3) 후면 Al 은 여전히 손실 -> 캡의 바닥(Al 흡수)은 남음. 완전 무손실은
%       후면도 DBR+투명전극이어야 하나 캐소드 전도성 문제로 별도 과제.
% ============================================================

%% ----- DBR 설계 파라미터 -----
lam0  = 600;          % DBR 중심파장 [nm] (오렌지 이미터 피크 근처로)
nHi   = 2.35;         % 고굴절 (TiO2). 실제 데이터 있으면 material.l_TiO2 로 교체
nLo   = 1.46;         % 저굴절 (SiO2). material.l_SiO2
Ndbr  = 4;            % DBR pair 수 (2~6 sweep 권장; 많을수록 고반사=고방향성=고손실민감)
useThinAgContact = true;   % true: 얇은 Ag 전극 1층 + DBR (현실적), false: 순수 DBR(광학상한)
tAgContact = 8;      % 얇은 Ag 전극 두께 [nm] (전도성 확보용, useThinAgContact=true 일 때)

%% ----- DBR 층 배열 만들기 (quarter-wave) -----
dHi = lam0/(4*nHi);   dLo = lam0/(4*nLo);          % 층 두께 [nm]
colHi = nHi*ones(401,1);  colLo = nLo*ones(401,1); % 무손실(k=0), 등방(no=ne)
dbr_cols = repmat([colHi colLo], 1, Ndbr);         % 401 x (2*Ndbr)
dbr_th   = repmat([dHi   dLo  ], 1, Ndbr);         % 1 x (2*Ndbr)

if useThinAgContact
    contact_col = material.l_Ag_McPeak;   contact_th = tAgContact;
else
    contact_col = zeros(401,0);           contact_th = [];   % 없음
end

%% ----- 스택 정의 (Ag 자리 -> [얇은Ag] + DBR) -----
% 기존:
%   no_bar=[air Al B3_o TCTA_B3_o TCTA_o TAPC_o  Ag         1.51];  th=[100 dETL 25 10 dHTL dAg]
% 교체:
%   no_bar=[air Al B3_o TCTA_B3_o TCTA_o TAPC_o  (Ag)+DBR   1.51];  th=[100 dETL 25 10 dHTL (tAg)+dbr]
no_bar = [ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO ...
          material.l_TCTA_o_JO material.l_TAPC_o_JO  contact_col dbr_cols  1.51*ones(401,1)];
ne_bar = [ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO ...
          material.l_TCTA_e_JO material.l_TAPC_e_JO  contact_col dbr_cols  1.51*ones(401,1)];
thickness = [100 dETL 25 10 dHTL  contact_th dbr_th];   % <-- Ag(dAg) 자리 교체

% 파장 슬라이싱은 기존 코드 그대로:
%   no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);  (아래 기존 줄 유지)
%   ne_bar=ne_bar(...);  layer_num=size(no_bar,2);  자동 반영됨.

%% [검증 절차]
%  1) 이 스택으로 CPS 돌려 EQE_sub, I_sub(θ) 확인:
%     - Ag 대비 I_sub(θ) 가 더 좁은 각도 링/피크로 모이는가? (방향성)
%     - EQE_sub(기판모드 효율) 가 Ag 얇은두께 대비 유지/향상되는가? (효율)
%  2) Ndbr = 2,3,4,5,6 sweep -> "방향성 vs 효율" 곡선이 Ag 곡선보다 위인지 확인.
%  3) 그다음 freeform 재지향까지 얹어 EQE_region 이 Ag 캡(~3%)을 넘는지.
