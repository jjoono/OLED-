% ============================================================
%  dbr_stack_snippet.m  -- 스택을 ITO(baseline) / Ag / DBR 로 바꿔 비교하는 조각
%
%  [권장 순서] STACK_MODE 를 'ITO' -> 'Ag' -> 'DBR' 로 바꿔가며 각각 CPS 돌려
%    EQE_sub(효율) 와 I_sub(θ)(방향성) 를 비교. 사다리:
%      ITO  : 금속미러 없음 = 약한 캐비티 -> 효율 높음, 방향성 ~0  (baseline 앵커)
%      Ag   : 손실 금속 캐비티 -> 방향성 조금, 두꺼울수록 효율 붕괴
%      DBR  : 무손실 유전체 캐비티(+ITO 전극) -> 방향성 유지하며 효율 유지 가능?
%
%  [사용] objFcn 의 no_bar/ne_bar/thickness 정의를 아래로 교체.
%    (dETL,dHTL 은 그대로 변수. Ag 모드의 dAg 도 그대로. DBR 은 Ndbr 로 조절.)
%
%  [!] ITO/DBR 은 비자성·등방 가정(no=ne). TiO2/SiO2 실제 n,k 있으면 교체.
%      후면 Al 은 여전히 손실이라 캡의 바닥은 남음(별도 과제).
% ============================================================

STACK_MODE = 'ITO';   % 'ITO' | 'Ag' | 'DBR'   <-- 이것만 바꿔가며 비교

% 공통 유기층 (air | Al | B3PyMPM | TCTA:B3(EML) | TCTA | TAPC | [상단] | substrate)
airc = ones(401,1);  subc = 1.51*ones(401,1);
o_common = [material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO];
e_common = [material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO];
th_common = [dETL 25 10 dHTL];        % B3PyMPM, EML, TCTA, TAPC 두께

switch STACK_MODE
case 'ITO'   % ===== baseline: ITO 전극만 (금속 상단미러 없음) =====
    dITO = 150;                        % ITO 두께 [nm] (고정; 원 코드값)
    top_o = material.l_ITO_SNU_temp;   top_e = material.l_ITO_SNU_temp;  top_th = dITO;

case 'Ag'    % ===== 현재: 얇은 Ag 캐비티 =====
    top_o = material.l_Ag_McPeak;      top_e = material.l_Ag_McPeak;     top_th = dAg;   % dAg 변수

case 'DBR'   % ===== 무손실 DBR (+ ITO 전극) =====
    lam0=600; nHi=2.35; nLo=1.46; Ndbr=4;   % Ndbr 2~6 sweep
    dHi=lam0/(4*nHi); dLo=lam0/(4*nLo);
    dITO=30;                            % 얇은 ITO 전극(전도성 확보), DBR 앞
    ito_o=material.l_ITO_SNU_temp;
    dbr_o = repmat([nHi*ones(401,1) nLo*ones(401,1)],1,Ndbr);
    dbr_th= repmat([dHi dLo],1,Ndbr);
    top_o = [ito_o dbr_o];  top_e = [ito_o dbr_o];  top_th = [dITO dbr_th];
end

no_bar    = [airc material.l_Al_JO o_common top_o subc];
ne_bar    = [airc material.l_Al_JO e_common top_e subc];
thickness = [100 th_common top_th];        % Al=100 + 유기 + 상단
% (이후 기존 코드의 파장 슬라이싱 / layer_num=size(no_bar,2) / CPS 호출 그대로)

%% [비교 지표]
%   - EQE_sub (기판모드 효율): ITO≈최고, Ag=중(흡수), DBR=? (무손실이라 유지 기대)
%   - I_sub(θ) 각도폭: ITO=넓음, Ag=약간좁음, DBR(Ndbr↑)=좁은 링 기대
%   -> "방향성 vs 효율" 평면에 세 점 찍으면, DBR 이 Ag 곡선보다 위(우상단)인지 확정.
