% ============================================================
%  plot_fig2a_selectivity.m
%
%  Fig. 2a 재작도. 기존 막대그래프(Lambertian x 0.56 기준)는 두 가지를 섞어
%  보여줘서 오해를 부른다: (i) 총 추출 차이, (ii) 배광(선택성) 차이.
%  여기서는 둘을 분리한다.
%
%  (a) 선택성 비교
%      회색 띠 = 총 EQE 상위 N_TOP 설계들의 선택성 10~90% 범위
%                = "효율이 높은 설계라면 어떤 목적으로 찾았든 갖게 되는 배광"
%                  (= 무조향 기준. 대조군을 따로 돌리지 않아도 되는 이유)
%      빨간 점  = 그 구간을 전용 목적함수로 최적화한 설계의 선택성
%      점선     = Lambertian 기준값
%  (b) 이득 분해
%      선택성 이득 / 총량 비 / 순이득 을 구간별로 나란히
%      순이득 = 전용 band EQE / (자연 선택성 x 최대 총 EQE)
%
%  입력: opt_4band_result.mat  (opt_4band_freeform.m 산출물)
%        EVAL_LOG = [x(1:13) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w]
%  출력: fig2a_selectivity.png
% ============================================================
clear;

RESULT_FILE = 'opt_4band_result.mat';
N_TOP       = 20;      % 자연 선택성을 뽑을 상위 총EQE 설계 수
PCT_LO      = 10;      % 띠의 하단/상단 백분위
PCT_HI      = 90;

D = load(RESULT_FILE);
EVAL_LOG   = D.EVAL_LOG;
band_eqe_hi= D.band_eqe_hi;
band_tot_hi= D.band_tot_hi;
if isfield(D,'BAND_NAMES'), BAND_NAMES = D.BAND_NAMES; else
    BAND_NAMES = {'0-20 deg','20-40 deg','40-60 deg','60-80 deg'}; end

S_LAMB = [0.117 0.296 0.337 0.220];
nB     = 4;

%% ===== 자연 선택성 (무조향 기준) =====
T = EVAL_LOG(:,14);
B = EVAL_LOG(:,15:18);
ok = isfinite(T) & T > 0 & all(isfinite(B),2);
T = T(ok);  B = B(ok,:);

[~, ord] = sort(T, 'descend');
top   = ord(1:min(N_TOP, numel(ord)));
Sdist = B(top,:) ./ T(top);            % 상위 설계들의 선택성 분포
S_lo  = prctile(Sdist, PCT_LO, 1);
S_hi  = prctile(Sdist, PCT_HI, 1);
S_med = median(Sdist, 1);
E_max = max(T);

%% ===== 전용 설계 =====
S_win = band_eqe_hi(1:nB)' ./ band_tot_hi(1:nB)';

gain_sel = S_win ./ S_med;                          % 선택성 이득
gain_tot = band_tot_hi(1:nB)' / E_max;              % 총량 비
gain_net = band_eqe_hi(1:nB)' ./ (S_med * E_max);   % 순이득

%% ===== 콘솔 요약 =====
fprintf('\n자연 선택성 (상위 %d개, %d~%d%%):\n', N_TOP, PCT_LO, PCT_HI);
for k = 1:nB
    fprintf('  %-10s %.3f - %.3f  (median %.3f, Lambertian %.3f)\n', ...
        BAND_NAMES{k}, S_lo(k), S_hi(k), S_med(k), S_LAMB(k));
end
fprintf('\n%-10s %8s %10s %9s %9s\n','band','S_win','선택성이득','총량비','순이득');
for k = 1:nB
    fprintf('%-10s %8.3f %10.3f %9.3f %9.3f\n', ...
        BAND_NAMES{k}, S_win(k), gain_sel(k), gain_tot(k), gain_net(k));
end

%% ===== 그림 =====
figure('Name','Fig 2a - selectivity vs natural band','Color','w', ...
       'Position',[100 100 980 400]);

% --- (a) 선택성 ---
subplot(1,2,1); hold on; box on;
xx = 1:nB;
for k = 1:nB
    % 자연 범위 띠
    patch([k-0.28 k+0.28 k+0.28 k-0.28], [S_lo(k) S_lo(k) S_hi(k) S_hi(k)], ...
          [0.80 0.80 0.82], 'EdgeColor', [0.55 0.55 0.58], 'FaceAlpha', 0.9);
    % Lambertian
    plot([k-0.34 k+0.34], [S_LAMB(k) S_LAMB(k)], 'k--', 'LineWidth', 1.1);
end
hW = plot(xx, S_win, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.85 0.20 0.20], ...
          'MarkerEdgeColor', 'k', 'LineStyle', 'none');
set(gca,'XTick',xx,'XTickLabel',BAND_NAMES,'XLim',[0.5 nB+0.5]);
ylabel('selectivity  S_j = EQE_{band} / EQE_{total}');
title(sprintf('(a) dedicated design vs natural spread (top %d)', N_TOP));
legend([hW], {'band-dedicated optimum'}, 'Location','northwest','FontSize',9);
text(0.55, max(S_hi)*0.98, sprintf('grey band: %d–%d%% of top-%d designs\ndashed: Lambertian', ...
     PCT_LO, PCT_HI, N_TOP), 'FontSize', 8, 'VerticalAlignment','top');
grid on;

% --- (b) 이득 분해 ---
subplot(1,2,2); hold on; box on;
bar(xx, [gain_sel(:), gain_tot(:), gain_net(:)]);
plot([0.5 nB+0.5], [1 1], 'k-', 'LineWidth', 1);
set(gca,'XTick',xx,'XTickLabel',BAND_NAMES,'XLim',[0.5 nB+0.5]);
ylabel('ratio');
title('(b) selectivity gain x total ratio = net gain');
legend({'selectivity gain  S_{win}/S_{nat}', 'total ratio  E_{win}/E_{max}', ...
        'net gain'}, 'Location','northoutside','Orientation','horizontal','FontSize',8);
grid on;

saveas(gcf, 'fig2a_selectivity.png');
fprintf('\nsaved -> fig2a_selectivity.png\n');
