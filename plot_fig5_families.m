% ============================================================
%  plot_fig5_families.m
%
%  Fig. 5 — MLA family 일반성 검사. family 당 3패널.
%    (a) 총 EQE 대 40-60도 band EQE 산점도 + 선형 fit  -> 준선형 붕괴
%    (b) 자연 배광(상위 N_TOP 설계 중앙값)을 볼록 기준·Lambertian 과 비교
%    (c) band 별 선택성-효율 상관 R 을 볼록 기준선과 나란히
%
%  [기준값] 볼록(주기 freeform) 계열, patch 25:
%      자연 배광 0.096 / 0.278 / 0.360 / 0.237
%      상관 R    +0.60 / +0.70 / +0.05 / -0.70
%
%  입력: 아래 FAMILIES 에 적힌 .mat 중 **존재하는 것만** 그린다.
%    - opt_4band_inverted_result.mat  (오목: opt_4band_inverted.m 산출)
%    - stress_random_result.mat       (무작위: stress_random_mla.m 산출)
%    - opt_4band_result.mat           (볼록 기준, 있으면 첫 행에 함께 표시)
%  출력: fig5_families.png
% ============================================================
clear;

S_CONVEX = [0.096 0.278 0.360 0.237];
S_LAMB   = [0.117 0.296 0.337 0.220];
R_CONVEX = [0.60 0.70 0.05 -0.70];
BANDS    = {'0-20','20-40','40-60','60-80'};
N_TOP    = 20;

FAMILIES = { ...
  'opt_4band_result.mat',          'convex (reference)'; ...
  'opt_4band_inverted_result.mat', 'inverted (concave)'; ...
  'stress_random_result.mat',      'randomly assembled' };

have = false(size(FAMILIES,1),1);
for k = 1:size(FAMILIES,1)
    have(k) = exist(FAMILIES{k,1}, 'file') == 2;
    if ~have(k)
        fprintf('[skip] %s 없음\n', FAMILIES{k,1});
    end
end
idx = find(have);
if isempty(idx), error('그릴 결과 파일이 하나도 없다.'); end

nF = numel(idx);
figure('Name','Fig 5 - generality across MLA families','Color','w', ...
       'Position',[60 60 1150 340*nF]);

fprintf('\n%-22s %8s | %s\n','family','E*','자연 배광 (0-20 / 20-40 / 40-60 / 60-80)');
for r = 1:nF
    k = idx(r);
    D = load(FAMILIES{k,1});
    L = D.EVAL_LOG;
    T = L(:,14);  B = L(:,15:18);
    ok = isfinite(T) & T > 0 & all(isfinite(B),2);
    T = T(ok);  B = B(ok,:);

    [~,ord] = sort(T,'descend');
    top  = ord(1:min(N_TOP,numel(ord)));
    Snat = median(B(top,:)./T(top), 1);
    Es   = max(T);

    Rsel = nan(1,4);
    for b = 1:4
        c = corrcoef(T, B(:,b)./T);  Rsel(b) = c(1,2);
    end

    fprintf('%-22s %8.4f | %.3f %.3f %.3f %.3f\n', ...
            FAMILIES{k,2}, Es, Snat);

    % --- (a) 붕괴 산점도 ---
    subplot(nF,3,(r-1)*3+1); hold on; box on;
    scatter(T, B(:,3), 12, [0.55 0.6 0.7], 'filled', 'MarkerFaceAlpha',0.55);
    pf = polyfit(T, B(:,3), 1);
    xx = linspace(min(T), max(T), 50);
    plot(xx, polyval(pf,xx), 'k-', 'LineWidth', 1.3);
    cc = corrcoef(T, B(:,3));
    xlabel('EQE_{total}'); ylabel('EQE_{40-60}');
    title(sprintf('(%c1) %s   R^2=%.3f', 'a'+r-1, FAMILIES{k,2}, cc(1,2)^2));
    grid on;

    % --- (b) 자연 배광 비교 ---
    subplot(nF,3,(r-1)*3+2); hold on; box on;
    bar(1:4, [Snat(:), S_CONVEX(:), S_LAMB(:)]);
    set(gca,'XTick',1:4,'XTickLabel',BANDS);
    ylabel('natural selectivity S_j');
    title(sprintf('(%c2) composition', 'a'+r-1));
    if r == 1
        legend({'this family','convex ref','Lambertian'}, ...
               'Location','northwest','FontSize',7);
    end
    grid on;

    % --- (c) 선택성-효율 상관 ---
    subplot(nF,3,(r-1)*3+3); hold on; box on;
    bar(1:4, [Rsel(:), R_CONVEX(:)]);
    plot([0.5 4.5],[0 0],'k-');
    set(gca,'XTick',1:4,'XTickLabel',BANDS,'XLim',[0.5 4.5]);
    ylabel('R( EQE_{total}, S_j )');
    title(sprintf('(%c3) drift', 'a'+r-1));
    if r == 1
        legend({'this family','convex ref'},'Location','southwest','FontSize',7);
    end
    grid on;
end

saveas(gcf, 'fig5_families.png');
fprintf('\nsaved -> fig5_families.png\n');
