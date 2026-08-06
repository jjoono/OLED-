% ============================================================
%  make_fig3_fom_map.m   --  Fig 3 (FoM 지도) 생성 및 검증
%
%  [입력]  pareto_front_result.mat  (pareto_front_freeform.m 이 저장)
%     EVAL_LOG = [ x(1:13) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w ]
%
%  [검증할 것]  마스터 관계식은 목표 입체각 Omega 만 바꾸면 모든 FoM 에 적용된다.
%    각 각도구간의 선택성이 Lambertian 값에 고정되어야 한다:
%        S_[th1,th2] = sin^2(th2) - sin^2(th1)
%      0-20  : 0.117      20-40 : 0.296
%      40-60 : 0.337      60-80 : 0.220
%    -> 시뮬 한 번으로 네 개의 독립 예측을 동시에 검증한다.
%       네 구간이 모두 예측값에 고정되면, 선택성이 40-60 에서만 우연히 맞은 것이
%       아니라 '설계 불변량'임이 확정된다 (원고 2.3 / 4.2절).
%
%  [출력]  fig3_fom_map.png, fig3_summary.txt
% ============================================================
clear;
D = load('pareto_front_result.mat');
LOG = D.EVAL_LOG;
nvar = numel(D.lb);

Et  = LOG(:, nvar+1);                 % EQE_total
B   = LOG(:, nvar+2 : nvar+5);        % [0-20, 20-40, 40-60, 60-80]
Ph  = LOG(:, nvar+6);                 % phase (1=random, 2=opt, 3=고정밀)

edges  = [0 20; 20 40; 40 60; 60 80];
names  = {'0–20°','20–40°','40–60°','60–80°'};
% Lambertian(에텐듀 보존) 예측 선택성
S_pred = sind(edges(:,2)).^2 - sind(edges(:,1)).^2;

ok = isfinite(Et) & Et > 0.05 & all(isfinite(B),2);   % 너무 낮은 설계는 노이즈 지배
fprintf('유효 설계 %d개 (전체 %d개 중)\n\n', sum(ok), size(LOG,1));

%% ---- 구간별 선택성 통계 ----
fid = fopen('fig3_summary.txt','w');
hdr = sprintf('%-8s | %8s | %14s | %8s | %8s\n', ...
    'band','예측 S','실측 S (mean±std)','편차','상관');
fprintf('%s', hdr);  fprintf(fid,'%s',hdr);
fprintf('%s\n', repmat('-',1,60));  fprintf(fid,'%s\n', repmat('-',1,60));

S_meas = nan(4,1);  S_std = nan(4,1);  R_corr = nan(4,1);
for b = 1:4
    s = B(ok,b) ./ Et(ok);
    S_meas(b) = mean(s,'omitnan');
    S_std(b)  = std(s,'omitnan');
    % EQE_total 과의 상관: 불변량이면 상관이 0에 가까워야 한다
    c = corrcoef(Et(ok), s);   R_corr(b) = c(1,2);
    line = sprintf('%-8s | %8.3f | %6.3f ± %.3f | %+7.1f%% | %+6.2f\n', ...
        names{b}, S_pred(b), S_meas(b), S_std(b), ...
        100*(S_meas(b)-S_pred(b))/S_pred(b), R_corr(b));
    fprintf('%s', line);  fprintf(fid,'%s',line);
end
fprintf(fid,'\n* 상관계수가 0 근처면 선택성이 EQE_total 과 무관 = 설계 불변량\n');
fclose(fid);
fprintf('\n* 상관계수가 0 근처면 선택성이 EQE_total 과 무관 = 설계 불변량\n');
fprintf('saved -> fig3_summary.txt\n');

%% ---- 그림 ----
figure('Name','Fig 3 — FoM map','Color','w','Position',[80 80 1200 780]);
cols = lines(4);

% (a) 마스터 관계식: 목표 입체각별 단일통과 상한 (해석적, 시뮬 불필요)
subplot(2,2,1);
n_s = 1.51;
th = linspace(0,90,200);
bound_cone = sind(th).^2 / n_s^2;                 % 정면 원뿔 [0, th]
plot(th, 100*bound_cone, 'k-', 'LineWidth', 2); hold on;
for b = 1:4
    bb = (sind(edges(b,2))^2 - sind(edges(b,1))^2)/n_s^2;
    plot(mean(edges(b,:)), 100*bb, 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', cols(b,:), 'MarkerEdgeColor','k');
end
yline(100/n_s^2,'--','전 반구 (1/n_s^2)','LabelHorizontalAlignment','left');
xlabel('\theta (deg)'); ylabel('단일통과 상한  P_\Omega/P_{sub}  (%)');
title('(a) 하나의 식, 여러 목표 입체각'); grid on;
legend({'정면 원뿔 [0,\theta]','각 구간 밴드'},'Location','northwest','FontSize',8);

% (b) 구간별 선택성 vs EQE_total  -> 불변성의 핵심 증거
subplot(2,2,2);
for b = 1:4
    scatter(Et(ok), B(ok,b)./Et(ok), 10, cols(b,:), 'filled', ...
        'MarkerFaceAlpha', 0.35); hold on;
    yline(S_pred(b), '--', 'Color', cols(b,:), 'LineWidth', 1.6);
end
xlabel('EQE_{total}'); ylabel('선택성  S = EQE_{band}/EQE_{total}');
title('(b) 선택성은 설계 불변량 (점선 = Lambertian 예측)'); grid on;
legend(names,'Location','east','FontSize',8);

% (c) 예측 vs 실측 (parity)
subplot(2,2,3);
errorbar(S_pred, S_meas, S_std, 'o', 'MarkerSize', 9, 'LineWidth', 1.4, ...
    'MarkerFaceColor',[.2 .4 .7]); hold on;
lim = [0 0.42];
plot(lim, lim, 'k--', 'LineWidth', 1.2);
for b = 1:4
    text(S_pred(b)+0.012, S_meas(b), names{b}, 'FontSize', 8);
end
xlim(lim); ylim(lim); axis square; grid on;
xlabel('예측  sin^2\theta_2 - sin^2\theta_1'); ylabel('실측 선택성');
title('(c) 네 구간 독립 검증');

% (d) 가능 / 금지 지도
subplot(2,2,4);
axis([0 1 0 1]); axis off; hold on;
rectangle('Position',[0.05 0.55 0.42 0.35],'FaceColor',[.85 .95 .90],'EdgeColor',[.2 .5 .4]);
text(0.26,0.86,'REACHABLE','HorizontalAlignment','center','FontWeight','bold','Color',[.15 .45 .35]);
text(0.26,0.74,sprintf('총 추출\n정면 증강\n밴드 집광 (S 고정)'),'HorizontalAlignment','center','FontSize',8);
rectangle('Position',[0.53 0.55 0.42 0.35],'FaceColor',[.98 .93 .82],'EdgeColor',[.7 .55 .2]);
text(0.74,0.86,'BOUNDED','HorizontalAlignment','center','FontWeight','bold','Color',[.6 .45 .15]);
text(0.74,0.74,sprintf('재활용 축적\n(각도 선별 필요,\n손실이 제한)'),'HorizontalAlignment','center','FontSize',8);
rectangle('Position',[0.05 0.10 0.90 0.35],'FaceColor',[.98 .88 .88],'EdgeColor',[.7 .3 .3]);
text(0.50,0.40,'FORBIDDEN','HorizontalAlignment','center','FontWeight','bold','Color',[.6 .2 .2]);
text(0.50,0.26,sprintf(['순 방위각 조향  (\\oint\\nablaz dA = 0)\n' ...
    '중앙이 어두운 듀얼뷰\n타일링 어레이의 각도 압축  (A_{out}/A_{src} = 1)']), ...
    'HorizontalAlignment','center','FontSize',8);
title('(d) 가능 / 금지 지도');

saveas(gcf,'fig3_fom_map.png');
fprintf('saved -> fig3_fom_map.png\n');
