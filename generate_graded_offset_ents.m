%% generate_graded_offset_ents.m
% 위치별로 다르게 적용할 "편심 grading" 렌즐릿 .ent 세트 생성.
% 큰 렌즈를 프레넬화한 것처럼, 배열 중심->가장자리로 갈수록 렌즐릿의 정점(편심)을
% 바깥으로 민다. 각 zone(반경 링)마다 하나의 .ent -> LightTools 에서 zone 별로
% 다른 unit-cell 로 물리거나(방법 A), 개별 엔티티로 배치(방법 B)한다.
%
% [원리] 발광원이 기판(gap 아래)에 있고, 어레이 위치 (X,Y) 의 렌즐릿은 소스의
% chief ray 를 atan(R/gap) 각도로 받는다(R=중심거리). 이를 목표방향으로 펴려면
% 렌즐릿마다 편심을 R 에 비례해 키워야 한다(=집단 collimation, 에텐듀 harvest).
% 여기선 정규화 편심 dx 를 zone 별로 선형 증가시켜 그 grading 을 이산 구현.
%
% [제작성] 편심돔은 여전히 raised-cosine + rim window(테두리 0) 로 타일링 가능.
% dx 가 커지면 한쪽이 셀을 벗어나므로 rim window 가 잘라준다(자연스런 클리핑).

templatePath = 'freeform_template_v2.ent';
Ra=1.2139; Rap=1.0; n=141; tbase=0.30; h=0.72;
Nzone = 6;                       % 반경 zone 개수 (= .ent 개수)
dx_max = 0.75;                   % 가장 바깥 zone 의 정규화 편심 (Rap 기준)
dxs = linspace(0, dx_max, Nzone);% zone 별 편심 (선형 grading; 필요시 R^1 등으로 교체)

g=linspace(-Ra,Ra,n); [X,Y]=meshgrid(g,g); r=hypot(X,Y);
fprintf('graded 편심 세트: %d zone, dx = %s\n', Nzone, mat2str(round(dxs,3)));
for k=1:Nzone
    dx=dxs(k);
    H = offset_dome(X,Y,Rap,h,dx);
    [~,im]=max(H(:)); apex=X(im);
    outPath = sprintf('ff_grad_z%02d.1.ent', k-1);
    write_freeform_ent(H,X,Y,n,tbase,templatePath,outPath);
    fprintf('  zone %d: dx=%.3f -> %s (apex_x=%+.3f)\n', k-1, dx, outPath, apex);
end
fprintf('완료. build_graded_mla.m 로 zone/위치에 매핑.\n');


%% ---- 편심 raised-cosine 돔 (정점을 +x 로 dx 만큼 이동) ----
function H = offset_dome(X,Y,Rap,h,dx)
    rho = hypot(X-dx, Y)/Rap;                 % 중심을 (dx,0) 으로 이동
    H = h.*0.5.*(1+cos(pi*min(max(rho,0),1)));
    % rim window: 셀 조리개(r<=Rap)에서 테두리 0 강제 (편심분 클리핑 포함)
    r=hypot(X,Y); W=ones(size(r)); rw=0.85*Rap;
    z=(r-rw)/(Rap-rw); m=(r>=rw)&(r<Rap);
    W(m)=1-(3*z(m).^2-2*z(m).^3); W(r>=Rap)=0;
    H=max(H.*W,0);
end

%% ---- .ent writer (v2: Z=+H, rear=-tbase, 원점 템플릿) ----
function write_freeform_ent(H,X,Y,n,tbase,templatePath,outPath)
    Z=H; Xv=X(:); Yv=Y(:); Zv=Z(:); N=n*n;
    tpl=fileread(templatePath);
    tok=regexp(tpl,'ORAStartData;([\s\S]*?)ORAEndData;','tokenExtents');
    s0=tok{1}(1); e0=tok{1}(2);
    buf=sprintf('0 1 %d %d 0 0 %d 0 0 0',n,n,N);
    for i=1:N, buf=[buf sprintf(' %.17g %.17g %.17g',Xv(i),Yv(i),Zv(i))]; end %#ok<AGROW>
    buf=[buf ' 0 0 4 CartesianMapper 1 0 0 0 0'];
    newtxt=[tpl(1:s0-1) char(10) buf char(10) tpl(e0+1:end)];
    newtxt=regexprep(newtxt, ...
        '(CSGLensSurfacePrimitive_1[\s\S]*?setPosition:  \{ 0\. 0\. )[-0-9.eE]+(  \} ;)', ...
        ['$1' num2str(-tbase,'%g') '$2'],'once');
    newtxt=regexprep(newtxt,'restoreSmoothResample: "Yes"','restoreSmoothResample: "No"','once');
    fid=fopen(outPath,'w'); fwrite(fid,newtxt); fclose(fid);
end
