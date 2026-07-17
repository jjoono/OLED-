%% generate_asym_max_ents.m
% "제작 가능한 클래스 안에서 비대칭 극대화" 렌즈 5종 .ent 생성.
% v2 규약(원점 템플릿, +Z 돌출, rear=-tbase, SmoothResample No) 그대로.
%
% [형상] raised-cosine 돔의 x축을 앞/뒤 비대칭 스케일:
%   앞쪽(+x)은 s_steep 로 압축 -> 가파른 전면 플랭크(정점 앞쏠림),
%   뒤쪽(-x)은 s_gentle 로 확장 -> 완만한 긴 후면 램프. 정점(apex)이 +x로 이동.
% [제작성 강제]
%   - rim_window: r>=0.85*Rap 에서 smootherstep 로 0 -> 테두리 높이 0 (타일링 가능)
%   - 최대 draft(경사) <= ~65도 (임프린트 이형 한계) 가 되도록 s_steep 튜닝
%   - 닫힌 돔(rim=0) 이므로 평균 기울기=0: 빔스티어링은 균일 프리즘이 아니라
%     "가파른/완만한 플랭크 대비 + TIR + 비선형 굴절" 에서 나옴.
%
% [DOF] 이 스크립트는 고정 5종을 뽑는 데모. BO 변수화하려면 (h,xc,s_steep,s_gentle)
%   4개를 파라미터로 노출하면 됨(아래 asym_dome 그대로 목적함수에서 호출).

templatePath = 'freeform_template_v2.ent';
Ra=1.2139; Rap=1.0; n=141; tbase=0.30; DRAFT_MAX=65;
g=linspace(-Ra,Ra,n); [X,Y]=meshgrid(g,g); r=hypot(X,Y);

designs = {
 'ff_x1_steepfront_mild',  0.72, 0.18, 1.35, 0.80
 'ff_x2_steepfront_med',   0.72, 0.26, 1.55, 0.72
 'ff_x3_teardrop',         0.72, 0.33, 1.75, 0.66
 'ff_x4_teardrop_strong',  0.72, 0.40, 1.80, 0.60
 'ff_x5_max',              0.72, 0.46, 1.92, 0.55 };

fprintf('%-24s %6s %7s %7s  moldable\n','name','draft','rim','apex_x');
for i=1:size(designs,1)
    name=designs{i,1}; h=designs{i,2}; xc=designs{i,3};
    s_steep=designs{i,4}; s_gentle=designs{i,5};
    H = asym_dome(X,Y,r,Rap,h,xc,s_steep,s_gentle);
    [gx,gy]=gradient(H,g,g); slope=atand(hypot(gx,gy));
    d=max(slope(r<=Rap)); rim=max(abs(H(r>0.97 & r<=Rap)));
    [~,im]=max(H(:)); apex_x=X(im);
    ok = (d<=DRAFT_MAX+0.5) && (rim<0.01);
    write_freeform_ent(H,X,Y,n,tbase,templatePath,[name '.1.ent']);
    fprintf('%-24s %5.1fd %7.4f %+7.3f   %s\n',name,d,rim,apex_x,string(ok));
end
fprintf('완료. LightTools 로 로드해 확인.\n');


%% ---- 형상 (BO 목적함수에서 그대로 재사용) ----
function H = asym_dome(X,Y,r,Rap,h,xc,s_steep,s_gentle)
    xr=X-xc;
    sx=s_gentle*ones(size(X)); sx(xr>0)=s_steep;
    Xe=xr.*sx;
    rho=sqrt(Xe.^2+Y.^2)/Rap;
    H=h.*0.5.*(1+cos(pi*min(max(rho,0),1)));
    % rim window: r>=0.85Rap 에서 1->0 (테두리 높이 0 강제 -> 타일링)
    W=ones(size(r)); rw=0.85*Rap;
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
