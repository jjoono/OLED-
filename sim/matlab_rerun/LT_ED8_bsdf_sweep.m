% LightTools BSDF sweep for Extended Data Fig. 8 (run on the PC with LightTools and the .lts open).
% Same COM calls as the original hexagonal-MLA BSDF script; only the parameter lists changed.
%   lens / substrate index 1.80
%   (a) aspect-ratio sweep : Scale = 1 (close-packed), StretchZ = AR = 0.02, 0.10:0.05:1.50
%   (b) fill-factor sweep  : AR = 0.5 and 1.0, FF = 0.8:-0.1:0.1  ->  Scale = sqrt(FF/0.9069)
%       ASSUMES 'Scale' shrinks the lens at a fixed pitch (FF = 0.9069*Scale^2).
%       If 'Scale' enlarges the pitch instead, use Scale = sqrt(0.9069/FF).
% Output: bsdf/<id>.mat (variable BSDF, 180 x 90, each column normalized to 1),
%         same layout as hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat,
%         plus ED8_BSDF_LightTools_all.mat with everything.

lt   = actxserver('LightTools.ltapi');
ltml = actxserver('ltcom64.LTAPI2');

n_MLA = 1.80;
FF0   = pi/(2*sqrt(3));                                  % 0.9069, close-packed

AR_a  = [0.02 0.10:0.05:1.50];                            % 30 runs
FF_b  = 0.8:-0.1:0.1;  AR_b = [0.5 1.0];                   % 16 runs

runs = {};                                               % id, Scale, StretchZ
for ar = AR_a
    runs(end+1,:) = {sprintf('AR%s', strrep(sprintf('%.2f',ar),'.','p')), 1, ar};
end
for ar = AR_b
    for ff = FF_b
        runs(end+1,:) = {sprintf('FF%s_AR%s', strrep(sprintf('%.4f',ff),'.','p'), strrep(sprintf('%.1f',ar),'.','p')), sqrt(ff/FF0), ar};
    end
end
nrun = size(runs,1);
if ~exist('bsdf','dir'), mkdir('bsdf'); end

% index of the lens/substrate material
List = ltml.LTDbList(lt,'lens_manager[1]','CONSTANT_INDEX');
Key  = ltml.LTListAtPos(lt,List,1);
ltml.LTDbSet(lt,Key,'Index',n_MLA);

BSDF_all = zeros(180,90,nrun); raw_all = zeros(180,90,nrun);
for r = 1:nrun
    id = runs{r,1}; scale = runs{r,2}; stretchZ = runs{r,3};
    if exist(fullfile('bsdf',[id '.mat']),'file'), fprintf('skip %s (exists)\n', id); continue; end   % lets you resume
    List = ltml.LTDbList(lt,'TEXTURE_ZONE_EXTENT[1]','TEXTURE_PARAMETER');
    Key  = ltml.LTListByName(lt,List,'Scale');    ltml.LTDbSet(lt,Key,'Value',scale);
    List = ltml.LTDbList(lt,'TEXTURE_ZONE_EXTENT[1]','TEXTURE_PARAMETER');
    Key  = ltml.LTListByName(lt,List,'StretchZ'); ltml.LTDbSet(lt,Key,'Value',stretchZ);
    raw = zeros(180,90);
    for i = 1:90
        List = ltml.LTDbList(lt,'LENS_MANAGER[1]','AIM_SPHERE');
        Key  = ltml.LTListAtPos(lt,List,1);
        ltml.LTDbSet(lt,Key,'Upper_Angle',i-1);
        ltml.LTDbSet(lt,Key,'Lower_Angle',i);
        ltml.LTBegin(lt); ltml.LTCmd(lt,'\V3D BeginAllSimulations'); ltml.LTEnd(lt);
        for j = 1:180
            raw(j,i) = ltml.LTDbGet(lt,'receiver[1].mesh[1]','CellFlux_UI',1,j);
        end
    end
    BSDF = raw ./ repmat(sum(raw,1),180,1);              % each incidence angle normalized to 1
    w = cosd(0:89).*sind(0:89); p = sum(sum(BSDF(1:90,:),1).*w)/sum(w);
    fprintf('%-16s Scale %.4f StretchZ %.2f | p(cos sin) = %.4f\n', id, scale, stretchZ, p);
    save(fullfile('bsdf',[id '.mat']),'BSDF','raw','n_MLA','scale','stretchZ');
    BSDF_all(:,:,r) = BSDF; raw_all(:,:,r) = raw;
end
% close-packed rows of the fill-factor sweep are the same geometries as AR 0.50 and AR 1.00
copyfile(fullfile('bsdf','AR0p50.mat'), fullfile('bsdf','FF0p9069_AR0p5.mat'));
copyfile(fullfile('bsdf','AR1p00.mat'), fullfile('bsdf','FF0p9069_AR1p0.mat'));
ids = runs(:,1); scales = cell2mat(runs(:,2)); stretchZs = cell2mat(runs(:,3));
save('ED8_BSDF_LightTools_all.mat','BSDF_all','raw_all','ids','scales','stretchZs','n_MLA');

% check against the BSDF used in the paper (n = 1.80, close-packed hemisphere = AR1p00)
S = load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat'); B0 = S.BSDF_MLA(:,:,11);
N = load(fullfile('bsdf','AR1p00.mat'));
fprintf('AR1p00 vs original slice 11: max |dB_T| = %.4f, max |dB| = %.4f\n', ...
    max(abs(sum(N.BSDF(1:90,:),1)-sum(B0(1:90,:),1))), max(abs(N.BSDF(:)-B0(:))));
