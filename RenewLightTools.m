function RenewLightTools()
% RENEWLIGHTTOOLS  배열(MLA) 시뮬 LightTools 인스턴스 1개를 강제 재시작 후 연결.
% 전역 ID_LT, ltml, ltloc 를 설정한다.
%
%   [변경] freeform .ent 를 MATLAB 이 직접 생성해 텍스처에 물리는 방식으로 바뀌면서
%   SweptEntity 인스턴스가 더 이상 필요 없다(과거엔 geometry 를 SaveLibrary 하려고
%   2개를 띄웠음). 이제 배열 모델(Lens_size_effect...lts) 1개만 띄운다.
%   ※ 경로/사용자명은 사용자 머신에 맞게 아래에서 수정.
global ID_LT ltml ltloc
lt_exe_path        = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
model_file_path_LT = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\Lens_size_effect_for_PSO_bump_modified_v1.1.lts';

fprintf('--- Restarting LightTools (array model) ---\n');
target_user = 'jhkim';
system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
pause(2);

system(sprintf('"%s" "%s" &', lt_exe_path, model_file_path_LT));
try
    ltml  = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패. 라이선스/설치 확인.');
end

find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);
[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    ID_LT = str2double(tokens{1}{1});
    fprintf('PID(LT)=%d\n', ID_LT);
else
    error('lt.exe(LT) 탐색 실패');
end
pause(5);
end
