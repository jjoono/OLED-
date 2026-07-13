function RenewLightTools()
% RENEWLIGHTTOOLS  LightTools 2개 인스턴스(SweptEntity + 배열모델) 강제 재시작 후 연결.
% 전역 ID_LT, ID_swept, ltml, ltloc 를 설정한다. (BO_Freeform3D_asym.m /
% test_freeform_geom.m 공용)
%   ※ 경로/사용자명은 사용자 머신에 맞게 아래에서 수정.
global ID_LT ID_swept ltml ltloc
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
model_file_path_swept = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\SweptEntity.2.lts';
model_file_path_LT    = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\Lens_size_effect_for_PSO_bump_modified_v1.1.lts';

fprintf('--- Restarting LightTools ---\n');
target_user = 'jhkim';
system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
pause(2);

system(sprintf('"%s" "%s" &', lt_exe_path, model_file_path_swept));
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
    ID_swept = str2double(tokens{1}{1});
    fprintf('PID(swept)=%d\n', ID_swept);
else
    error('lt.exe(swept) 탐색 실패');
end

system(sprintf('"%s" "%s" &', lt_exe_path, model_file_path_LT));
[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    ID_LT = str2double(tokens{3}{1});
    fprintf('PID(LT)=%d\n', ID_LT);
else
    error('lt.exe(LT) 탐색 실패');
end
pause(5);
end
