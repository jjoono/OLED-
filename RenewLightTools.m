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
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);

system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
% [속도] 고정 pause(2) 대신, 프로세스가 실제로 사라질 때까지만 폴링 대기
% (헤드리스 실행파일이 없어 GUI 실행 자체는 못 줄이지만, 이 고정대기는 줄일 수 있음)
tKillTimeout = 10;  tKill0 = tic;
while toc(tKill0) < tKillTimeout
    [~, cmdout] = system(find_cmd);
    if ~contains(cmdout, 'lt.exe'), break; end
    pause(0.3);
end

% [속도] "start /min" 으로 최소화 실행 -> 창 렌더링 비용 절감(완전 headless는 아님)
system(sprintf('start /min "" "%s" "%s"', lt_exe_path, model_file_path_LT));
try
    ltml  = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패. 라이선스/설치 확인.');
end

[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    ID_LT = str2double(tokens{1}{1});
    fprintf('PID(LT)=%d\n', ID_LT);
else
    error('lt.exe(LT) 탐색 실패');
end

% [속도] 고정 pause(5) 대신, COM 이 실제로 명령을 받을 준비가 될 때까지만 폴링
% (모델 로딩이 5초보다 빨리 끝나면 그만큼 즉시 진행; 상한 20초로 무한대기 방지)
tReadyTimeout = 20;  tReady0 = tic;  ready = false;
while toc(tReady0) < tReadyTimeout
    try
        lt = ltloc.GetLTAPI(ID_LT);
        ltml.LTCmd(lt, 'Message "Check Connection"');
        ready = true;
        break;
    catch
        pause(0.5);
    end
end
if ~ready
    fprintf('[경고] %d초 내 COM 준비 확인 실패. 계속 진행하되 이후 호출이 실패할 수 있음.\n', tReadyTimeout);
end
end
