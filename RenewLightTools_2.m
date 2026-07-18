function RenewLightTools_2(path1, path2)
% RENEWLIGHTTOOLS_2  LightTools 2개 인스턴스(SweptEntity + 배열모델) 강제 재시작 후 연결.
% 전역 ID_LT, ID_swept, ltml, ltloc 를 설정한다. (BO_Freeform3D_asym.m /
% test_freeform_geom.m 공용)
%
%   [속도 개선 - 원본 대비 변경점]
%    (1) "start /min" 으로 실행 -> 창을 최소화 상태로 띄워 렌더링 비용 절감.
%        (LightTools 는 별도 배치/headless 실행파일이 없음을 확인함 - 완전
%         무음실행은 지원하지 않으므로 이게 현실적 상한.)
%    (2) 고정 pause(2)/pause(5) 를 "실제로 준비될 때까지 폴링 + 상한시간"으로
%        교체. 로딩이 예상보다 빨리 끝나면 그만큼 즉시 다음 단계로 진행하고,
%        상한(타임아웃)을 둬서 문제 상황에서 무한대기하지 않는다.
%    (3) PID 파싱 로직(tokens{1}{1}, tokens{3}{1})은 원본 그대로 유지 -
%        tasklist CSV 출력 형식에 의존하므로 검증된 인덱싱을 바꾸지 않음.
%
%   ※ 경로/사용자명은 사용자 머신에 맞게 아래에서 수정.
global ID_LT ID_swept ltml ltloc
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
LT1 = path1;  LT2 = path2;

fprintf('--- Restarting LightTools ---\n');
target_user = 'jhkim';
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);

% ===== 0) 기존 lt.exe 강제 종료 + 완전히 사라질 때까지 폴링 =====
system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
tKillTimeout = 10;  t0 = tic;
while toc(t0) < tKillTimeout
    [~, cmdout] = system(find_cmd);
    if ~contains(cmdout, 'lt.exe'), break; end
    pause(0.3);
end

% ===== 1) 첫 번째 인스턴스(SweptEntity) 실행 =====
system(sprintf('start /min "" "%s" "%s"', lt_exe_path, LT1));
try
    ltml  = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패. 라이선스/설치 확인.');
end

% lt.exe 프로세스가 최소 1개 뜰 때까지 폴링 (원본은 즉시 1회 조회만 했음;
% start /min 은 반환이 더 빨라 아직 안 떴을 수 있어 안전하게 폴링으로 보강)
tFind1Timeout = 20;  t1 = tic;  found1 = false;
while toc(t1) < tFind1Timeout
    [status, cmdout] = system(find_cmd);
    if status == 0 && contains(cmdout, 'lt.exe')
        found1 = true;
        break;
    end
    pause(0.3);
end
if ~found1
    error('lt.exe(swept) 탐색 실패');
end
tokens = regexp(cmdout, '"(\d+)"', 'tokens');
ID_swept = str2double(tokens{1}{1});
fprintf('PID(swept)=%d\n', ID_swept);

% ===== 2) 두 번째 인스턴스(배열모델) 실행 =====
system(sprintf('start /min "" "%s" "%s"', lt_exe_path, LT2));

% 두 번째 lt.exe 행까지 tasklist 에 나타날 때까지 폴링 (tokens{3} 를 쓰려면
% 최소 2개 행이 파싱돼야 하므로, token 개수로 "두 번째 프로세스 등장"을 판정)
tFind2Timeout = 20;  t2 = tic;  found2 = false;
while toc(t2) < tFind2Timeout
    [status, cmdout] = system(find_cmd);
    if status == 0 && contains(cmdout, 'lt.exe')
        tokens = regexp(cmdout, '"(\d+)"', 'tokens');
        if numel(tokens) >= 3
            found2 = true;
            break;
        end
    end
    pause(0.3);
end
if ~found2
    error('lt.exe(LT) 탐색 실패');
end
ID_LT = str2double(tokens{3}{1});
fprintf('PID(LT)=%d\n', ID_LT);

% ===== 3) COM 준비 확인 (고정 pause(5) 대체) =====
% 두 인스턴스 모두 실제로 LTCmd 를 받을 수 있는 상태가 될 때까지 폴링.
% 하나라도 준비 안 되면 그 인스턴스만 재시도하고, 상한시간 도달 시 경고만
% 남기고 진행(이후 호출에서 자연스럽게 실패하며 상위 재시작 루틴이 처리).
tReadyTimeout = 20;
for whichID = [ID_swept, ID_LT]
    tR = tic;  ready = false;
    while toc(tR) < tReadyTimeout
        try
            lt = ltloc.GetLTAPI(whichID);
            ltml.LTCmd(lt, 'Message "Check Connection"');
            ready = true;
            break;
        catch
            pause(0.5);
        end
    end
    if ~ready
        fprintf('[경고] PID=%d 가 %d초 내 COM 준비 확인 실패. 계속 진행하되 이후 호출이 실패할 수 있음.\n', ...
            whichID, tReadyTimeout);
    end
end
end
