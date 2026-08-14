% ============================================================
%  night_run.m — 자기 전 한 번 실행. 세 단계를 순서대로 돌린다.
%
%   1) warmstart_from_hemisphere  : arm 4 (60-80 deg) — Fig 2f / Table S6 완성.
%                                   merge 로 arm 1,2,3,5 는 보존된다.  (~2 h)
%   2) reeval_confirm_2040        : 경계선(t=2.1)이던 20-40 deg 를 반복 5회로
%                                   확정. 탐색 없음.                   (~2.5 h)
%   3) check_patch_convergence    : 최고 설계를 patch 15/25/35 로 재평가,
%                                   "25 면 충분한가" 에 답.            (~2 h)
%
%  합계 약 6.5 h. 각 단계는 자체 저장/로그를 가지므로 중간에 죽어도
%  앞 단계 결과는 남고, 다시 실행하면 1) 은 merge 로 이어붙는다.
%  단계를 빼려면 해당 줄을 주석 처리.
% ============================================================
try, warmstart_from_hemisphere;  catch ME, fprintf(2,'%s\n',getReport(ME,'extended')); end
try, reeval_confirm_2040;        catch ME, fprintf(2,'%s\n',getReport(ME,'extended')); end
try, check_patch_convergence;    catch ME, fprintf(2,'%s\n',getReport(ME,'extended')); end
fprintf('\n[night_run] 전 단계 종료. *_result.mat 3개와 로그를 확인할 것.\n');
