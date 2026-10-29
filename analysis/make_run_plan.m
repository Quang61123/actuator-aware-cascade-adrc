function T = make_run_plan()
%MAKE_RUN_PLAN Fixed list of all redesign SITL runs, written BEFORE any
%   run (run_plan.csv). Every listed run is executed once and reported,
%   failed runs included; no run is repeated selectively.
%   Controller ids of the PX4 gain table (create_sitl_controller):
%     1 PID-Budget, 2 AA-Budget, 4 AA-HalfGain, 9 AA-Budget-S, 10 AA-HalfGain-S
%     (unchanged), 15 PID-Budget-S (new), 16 A1-Commanded (new),
%     17-19 PID-E, AA-Budget-E, AA-HalfGain-E (iteration 0),
%     20-22 reserved for a second -E iteration, which was not needed (iteration_decision.txt);
%     the compiled gain table assigns ids 20-22 to A1-Commanded and RO-HalfGain (ratio >= 1) and PID-Budget (b*).
%   Phase P1:
%     E0H  hover 60 s, ids 17-19, seeds 96001-96010
%     PSS  S6 40 s, id 15, seeds 91001-91030
%     PSH  hover 60 s, id 15, seeds 96001-96010
%     A1S  S6 40 s, id 16, seeds 91001-91010
%     LH   hover 300 s, ids 1 2 4 15 9 10, fit seeds 97001-97005 and
%          held-out seeds 97006-97010
%   Phase P1b_unused (NOT to be run: ids 20-22 were later reassigned, see above):
%     E1H_unused  hover 60 s, ids 20-22, seeds 96001-96010
%   Phase P2 (final -E gains, "EF" = 17-19):
%     EFS  S6, seeds 91001-91030; EF5 S5, seeds 91001-91010;
%     EFH  hover 60 s, seeds 96001-96010 (= E0H, not rerun);
%     EFL  hover 300 s, seeds 97001-97010 (5 fit + 5 held-out)
here = fileparts(mfilename('fullpath')); RA = fullfile(here,'..','results','analysis');
rows = cell(0, 9);
    function a(ph, grp, ids, scen, seeds, Tend, win, split)
        for s = seeds, for id = ids
            sp = ""; if ~isempty(split), if s <= split, sp = "fit"; else, sp = "heldout"; end, end
            rows(end+1,:) = {ph, grp, string(id), scen, s, Tend, win(1), win(2), sp}; %#ok<AGROW>
        end, end
    end
a("P1","E0H",17:19,"hover",96001:96010,60,[20 60],[]);
a("P1","PSS",15,"S6",91001:91030,40,[10 40],[]);
a("P1","PSH",15,"hover",96001:96010,60,[20 60],[]);
a("P1","A1S",16,"S6",91001:91010,40,[10 40],[]);
a("P1","LH",[1 2 4 15 9 10],"hover",97001:97010,300,[20 300],97005);
a("P1b_unused","E1H_unused",20:22,"hover",96001:96010,60,[20 60],[]);
for k = 1:3
    a("P2","EFS","EF"+k,"S6",91001:91030,40,[10 40],[]);
    a("P2","EF5","EF"+k,"S5",91001:91010,45,[5 45],[]);
    a("P2","EFL","EF"+k,"hover",97001:97010,300,[20 300],97005);
end
T = cell2table(rows, 'VariableNames', {'Phase','Group','CtrlId','Scenario','Seed','Tend_s','WinStart_s','WinEnd_s','Split'});
T.RunID = (1:height(T))';
T = T(:, [end 1:end-1]);
writetable(T, fullfile(RA,'run_plan.csv'));
end
