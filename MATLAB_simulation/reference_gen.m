function r = reference_gen(t, scen)
%REFERENCE_GEN scen 5 = lemniscate, all other scenarios = set-point.
%#codegen
if scen == 5, r = reference_gen_base(t, 5); else, r = reference_gen_base(t, 1); end
end


