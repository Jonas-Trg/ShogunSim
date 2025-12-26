function F_target = ShSim_maxBE(tvehicle, speed, F_slope, F_res)
% ShSim_maxBE
% 0.1.2.3
%
% Ver       Date        Sign    Descr
% 0.1.2.2   2025-12-19  JR      Cleaning code.
% 0.1.2.3   2025-12-19  JR      Davies_a, _b and cn made single value
%
F_target = -tvehicle.retard_max * tvehicle.totMass;
if tvehicle.utiliseEBrake && speed > 0
    F_target = max([-tvehicle.retard_min * tvehicle.retard_min_speed / speed * tvehicle.totMass, F_target]);
end
switch tvehicle.RetardMethod
    case 1 % constant retardation with davies a
        F_target = F_target + tvehicle.davies_a;
    case 2 % constant retardation with davies all
        force_davies = (tvehicle.davies_a + tvehicle.davies_b * speed * 3.6 + tvehicle.davies_c1 * (speed * 3.6)^2);
        F_target = F_target + force_davies;
    case 3 % contant retardation
        F_target = -tvehicle.retard_max * tvehicle.totMass + F_res + F_slope;
end
% Include Total brake curve
x = find(tvehicle.TotBrake_v > speed, 1, 'first') - 1;
if ~isempty(x)
    d = (speed - tvehicle.TotBrake_v(x)) / (tvehicle.TotBrake_v(x + 1) - tvehicle.TotBrake_v(x));
    if x > 1
        P = d * tvehicle.TotBrake_p(x + 1) + (1 - d) * tvehicle.TotBrake_p(x);
        Brk_F = P / speed;
    else
        Brk_F = d * tvehicle.TotBrake_f(x + 1) + (1 - d) * tvehicle.TotBrake_f(x);
    end
    F_target = max([-Brk_F, F_target]);
end  % ShSim_maxBE
