function vehicle = ShSim_CalcMassRes(vehicle)
% ShSim_CalcMassRes
% 0.1.1.1
%
% vehicle.haul_mass
% vehicle.haul_length
% vehicle.haul_push
% vehicle.haul_a_dash
% vehicle.haul_a_actual
% vehicle.haul_b1_dash
% vehicle.haul_b2
% vehicle.haul_c_front
% vehicle.haul_c_len
% vehicle.haul_c_front2
% vehicle.haul_curve_cr0
% vehicle.haul_curve_cr1
% vehicle.haul_curve_cr2
% vehicle.load_cond
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-26  JR      Code clean up
%

if isempty(vehicle.ws_mass{vehicle.ws_selection})
    ws_mass0 = 0;
else
    ws_mass0 = vehicle.ws_mass{vehicle.ws_selection} * (vehicle.nAxDrive + vehicle.nAxTrail) * 2;
end
if isempty(vehicle.ws_rotMass{vehicle.ws_selection})
    ws_rotMass0 = 0;
else
    ws_rotMass0 = vehicle.ws_rotMass{vehicle.ws_selection} * (vehicle.nAxDrive + vehicle.nAxTrail) * 2;
end
vehicle.statmassLoad = vehicle.TareMassUnit + [vehicle.LoadMassUnit{:}] + ws_mass0;
vehicle.totMass = vehicle.statmassLoad + vehicle.RotationMassUnit + ws_rotMass0;

vehicle.davies_a = max(vehicle.statmassLoad * vehicle.davies_a_dash, vehicle.davies_a_actual) / 1000;
vehicle.davies_b = (vehicle.davies_b1_dash * vehicle.statmassLoad + vehicle.davies_b2) / 1000;
vehicle.davies_c1 = (vehicle.davies_c_front + vehicle.davies_c_len * vehicle.unitLength * vehicle.nUnits) / (vehicle.nUnits * 1000);
% vehicle.davies_c1 = (vehicle.davies_c_front + vehicle.davies_c_len * vehicle.unitLength * (1:10)) ./ ((1:10) * 1000);
vehicle.davies_c2 = (vehicle.c_front2 + vehicle.davies_c_len * vehicle.unitLength * vehicle.nUnits) / (vehicle.nUnits * 1000);
% vehicle.davies_c2 = (vehicle.c_front2 + vehicle.davies_c_len * vehicle.unitLength * (1:10)) ./ ((1:10) * 1000);
return % ShSim_CalcMassRes

