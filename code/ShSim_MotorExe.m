function [uph, iph, loss, slip, pf] = ShSim_MotorExe(nRPM, Tq, motor)
% ShSim_MotorExe
% 0.1.1.1
%
%   Description
%       With a motor map, this function returns the data from the operating
%       defined point.
%
%   Input signals
%       n - motor speed in rpm
%       t - motor torue in Nm
%       motor - motor map struct containing:
%           rpmBaseMotor
%           rpmBaseGenerator
%           TorqueMT, TorqueGT - Torque point in matrix below base speed (rpmBase...) 
%           PowerMP, PowerGP - Power input (n*t) above base speed in matrix
%           SpeedMT, SpeedGT, SpeedMP, SpeedGP - speed points in matrix
%           UPhaseMT,..GT,..MP,..GP - U-phase return value matrix
%           IPhaseMT,..GT,..MP,..GP - I-phase return value matrix
%           TotLossMT,..GT,..MP,..GP - total losses in kW value matrix
%           PFactMT,..GT,..MP,..GP - power factor return value matrix
%   Output signals
%       u - Motor voltage in V
%       iph - Motor current in Arms
%       loss - power loss in W
%       slip - motor slip
%       pf - power factor
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-09  JR      New coding design (spaces) plus motor loss in Watts
%
if any(Tq > motor.TorqueMT(1, end))
    b = Tq > motor.TorqueMT(1, end);
    fprintf(['Warning: Max torque exceeded: T_min=' num2str(max(Tq)) ' Nm. Limit=' num2str(motor.TorqueMT(1,end),'%4.1f') ' Nm.\n']);
    Tq(b) = motor.TorqueMT(1, end);
elseif any(Tq < motor.TorqueGT(1, end))
    b = Tq < motor.TorqueGT(1, end);
    fprintf(['Warning: Max braking torque exceeded: T_min=' num2str(min(Tq)) ' Nm. Limit=' num2str(motor.TorqueGT(1,end),'%4.1f') ' Nm.\n']);
    Tq(b) = motor.TorqueGT(1, end);
end
if any(Tq .* nRPM / 60 * 2 * pi > motor.PowerMP(1, end))
    b = find(Tq .* nRPM / 60 * 2 * pi > motor.PowerMP(1, end));
    fprintf(['Warning: Max motor power exceeded: P_Rail=' num2str(max(Tq.*nRPM/60*2*pi)/1000) 'kW. Limit=' num2str(motor.PowerMP(1,end)/1000,'%4.1f') ' kW.\n']);
    Tq(b) = motor.PowerMP(1, end) * 60 ./ (nRPM(b) * 2 * pi);
elseif any(Tq .* nRPM / 60 * 2 * pi < motor.PowerGP(1, end))
    b = find(Tq .* nRPM / 60 * 2 * pi < motor.PowerGP(1,end));
    fprintf(['Warning: Max generator power exceeded: P_Rail=' num2str(min(Tq.*nRPM/60*2*pi)/1000) 'kW. Limit=' num2str(motor.PowerGP(1,end)/1000,'%4.1f') ' kW.\n']);
    Tq(b) = motor.PowerGP(1, end) * 60 ./ (nRPM(b) * 2 * pi);
end
if any(nRPM > motor.SpeedMP(end, end))
    b = nRPM > motor.SpeedMP(end, end);
    fprintf(['Warning: Max motor speed exceeded: n_rpm=' num2str(max(nRPM),'%3.1f') ' rpm. Limit=' num2str(motor.SpeedMP(end,end),'%3.1f') ' rpm.\n']);
    nRPM(b) = motor.SpeedMP(end, end);
end

b = find(nRPM <= motor.rpmBaseMotor & Tq >= 0);
if ~isempty(b)
    uph(b,1)  = int2(motor.TorqueMT, motor.SpeedMT, motor.UPhaseMT,  Tq(b), nRPM(b));
    iph(b,1)  = int2(motor.TorqueMT, motor.SpeedMT, motor.IPhaseMT,  Tq(b), nRPM(b));
    loss(b,1) = int2(motor.TorqueMT, motor.SpeedMT, motor.TotLossMT, Tq(b), nRPM(b));
    pf(b,1)   = int2(motor.TorqueMT, motor.SpeedMT, motor.PFactMT,   Tq(b), nRPM(b));
    slip(b,1) = int2(motor.TorqueMT, motor.SpeedMT, motor.SlipMT,    Tq(b), nRPM(b));
end
b=find(nRPM>motor.rpmBaseMotor & Tq>=0);
if ~isempty(b)
    uph(b,1)  = int2(motor.PowerMP, motor.SpeedMP, motor.UPhaseMP,  Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    iph(b,1)  = int2(motor.PowerMP, motor.SpeedMP, motor.IPhaseMP,  Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    loss(b,1) = int2(motor.PowerMP, motor.SpeedMP, motor.TotLossMP, Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    pf(b,1)   = int2(motor.PowerMP, motor.SpeedMP, motor.PFactMP,   Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    slip(b,1) = int2(motor.PowerMP, motor.SpeedMP, motor.SlipMP,    Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
end
b=find(nRPM<=motor.rpmBaseGenerator & Tq<0);
if ~isempty(b)
    uph(b,1)  = int2(motor.TorqueGT, motor.SpeedGT, motor.UPhaseGT,  Tq(b), nRPM(b));
    iph(b,1)  = int2(motor.TorqueGT, motor.SpeedGT, motor.IPhaseGT,  Tq(b), nRPM(b));
    loss(b,1) = int2(motor.TorqueGT, motor.SpeedGT, motor.TotLossGT, Tq(b), nRPM(b));
    pf(b,1)   = int2(motor.TorqueGT, motor.SpeedGT, motor.PFactGT,   Tq(b), nRPM(b));
    slip(b,1) = int2(motor.TorqueGT, motor.SpeedGT, motor.SlipGT,    Tq(b), nRPM(b));
end
b=find(nRPM>motor.rpmBaseGenerator & Tq<0);
if ~isempty(b)
    uph(b,1)  = int2(motor.PowerGP, motor.SpeedGP, motor.UPhaseGP,  Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    iph(b,1)  = int2(motor.PowerGP, motor.SpeedGP, motor.IPhaseGP,  Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    loss(b,1) = int2(motor.PowerGP, motor.SpeedGP, motor.TotLossGP, Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    pf(b,1)   = int2(motor.PowerGP, motor.SpeedGP, motor.PFactGP,   Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
    slip(b,1) = int2(motor.PowerGP, motor.SpeedGP, motor.SlipGP,    Tq(b) .* nRPM(b) / 60 * 2 * pi, nRPM(b));
end
return % ShSim_MotorExe