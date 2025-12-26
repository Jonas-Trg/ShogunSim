function [TE_achieved, ULine, ULink, ILine, BRPower, iph, IS, mSlip, swFreq, ESPower, PropLoss]=...
    ShSim_DC_TE(TE_demand, vehicle, speed, routedata, sim_config, limitation)
% ShSim_DC_TE
% 0.1.1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DC Traction mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Description
%       This function is called several times for every sample in traction
%       mode with in DC mode.
%
%   Input data
%       TE_demand - effort reference in kN
%       vehicle - struct containing at least:
%       speed - in m/s
%       routedata*, if used it shall contains:
%           USupply - Line voltage definition from line
%           * Can be left empty, i.e. routedata=[];
%       simulation
%           useSubStations - defines whether to use voltage from line supply or
%                            vehicle uNorm
%       limitation, struct containing at least:
%           thermEffort
%           thermPower
%
%   Output data:
%       TE - achieved effort on wheel
%       ULine - Line voltage
%       ULink - DC-link voltage
%       ILine - Line current
%       BRPower - Brake resistor power
%       iMotor - Motor current (future provision)
%       mSlip - slip in the machine (future provision)
%       fMode - operation mode (not used)
%       ESPower - power used in the Energy Saver
%       PropLoss - The actual loss in the propulsion system - total on train level
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-09  JR      New coding design (spaces) plus motor loss in Watts
% 0.1.1.2   2025-12-22  JR      Removed 'hasEfficiency'
%
ESPower = 0;
BRPower = 0;
iph = 0;
IS = 0;
mSlip = 0;

if vehicle.ES_Op_Mode > 0
    limitation.ES_ChargeLevel = 1;
else
    limitation.ES_ChargeLevel = 0;
end

if isempty(routedata)
    ULine = vehicle.MinUL;
else
    if sim_config.useSubStations
        ULine = routedata.USupply;
    else
        ULine = vehicle.MinUL;
    end
end

ZTot = vehicle.rLineInd / vehicle.nHVSystems;

% Max tractive effort
TE_achieved = ShSim_maxTE(TE_demand, vehicle, speed, limitation);
if ZTot > 0
    p = -ULine / ZTot;
    q = TE_achieved * speed * 1000 / ZTot / 0.9; % Power factor 0.9 must be explained
    ILine = -p / 2 - sqrt(p^2 / 4 - q);
    ULink0 = ULine - ILine * ZTot;
else
    ULink0 = ULine;
    ULink = ULine;
end

olddiff = inf;
doCalc = 2; % 2 = first loop, 1 = tuning iline from above, -1 = tuning iline from below
iLineLim = vehicle.ILDriveMax;
ILine = inf;
while doCalc ~= 0
    if TE_achieved == 0 && speed == 0
        pDcLink = 0;
        swFreq = 0;
    else
        % Drives and Motor
        nMotorRPM = vehicle.gearRatio * speed / (vehicle.wheelsize{vehicle.ws_selection} * pi) * 60; %revs per min

        % torque per motor in Nm
        tq = TE_achieved * vehicle.wheelsize{vehicle.ws_selection} * 1000 / (vehicle.gearRatio * vehicle.nAxDrive * 2);

        % Gearbox losses in Nm
        tq = tq + vehicle.gearLoss0 + abs(tq * vehicle.gearLoss1) + nMotorRPM * vehicle.gearLoss2;

        % Traction motor
        switch vehicle.motorMethod % 1 = fixed efficiency, 2 = motor map
            case 1 % fixed efficiency
                PInverter = tq * nMotorRPM * 2 * pi / 60; % W
                PInverter = PInverter + abs(PInverter * vehicle.TM_FixLoss);
                uph = min([(speed + 1) / vehicle.vBasePwr, 1]) * ULink * sqrt(6) / pi;
                iph = TE_achieved * speed * 1000/(sqrt(3) * uph * 0.9) / vehicle.nAxDrive;
                freq = nMotorRPM * vehicle.nPoles / 120;

            case 2 % motor map
                [uph, iph, mloss, mSlip, pf] = ShSim_MotorExe(nMotorRPM, tq, vehicle.motorMap); %#ok<ASGLU>
                PInverter = (nMotorRPM .* tq / 60 * 2 * pi + mloss); %W
                freq = nMotorRPM * vehicle.nPoles / 120;
        end

        % Inverter cable
        if vehicle.R_MotCable > 0
            PInverter = (PInverter + 3 * vehicle.R_MotCable * iph^2) * vehicle.nMotors;
        else
            PInverter = PInverter * vehicle.nMotors;
        end
        iInv = iph * vehicle.nMotors;

        % Inverter
        if isempty(vehicle.INV_FixLoss)
            b = find(freq > vehicle.FreqInvMotor, 1, 'last');
            if isempty(b)
                swFreq = vehicle.asyncSwFqInvMotor;
            else
                swFreq = freq * vehicle.nModesInvMotor(b);
            end
            InvLoss = vehicle.INV_RLoss * iInv^2 + vehicle.INV_FLoss * swFreq + vehicle.INV_FRLoss * swFreq * iInv^2;
            pDcLink = (PInverter + InvLoss) * vehicle.nINV; % Total DC-link power in the unit
        else
            swFreq = 0;
            pDcLink=(PInverter + abs(PInverter * vehicle.INV_FixLoss)) * vehicle.nINV;
        end
    end

    % DC-link - pDcLink is the DC-link power per HV system in Watt
    pLink0 = (pDcLink * vehicle.nHVSystems + vehicle.AuxPower / vehicle.APS_efficiency * 1000);

    % Line filter
    if ZTot > 0
        q = pLink0 / ZTot;
        ILine = -p / 2 - sqrt(p^2 / 4 - q);
        ULink = ULine - ILine * ZTot;
    else
        ULink = ULine;
        ILine = pLink0 / ULink;
    end

    if abs(1 - ULink / ULink0) > 0.005
        ULink0 = ULink;
    else
        if doCalc == 2 % first loop
            if ILine > iLineLim
                doCalc = 1;
            else
                doCalc = 0; % Done!
            end
        end
        diff = 1 - abs(iLineLim / ILine);
        if diff < 0.001
            doCalc = 0;
        else
            if diff < olddiff
                TE_achieved = ShSim_maxTE(TE_achieved * iLineLim / ILine, vehicle, speed, limitation);
            else
                TE_achieved = ShSim_maxTE(TE_achieved * (1 - diff * 0.5), vehicle, speed, limitation);
            end
            olddiff = diff;
        end
    end
end
PropLoss = ULine * ILine / 1000 - (TE_achieved * speed + vehicle.AuxPower / vehicle.APS_efficiency);
return % ShSim_DC_TE
