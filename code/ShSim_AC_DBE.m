function [TE, ULine, ULink, ILine, BRPower, iph, IS, mSlip, swFreq, ESPower, PropLoss] =...
    ShSim_AC_DBE(DBE_demand, vehicle, speed, routedata, sim_config, limitation)
% ShSim_AC_DBE
% 0.1.1.3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   AC Braking mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Description
%       This function is called several times for every sample in braking
%       mode during AC operation.

%   Input data
%       TE_demand - effort reference in kN
%       vehicle - struct containing at least:
%       speed - in m/s
%       routedata*, if used it shall contains:
%           URegen - defines the line voltage in braking
%           * Can be left empty, i.e. routedata=[];
%       simulation
%           useSubStations - defines whether to use voltage from line supply or
%                           vehicle.U_Line_Nom
%       limitation, struct containing at least:
%           thermEffort
%           thermPower
%   
%   Output data:
%       TE - achieved effort on wheel - total on train level
%       ULine - Line voltage
%       ULink - DC-link voltage
%       ILine - Line current - total on train level
%       BRPower - Brake resistor power - total on train level
%       iMotor - Motor current (future provision) - per motor
%       mSlip - slip in the machine (future provision)
%       fMode - operation mode (not used)
%       ESPower - power used in the Energy Saver (not used in AC mode at
%                 the moment)
%       PropLoss - The actual loss in the propulsion system - total on train level
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-09  JR      New coding design (spaces) plus motor loss in Watts
% 0.1.1.2   2025-12-18  JR      code error correction
% 0.1.1.3   2025-12-22  JR      removing 'hasEfficiency' + code cleaning
%
%
% ESPower=0;
BRPower = 0;
iph = 0;
mSlip = 0;

if isempty(routedata)
    if isinf(vehicle.OVClevel)
        ULine = vehicle.U_Line_Nom;
    else
        ULine = vehicle.OVClevel;
    end
    routedata.Regen_current = inf;
    routedata.Regen_percent = 1;
else
    if sim_config.useSubStations
        ULine = routedata.URegen;
    elseif isinf(vehicle.OVClevel)
        ULine = vehicle.U_Line_Nom;
    else
        ULine = vehicle.OVClevel;
    end
end

if isempty(vehicle.Trafo_FixLoss)
    baseTrafoLoss = (vehicle.trafoIdleLoss + vehicle.trafoHarmonicLoss) * vehicle.nHVSystems;
end


% Dynamic brake demand
TE = ShSim_maxDBE(DBE_demand, vehicle, speed, limitation);
% Actions
% incorporate:
% routedata.Regen_percent
% routedata.Regen_current
% pull-out torque
%
ILine = min(TE * speed * 1000 / ULine, -eps);
doCalc = 2; % 2 = first loop, 1 = tuning iline from above, -1 = tuning iline from below
iLineLim = max([-vehicle.ILBrakeMax, -routedata.Regen_current ILine * routedata.Regen_percent]);
ILine = -inf;
while doCalc ~= 0
    % Drives and Motor
    nMotorRPM = vehicle.gearRatio * speed / (vehicle.wheelsize{vehicle.ws_selection} * pi) * 60; % revs per min

    % torque per motor
    tq = TE * vehicle.wheelsize{vehicle.ws_selection} * 1000 / (vehicle.gearRatio * vehicle.nAxDrive * 2);

    % Gearbox losses
    tq = tq + vehicle.gearLoss0 + abs(tq * vehicle.gearLoss1) + nMotorRPM * vehicle.gearLoss2;

    % traction motor
    switch vehicle.motorMethod % 1 = fixed efficiency, 2 = motor map
        case 1 % fixed efficiency
            PInverter = tq * nMotorRPM * pi * 2 / 60;
            PInverter = PInverter + abs(PInverter * vehicle.TM_FixLoss);
            uph = min([speed / vehicle.vBasePwr, 1]) * vehicle.uDcLink * sqrt(6) / pi;
            iph = -TE * speed * 1000 / (sqrt(3) * uph * 0.9) / vehicle.nAxDrive;
            freq = nMotorRPM * vehicle.nPoles / 120;
        case 2 % motor map
            [uph, iph, mloss, mSlip, pf] = ShSim_MotorExe(nMotorRPM, tq, vehicle.motorMap); %#ok<ASGLU>
            PInverter = nMotorRPM .* tq / 60 * 2 * pi + mloss;
            freq = nMotorRPM * vehicle.nPoles / 120;
    end

    % Inverter cable
    if vehicle.R_MotCable > 0
        PInverter=(PInverter + 3 * vehicle.R_MotCable * iph^2) * vehicle.nMotors;
    else
        PInverter = PInverter * vehicle.nMotors;
    end
    iInv = iph * vehicle.nMotors;
    % Inverter
    if isempty(vehicle.INV_FixLoss)
        p = find(freq > vehicle.FreqInvGenerator, 1, 'last');
        if isempty(p)
            swFreq = vehicle.asyncSwFqInvGenerator;
        else
            swFreq = freq * vehicle.nModesInvGenerator(p);
        end
        InvLoss = vehicle.INV_RLoss * iInv^2 + vehicle.INV_FLoss * swFreq + vehicle.INV_FRLoss * swFreq * iInv^2;
        pDcLink = (PInverter + InvLoss) * vehicle.nINV; % Total DC-link power in the unit
    else
        swFreq = 0;
        pDcLink = (PInverter + abs(PInverter * vehicle.INV_FixLoss)) * vehicle.nINV;
    end

    % DC-link - pDcLink is the DC-link power per HV system
    pConv0 = (pDcLink + (vehicle.AuxPower / vehicle.APS_efficiency) * 1000 / vehicle.nHVSystems) / vehicle.nCONV + BRPower / vehicle.nCONV;

    % Converter
    if isempty(vehicle.CONV_FixLoss)
        secVoltage = ULine / vehicle.MT_ratio;
        p = -secVoltage^2 / (vehicle.CONV_RLoss + vehicle.CONV_FRLoss * vehicle.swFqConv);
        q = -(pConv0 + vehicle.CONV_FLoss * vehicle.swFqConv) * p;
        pConv = -p / 2 - sqrt((p / 2)^2 - q);
        ConvLoss = pConv - pConv0; %#ok<NASGU>
        IS = pConv / secVoltage;
        PTrafo = pConv * vehicle.nCONV;
    else
        PTrafo = (pConv0 + abs(pConv0 * vehicle.CONV_FixLoss)) * vehicle.nCONV;
        IS = PTrafo * ULine / vehicle.MT_ratio;
    end

    % Converter cable
    if vehicle.R_ConCable > 0
        cable_loss = vehicle.R_ConCable * IS^2;
        PTrafo = (PTrafo + 2 * cable_loss * vehicle.nCONV) * vehicle.nHVSystems;
    else
        PTrafo = PTrafo * vehicle.nHVSystems;
    end

    % Trafo
    if isempty(vehicle.Trafo_FixLoss)
        R_trafo = ((vehicle.TrafoTemp - 20) * 0.00404 + 1) * vehicle.rTrafoWind; % ohm equivalent - re-calc to actual temperature
        TrafoLoss = R_trafo * (PTrafo / ULine)^2 + abs(PTrafo / ULine) * vehicle.trafoIronLoss + baseTrafoLoss;
        pLine = (PTrafo + TrafoLoss); % Total active line power per unit
    else
        pLine = (PTrafo + abs(PTrafo * vehicle.Trafo_FixLoss));

    end

    % Line current
    ILine = pLine / ULine;

    if doCalc == 2 % first loop
        if ILine < iLineLim % we are exceeding the line current limit in braking
            doCalc = 1;
        else
            doCalc = 0;
        end
    end
    if doCalc
        if abs(1 - iLineLim / ILine) < 0.001
            doCalc = 0;
        elseif vehicle.OVCMethod
            BRPower = max([BRPower + (iLineLim - ILine) * ULine, 0]);
        else
            TE = ShSim_maxDBE(max([min([TE * iLineLim / ILine, 0]) DBE_demand]), vehicle, speed, limitation);
        end
    end
end
ULink = vehicle.uDcLink;
PropLoss = ULine * ILine / 1000 - TE * speed - vehicle.AuxPower / vehicle.APS_efficiency - BRPower / 1000;

if vehicle.ES_Op_Mode > 0 && BRPower > 0
    ESPower = -min([BRPower * 1000, vehicle.ES_PowerIn * vehicle.nHVSystems]);
    if limitation.ES_ChargeLevel < 1
        resESTime = vehicle.ES_capacity * (1 - limitation.ES_ChargeLevel) * 3600 / vehicle.ES_PowerIn;
        if resESTime < 1
            ESPower = ESPower * resESTime;
        end
    else
        ESPower = 0;
    end
    BRPower = BRPower + ESPower / 1000;
else
    ESPower = 0;
end
return