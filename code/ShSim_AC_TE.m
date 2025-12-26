function [TE, ULine, ULink, ILine, BRPower, iph, IS, mSlip, swFreq, ESPower, PropLoss] =...
    ShSim_AC_TE(TE_demand, vehicle, speed, routedata, sim_config, limitation)
% ShSim_AC_TE
% 0.1.1.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   AC Traction mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Description
%       This function is called several times for every sample in traction
%       or coasting mode in AC mode operation.
%
%   Input data
%       TE_demand - effort reference in kN per unit
%       vehicle - struct containing at least:
%       speed - in m/s
%       routedata*, if used it shall contains:
%           USupply - Line voltage definition from line
%           * Can be left empty, i.e. routedata=[];
%       simulation
%           useSubStations - defines whether to use voltage from line supply or
%                           vehicle.U_Line_Nom
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
%       ESPower - power used in the Energy Saver (not used in AC mode at
%                 the moment)
%       PropLoss - The actual loss in the propulsion system - total on train level
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-09  JR      New coding design (spaces) plus motor loss in Watts
% 0.1.1.2   2025-12-22  JR      Removed 'hasEfficiency'
%
ESPower = 0;
BRPower = 0;
iph = 0;
mSlip = 0;

if isempty(routedata)
    if vehicle.MinUL > 0
        ULine = vehicle.MinUL;
    else
        ULine = vehicle.U_Line_Nom;
    end
else
    if sim_config.useSubStations
        ULine = routedata.USupply;
    elseif vehicle.MinUL > 0
        ULine = vehicle.MinUL;
    else
        ULine = vehicle.U_Line_Nom;
    end
end

if isempty(vehicle.Trafo_FixLoss)
    baseTrafoLoss = (vehicle.trafoIdleLoss + vehicle.trafoHarmonicLoss) * vehicle.nHVSystems;
    if speed<1 && TE_demand == 0 && 0
        baseTrafoLoss = vehicle.trafoIdleLoss; %#ok<UNRCH>
    end
end


TE = ShSim_maxTE(TE_demand, vehicle, speed, limitation);

doCalc = 2; % 2 = first loop, 1 = tuning iline from above, -1 = tuning iline from below
iLineLim = vehicle.ILDriveMax;
ILine = inf;
while doCalc ~= 0
    IL_old = ILine; % save previous line current
    if TE == 0 && speed == 0 % Inverter is OFF
        pDcLink = 0;
        swFreq = 0;
    else
        % Drives and Motor
        nMotorRPM = vehicle.gearRatio*speed / (vehicle.wheelsize{vehicle.ws_selection} * pi) * 60; % revs per min

        % tq is the torque per motor
        tq = TE * vehicle.wheelsize{vehicle.ws_selection} * 1000 / (vehicle.gearRatio * vehicle.nAxDrive * 2);
        
        % Gearbox losses in Nm
        tq = tq + vehicle.gearLoss0 + abs(tq * vehicle.gearLoss1) + nMotorRPM * vehicle.gearLoss2;

        switch vehicle.motorMethod % 1 = fixed efficiency, 2 = motor map
            case 1 % fixed efficiency
                PInverter = tq * nMotorRPM * 2 * pi / 60;
                PInverter = PInverter + abs(PInverter * vehicle.TM_FixLoss);
                uph = min([(speed + 1) / vehicle.vBasePwr, 1]) * vehicle.uDcLink * sqrt(6) / pi;
                iph = TE * speed * 1000 / (sqrt(3) * uph * 0.9) / vehicle.nAxDrive;
                freq = nMotorRPM * vehicle.nPoles / 120;
                
            case 2 % motor map
                [uph, iph, mloss, mSlip, pf] = ShSim_MotorExe(nMotorRPM, tq, vehicle.motorMap); %#ok<ASGLU>
                PInverter = nMotorRPM .* tq / 60 * 2 * pi + mloss;
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
            p = find(freq > vehicle.FreqInvMotor, 1, 'last');
            if isempty(p)
                swFreq = vehicle.asyncSwFqInvMotor;
            else
                swFreq = freq * vehicle.nModesInvMotor(p);
            end
            InvLoss = vehicle.INV_RLoss * iInv^2 + vehicle.INV_FLoss * swFreq + vehicle.INV_FRLoss * swFreq * iInv^2;
            pDcLink = (PInverter + InvLoss) * vehicle.nINV; % Total DC-link power in the unit
        else
            swFreq = 0;
            pDcLink = (PInverter + abs(PInverter * vehicle.INV_FixLoss)) * vehicle.nINV;
        end
    end

    % DC-link - pDcLink is the DC-link power per HV system
    pConv0 = (pDcLink + vehicle.AuxPower / vehicle.APS_efficiency * 1000 / vehicle.nHVSystems) / vehicle.nCONV;

    % Converter
    if isempty(vehicle.CONV_FixLoss)
        secVoltage = ULine / vehicle.MT_ratio;
        p = -secVoltage^2 / (vehicle.CONV_RLoss + vehicle.CONV_FRLoss * vehicle.swFqConv);
        q = -(pConv0 + vehicle.CONV_FLoss * vehicle.swFqConv) * p;
        pConv = -p / 2 -sqrt((p / 2)^2 - q);
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

    % Transformer
    if isempty(vehicle.Trafo_FixLoss)
        R_trafo = ((vehicle.TrafoTemp - 20) * 0.00404 + 1) * vehicle.rTrafoWind; % ohm equivalent - re-calc to actual temperature
        TrafoLoss = R_trafo * (PTrafo / ULine)^2 + abs(PTrafo / ULine) * vehicle.trafoIronLoss + baseTrafoLoss;
        pLine = (PTrafo + TrafoLoss); % Total active line power per unit
    else
        pLine = (PTrafo + abs(PTrafo * vehicle.Trafo_FixLoss)); % test
    end

    ILine = pLine / ULine;

    if doCalc == 2 % first loop
        if ILine > iLineLim
            doCalc = 0.1;
        else
            doCalc = 0;
        end
    elseif doCalc > 0 %we are in reducing mode
        if ILine < iLineLim
            doCalc = -doCalc / 2;
        elseif IL_old - ILine < 0.1
            doCalc = 0;
        end
    else
        if ILine > iLineLim
            doCalc =- doCalc / 2;
        end
    end
    if abs(doCalc) < 0.001
        doCalc = 0;
    else
        TE = ShSim_maxTE(TE * (1 - doCalc), vehicle, speed, limitation);
    end
end

ULink = vehicle.uDcLink;
PropLoss = ULine * ILine / 1000 - TE * speed - vehicle.AuxPower / vehicle.APS_efficiency;
return % ShSim_AC_TE