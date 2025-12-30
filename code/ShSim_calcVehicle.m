function vcl = ShSim_calcVehicle(vehicle_data)
% vcl = vehicle_data;
% Architecture
vcl.class = vehicle_data.class;
% vcl.Purpose = vehicle_data.Purpose;
vcl.nUnits = vehicle_data.nUnits;
vcl.nHVSystems = vehicle_data.nHVSystems;
vcl.nCONV = vehicle_data.nCONV;
vcl.nAPUs = vehicle_data.nAPUs;
vcl.nINV = vehicle_data.nINV;
vcl.nMotors = vehicle_data.nMotors;
vcl.nAxDrive = vehicle_data.nAxDrive;
vcl.nAxTrail = vehicle_data.nAxTrail;
vcl.notes = vehicle_data.notes;
ws = vehicle_data.wheelsize;
for i1 = 1:size(ws, 2)
    ws{i1} = ws{i1} / 1000; % wheel diameter in metre
end
% Dynamics
switch vehicle_data.MaxSpdUnit
    case 1 % km/h
        spd_scale = 3.6;
    case 2 % mph
        spd_scale = 2.2369;
    case 3 % m/s
        spd_scale = 1;
    case 4 % Hz
        spd_scale = vehicle_data.gearRatio * vehicle_data.nPoles / (ws{vehicle_data.ws_selection} * pi * 2);
end
vcl.vMax = vehicle_data.vMax / spd_scale;
vcl.unitLength = replaceifempty(vehicle_data.vMax, 0);
vcl.max_acc = replaceifempty(vehicle_data.max_acc, inf);
vcl.jerkrateT = replaceifempty(vehicle_data.jerkrateT, inf);
vcl.AccMethod = vehicle_data.AccMethod - 1;
vcl.retard_max = replaceifempty(vehicle_data.retard_max, inf);
vcl.jerkrateB = replaceifempty(vehicle_data.jerkrateB, inf);
if isempty(vehicle_data.retard_min) || isempty(vehicle_data.minRetard_v)
    vcl.retard_min = inf;
    vcl.retard_min_speed = inf;
else
    vcl.retard_min = vehicle_data.retard_min;
    vcl.retard_min_speed = vehicle_data.minRetard_v / spd_scale;
end
vcl.residual_retard = replaceifempty(vehicle_data.residual_ret, 0);
vcl.RetardMethod = vehicle_data.RetardMethod;
vcl.utiliseEBrake = vehicle_data.utilEBrake;
tare = vehicle_data.TareMassUnit;
rot = replaceifempty(vehicle_data.RotMassUnit, 0);
if ~isempty(vehicle_data.TotBrake)
    for i1 = 1:size(vehicle_data.TotBrake)
        vcl.TotBrake_v(i1) = vehicle_data.TotBrake{i1, 1} / spd_scale;
        vcl.TotBrake_f(i1) = vehicle_data.TotBrake{i1, 2};
        vcl.TotBrake_p(i1) = vcl.TotBrake_v(i1) * vcl.TotBrake_f(i1);
    end
else
    vcl.TotBrake_v = [0, vcl.vMax];
    vcl.TotBrake_f = [inf, inf];
    vcl.TotBrake_p = [0, 0];
end

vcl.TareMassUnit = tare;
vcl.RotationMassUnit = rot;
n = size(vehicle_data.LoadMass, 1);
for i1 = 1:n
    vcl.LoadMassUnit{i1} = vehicle_data.LoadMass{i1, 2} / 1000; % change to ton
    vcl.LoadMassName{i1} = vehicle_data.LoadMass{i1, 1};
    vcl.LoadMassDef{i1} = vehicle_data.LoadMass{i1, 3};
end
vcl.LoadCondition = vehicle_data.LoadCond;
vcl.statmassLoad = 0;
vcl.totMass = 0;

vcl.davies_a_actual = replaceifempty(vehicle_data.dvs_a_actual, 0);
vcl.davies_a_dash = replaceifempty(vehicle_data.dvs_a_dash, 0);
vcl.davies_a = 0;
vcl.davies_b1_dash = replaceifempty(vehicle_data.dvs_b1_dash, 0);
vcl.davies_b2 = replaceifempty(vehicle_data.dvs_b2, 0);
vcl.davies_b = 0;
vcl.davies_c_front = replaceifempty(vehicle_data.dvs_c_front, 0);
vcl.davies_c_len = replaceifempty(vehicle_data.dvs_c_len, 0);
vcl.davies_c1 = 0;
vcl.c_front2 = replaceifempty(vehicle_data.dvs_c_front2, 0);
vcl.davies_c2 = 0;
vcl.curve_cr0 = replaceifempty(vehicle_data.curve_cr0, 0);
vcl.curve_cr1 = replaceifempty(vehicle_data.curve_cr1, 0);
vcl.curve_cr2 = replaceifempty(vehicle_data.curve_cr2, 0);

% Traction system
sel = vehicle_data.SpdSelectTraction;
switch sel
    case 1 % km/h
        spd_scale = 3.6;
    case 2 % mph
        spd_scale = 2.2369;
    case 3 % m/s
        spd_scale = 1;
    case 4 % Hz
        spd_scale = vehicle_data.gearRatio * vehicle_data.nPoles * 1000 / (vehicle_data.wheelsize{vehicle_data.ws_selection} * pi * 2);
end

vcl.maxTE = vehicle_data.maxTE;
vcl.maxDBE = vehicle_data.maxDBE;
vcl.maxPwrTEbase = vehicle_data.maxPwrTEbase;
vcl.maxPwrDBEbase = vehicle_data.maxPwrDBEbase;
vcl.maxRatePwrTE = replaceifempty(vehicle_data.maxRatePwrTE, inf);
vcl.maxRatePwrDBE = replaceifempty(vehicle_data.maxRatePwrDBE, inf);

vcl.vBasePwr = vcl.maxPwrTEbase / vcl.maxTE;
vcl.vBasePwr2 = replaceifempty(vehicle_data.vBasePwr2 / spd_scale, inf);
vcl.vBaseBrk = vcl.maxPwrDBEbase / vcl.maxDBE;
vcl.vBaseBrk2 = replaceifempty(vehicle_data.vBaseBrk2 / spd_scale, inf);
vcl.dynBrkFadeStart = replaceifempty(vehicle_data.dynBrkFadeStart / spd_scale, 0);
vcl.dynBrkFadeEnd = replaceifempty(vehicle_data.dynBrkFadeEnd / spd_scale, 0);
if ~isempty(vehicle_data.PwrTE)
    n = size(vehicle_data.PwrTE,1);
    for i1 = 1:n
        vcl.PwrTE_v(i1) = vehicle_data.PwrTE{i1, 1} / spd_scale;
        vcl.PwrTE_f(i1) = vehicle_data.PwrTE{i1, 2};
        vcl.PwrTE_p(i1) = vcl.PwrTE_v(i1) * vcl.PwrTE_f(i1);
    end
else
    vcl.PwrTE_v = [0, vcl.vMax];
    vcl.PwrTE_f = [inf, inf];
    vcl.PwrTE_p = [0, 0];
end
if ~isempty(vehicle_data.PwrDBE)
    n = size(vehicle_data.PwrDBE,1);
    for i1 = 1:n
        vcl.PwrDBE_v(i1) = vehicle_data.PwrDBE{i1, 1} / spd_scale;
        vcl.PwrDBE_f(i1) = vehicle_data.PwrDBE{i1, 2};
        vcl.PwrDBE_p(i1) = vcl.PwrDBE_v(i1) * vcl.PwrDBE_f(i1);
    end
else
    vcl.PwrDBE_v = [0, vcl.vMax];
    vcl.PwrDBE_f = [inf, inf];
    vcl.PwrDBE_p = [0, 0];
end
vcl.customTERef = vehicle_data.customTERef;

% Energy
if vehicle_data.ULSelectionIdx >=3
    vcl.U_Line_Nom = vehicle_data.U_Line_Nom * 1000;
else
    vcl.U_Line_Nom = vehicle_data.U_Line_Nom;
end
if isempty(vehicle_data.uDcLink) || isempty(vehicle_data.uFreq)
    vcl.uDcLink = inf;
    vcl.uFreq = 0;
else
    vcl.uDcLink = vehicle_data.uDcLink;
    vcl.uFreq = vehicle_data.uFreq;
end
vcl.AuxPowerAC = replaceifempty(vehicle_data.AuxPowerAC, 0);
vcl.BACPower = replaceifempty(vehicle_data.BACPower, 0);
vcl.AuxPower = vcl.AuxPowerAC + vcl.BACPower ;
vcl.APS_efficiency = replaceifempty(vehicle_data.APS_efficiency / 100, 1);
vcl.AuxPF = replaceifempty(vehicle_data.AuxPF, 0);
ILDriveMax_train = replaceifempty(vehicle_data.ILDriveMax_train, inf) ;
ILDriveMax_unit = replaceifempty(vehicle_data.ILDriveMax_unit, inf);
ILDriveMax_system = replaceifempty(vehicle_data.ILDriveMax_system, inf);
vcl.ILDriveMax = min([ILDriveMax_train / vcl.nUnits, ILDriveMax_unit, ILDriveMax_system * vcl.nHVSystems]) + eps;
ILBrakeMax_train = replaceifempty(vehicle_data.ILBrakeMax_train, inf);
ILBrakeMax_unit = replaceifempty(vehicle_data.ILBrakeMax_unit, inf);
ILBrakeMax_system = replaceifempty(vehicle_data.ILBrakeMax_system, inf);
vcl.ILBrakeMax = min([ILBrakeMax_train / vcl.nUnits, ILBrakeMax_unit, ILBrakeMax_system * vcl.nHVSystems]) + eps;
vcl.MinUL = replaceifempty(vehicle_data.MinUL, vcl.U_Line_Nom);
vcl.OVClevel = replaceifempty(vehicle_data.OVClevel, vcl.U_Line_Nom);
vcl.OVCMethod = vehicle_data.OVCMethod;
vcl.R_ConCable = replaceifempty(vehicle_data.R_ConCable / 1000, 0);
vcl.rLineInd = replaceifempty(vehicle_data.rLineInd / 1000, 0);
vcl.MT_ratio = replaceifempty(vehicle_data.MT_ratio, 0);
vcl.Trafo_FixLoss = replaceifempty(vehicle_data.Trafo_FixLoss / 100, 0);
vcl.trafoIdleLoss = replaceifempty(vehicle_data.trafoIdleLoss * 1000, 0);
vcl.trafoHarmonicLoss = replaceifempty(vehicle_data.trafoHarmonicLoss * 1000, 0);
vcl.trafoIronLoss = replaceifempty(vehicle_data.trafoIronLoss / 100, 0);
vcl.rTrafoPrim = replaceifempty(vehicle_data.rTrafoPrim / 1000, 0);
vcl.rTrafoSec = replaceifempty(vehicle_data.rTrafoSec, 0) / 1000;
vcl.TrafoTemp = replaceifempty(vehicle_data.TrafoTemp, 0);
vcl.rTrafoWind = (vcl.rTrafoPrim + vcl.rTrafoSec * vcl.MT_ratio / vcl.nCONV) / vcl.nHVSystems;
vcl.ES_Op_Mode = vehicle_data.ES_Op_Mode - 1;
vcl.ES_capacity = replaceifempty(vehicle_data.ES_capacity * 1000, 0);
vcl.ES_PowerOut = replaceifempty(vehicle_data.ES_PowerOut * 1000, 0);
vcl.ES_PowerIn = replaceifempty(vehicle_data.ES_PowerIn * 1000, 0);
vcl.ES_BoostILimit = replaceifempty(vehicle_data.ES_BoostILimit, 0);
vcl.ES_FadeTime = replaceifempty(vehicle_data.ES_FadeTime, 0);

% Drives
vcl.gearRatio = vehicle_data.gearRatio;
vcl.gearLoss0 = replaceifempty(vehicle_data.gearLoss0, 0);
vcl.gearLoss1 = replaceifempty(vehicle_data.gearLoss1, 0);
vcl.gearLoss2 = replaceifempty(vehicle_data.gearLoss2, 0);
vcl.wheelsize = ws;
vcl.ws_mass = vehicle_data.ws_mass;
vcl.ws_rotMass = vehicle_data.ws_rotMass;
for i1 = 1:3
    vcl.wheelsize{i1} = replaceifempty(vcl.wheelsize{i1}, vcl.wheelsize{1});
    vcl.ws_mass{i1} = replaceifempty(vcl.ws_mass{i1}, 0);
    vcl.ws_rotMass{i1} = replaceifempty(vcl.ws_rotMass{i1}, 0);
end
vcl.ws_selection = vehicle_data.ws_selection;
vcl.hasMotorMap = vehicle_data.hasMotorMap;
vcl.motorMethod = vehicle_data.motorMethod;
vcl.motorMap = vehicle_data.motorMap;
vcl.nPoles = vehicle_data.nPoles;
vcl.TM_FixLoss = replaceifempty(vehicle_data.TM_FixLoss, 0) / 100;

if vcl.hasMotorMap
    vcl.Speed0 = [vcl.motorMap.rpmBaseMotor, vcl.motorMap.rpmBaseGenerator] / vcl.gearRatio * vcl.wheelsize{vcl.ws_selection} * pi;
else
    vcl.Speed0=[vcl.vBasePwr, vcl.vBaseBrk];
end
vcl.R_MotCable  = replaceifempty(vehicle_data.R_MotCable , 0) / 1000;

% Converters
vcl.swFqConv = replaceifempty(vehicle_data.swFqConv, 0);
vcl.CONV_RLoss = replaceifempty(vehicle_data.CONV_RLoss, 0) / 1000;
vcl.CONV_FLoss = replaceifempty(vehicle_data.CONV_FLoss, 0);
vcl.CONV_FRLoss = replaceifempty(vehicle_data.CONV_FRLoss, 0) / 1e6;
vcl.CONV_FixLoss = vehicle_data.CONV_FixLoss / 100;

vcl.INV_RLoss = replaceifempty(vehicle_data.INV_RLoss, 0) / 1000;
vcl.INV_FLoss = replaceifempty(vehicle_data.INV_FLoss, 0);
vcl.INV_FRLoss = replaceifempty(vehicle_data.INV_FRLoss, 0) / 1e6;
vcl.INV_FixLoss = vehicle_data.INV_FixLoss / 100;
vcl.asyncSwFqInvMotor = replaceifempty(vehicle_data.asyncSwFqInvMotor, 0);
vcl.asyncSwFqInvGenerator = replaceifempty(vehicle_data.asyncSwFqInvGenerator, 0);
n = size(vehicle_data.MotorPulseData, 1);
if n > 1
    vcl.nModesInvMotor(n - 1, 1) = 0;
    vcl.FreqInvMotor(n - 1, 1) = 0;
    for i1 = 2:n
        vcl.nModesInvMotor(i1 - 1, 1) = str2double(vehicle_data.MotorPulseData{i1, 1});
        vcl.FreqInvMotor(i1 - 1, 1) = vehicle_data.MotorPulseData{i1, 2};
    end
else
    vcl.nModesInvMotor = [];
    vcl.FreqInvMotor = [];
end
% vcl.nModesInvMotor = [vehicle_data.MotorPulseData{:, 1}];
% vcl.FreqInvMotor = [];
n = size(vehicle_data.GeneratorPulseData, 1);
if n > 1
    vcl.nModesInvGenerator(n - 1, 1) = 0;
    vcl.FreqInvGenerator(n - 1, 1) = 0;
    for i1 = 2:n
        vcl.nModesInvGenerator(i1 - 1, 1) = str2double(vehicle_data.GeneratorPulseData{i1, 1});
        vcl.FreqInvGenerator(i1 - 1, 1) = vehicle_data.GeneratorPulseData{i1, 2};
    end
else
    vcl.nModesInvGenerator = [];
    vcl.FreqInvGenerator = [];
end

% vcl.nModesInvGenerator = vehicle_data.GeneratorPulseData;
% vcl.FreqInvGenerator = [];

% Notching
vcl.notchMethod = vehicle_data.notchMethod;
vcl.nNotchT = replaceifempty(vehicle_data.nNotchT, 0);
vcl.nNotchB = replaceifempty(vehicle_data.nNotchB, 0);
vcl.notchMinTime = replaceifempty(vehicle_data.notchMinTime, 0);
vcl.notchrate = replaceifempty(vehicle_data.notchrate, 0);
switch vehicle_data.notchSpdSel
    case 1 % km/h
        spd_scale = 3.6;
    case 2 % mph
        spd_scale = 2.2369;
    case 3 % m/s
        spd_scale = 1;
end
if size(vehicle_data.notchDeltaV, 2) < 2
    vcl.notchDeltaV = [0, 0];
else
    vcl.notchDeltaV = vehicle_data.notchDeltaV / spd_scale;
end
if size(vehicle_data.notchMinV, 2) < 2
    vcl.notchMinV = [0, 0];
else
    vcl.notchMinV = vehicle_data.notchMinV / spd_scale;
end
vcl.notchDeltaVAcc = replaceifempty(vehicle_data.notchDeltaVAcc, 0) / spd_scale;
vcl.notchMaxAcc = replaceifempty(vehicle_data.notchMaxAcc, 0);

% Thermal parameters
vcl.AmbientTemp = vehicle_data.AmbientTemp;
vcl.IPrimTrafoCont = replaceifempty(vehicle_data.IPrimTrafoCont, inf);
vcl.IPrimTrafoContTC = replaceifempty(vehicle_data.IPrimTrafoContTC, 3600);

vcl.I2ndTrafoCont = replaceifempty(vehicle_data.I2ndTrafoCont, inf);
vcl.I2ndTrafoContTC = replaceifempty(vehicle_data.I2ndTrafoContTC, 3600);
vcl.IDCFilterCont = replaceifempty(vehicle_data.IDCFilterCont, inf);
vcl.IDCFilterContTC = replaceifempty(vehicle_data.IDCFilterContTC, 3600);
vcl.IphMCMMax = replaceifempty(vehicle_data.IphMCMMax, inf);
vcl.IphMCMContSin = replaceifempty(vehicle_data.IphMCMContSin, inf);
vcl.IphMCMContHex = replaceifempty(vehicle_data.IphMCMContHex, inf);
vcl.IphMCMContTC = replaceifempty(vehicle_data.IphMCMContTC, 3600);
vcl.IphLCMMax = replaceifempty(vehicle_data.IphLCMMax, inf);
vcl.IphLCMCont = replaceifempty(vehicle_data.IphLCMCont, inf);
vcl.IphLCMContTC = replaceifempty(vehicle_data.IphLCMContTC, 3600);
vcl.IstatTMMax = replaceifempty(vehicle_data.IstatTMMax, inf);
vcl.IstatTMCont = replaceifempty(vehicle_data.IstatTMCont, inf);
vcl.IstatTMContTC = replaceifempty(vehicle_data.IstatTMContTC, 3600);
vcl.TMventType = vehicle_data.TMventType;
vcl.BRThermC = replaceifempty(vehicle_data.BRThermC, inf);
vcl.BRThermTC_FanOnLow = replaceifempty(vehicle_data.BRThermTC_FanOnLow, 3600);
vcl.BRThermTC_FanOnHigh = replaceifempty(vehicle_data.BRThermTC_FanOnHigh, 3600);
vcl.BRThermTC_FanOff = replaceifempty(vehicle_data.BRThermTC_FanOff, 3600);
vcl.BRThermTC_Tempfactor = replaceifempty(vehicle_data.BRThermTC_Tempfactor, 100) / 100;
vcl.BRFanOnTemp = replaceifempty(vehicle_data.BRFanOnTemp, inf);
vcl.BRFanHighSpeedTemp = replaceifempty(vehicle_data.BRFanHighSpeedTemp, inf);
vcl.BRFanOffTemp = replaceifempty(vehicle_data.BRFanOffTemp, inf);

switch vehicle_data.BRFanSpeedSel
    case 1 % km/h
        spd_scale = 3.6;
    case 2 % mph
        spd_scale = 2.2369;
    case 3 % m/s
        spd_scale = 1;
    case 4 % Hz
        spd_scale = vehicle_data.gearRatio * vehicle_data.nPoles/(vehicle_data.wheelsize{vehicle_data.ws_selection} * pi * 2);
end
vcl.BRFanOffSpeed = replaceifempty(vehicle_data.BRFanOffSpeed, inf) / spd_scale;
vcl.BRTemp_High = replaceifempty(vehicle_data.BRTemp_High, inf);
vcl.BRTemp_Max = replaceifempty(vehicle_data.BRTemp_Max, inf);


return % recalcfields

function val = replaceifempty(fld, rpl)
if isempty(fld)
    val = rpl;
else
    val = fld;
end
return % replaceifempty