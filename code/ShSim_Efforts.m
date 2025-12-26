function [TE_achieved, ULine, ULink, ILine, BRPower, Iph, IS, mSlip, fMode, ESPower, PropLoss]=...
    ShSim_Efforts(TE_demand, vehicle, speed, routedata, sim_config, limitations)
%ShSim_Efforts
% 0.1.1.0
%
%   Description
%       This function is called several times for every sample when
%       executing the route calculation.
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
%           vehicle uNorm
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
%   
if vehicle.ES_Op_Mode>0
    limitations.ES_ChargeLevel=1;
else
    limitations.ES_ChargeLevel=0;
end
if sim_config.ThermReduction==0 || ~isfield(limitations,'thermPower')
    limitations.thermPower=1;
    limitations.thermEffort=1;
end
if vehicle.uFreq==0 %DC mode
    if TE_demand>=0
        [TE_achieved, ULine, ULink, ILine, BRPower, Iph, IS, mSlip, fMode, ESPower, PropLoss]=...
            ShSim_DC_TE(TE_demand, vehicle, speed, routedata, sim_config, limitations);
    else
        [TE_achieved,ULine,ULink,ILine,BRPower,Iph,IS,mSlip,fMode,ESPower,PropLoss]=...
            ShSim_DC_DBE(TE_demand, vehicle, speed, routedata, sim_config, limitations);
    end
else
    if TE_demand>=0
        [TE_achieved,ULine,ULink,ILine,BRPower,Iph,IS,mSlip,fMode,ESPower,PropLoss]=ShSim_AC_TE(TE_demand, vehicle, speed, routedata, sim_config, limitations);
    else
        [TE_achieved,ULine,ULink,ILine,BRPower,Iph,IS,mSlip,fMode,ESPower,PropLoss]=ShSim_AC_DBE(TE_demand, vehicle, speed, routedata, sim_config, limitations);
    end
end
return % ShSim_Efforts