function r=ShSim_sum(in)
% ShSim_sum
% 0.1.1.0
n=size(in,2);
r=in;
temp=r(end);
temp.timedata = [];
r(n+1)=temp;
r(n+1).depart=r(1).depart;
r(n+1).from_station=r(1).from_station;
for k=1:n-1 %it should be n-1 here. last record is already in r(end)
    r(n+1).runtime           = r(n+1).runtime           + r(k).runtime;
    r(n+1).distance          = r(n+1).distance          + r(k).distance;
    r(n+1).w_traction        = r(n+1).w_traction        + r(k).w_traction;
    r(n+1).w_regen           = r(n+1).w_regen           + r(k).w_regen;
    r(n+1).w_Aux             = r(n+1).w_Aux             + r(k).w_Aux;
    r(n+1).w_railPos         = r(n+1).w_railPos         + r(k).w_railPos;
    r(n+1).w_railNeg         = r(n+1).w_railNeg         + r(k).w_railNeg;
    r(n+1).w_resistance      = r(n+1).w_resistance      + r(k).w_resistance;
    r(n+1).w_friction        = r(n+1).w_friction        + r(k).w_friction;
    r(n+1).w_friction2       = r(n+1).w_friction2       + r(k).w_friction2;
    r(n+1).w_PropLoss        = r(n+1).w_PropLoss        + r(k).w_PropLoss;
    r(n+1).i_rms_sqsum       = r(n+1).i_rms_sqsum       + r(k).i_rms_sqsum;
    r(n+1).i_rms_num         = r(n+1).i_rms_num         + r(k).i_rms_num;
    r(n+1).i_rms2_sqsum      = r(n+1).i_rms2_sqsum      + r(k).i_rms2_sqsum;
    r(n+1).i_rms2_num        = r(n+1).i_rms2_num        + r(k).i_rms2_num;
    r(n+1).t_coast           = r(n+1).t_coast           + r(k).t_coast;
    r(n+1).d_coast           = r(n+1).d_coast           + r(k).d_coast;
    r(n+1).ThermRedEffortSum = r(n+1).ThermRedEffortSum + r(k).ThermRedEffortSum;
    r(n+1).ThermRedPowerSum  = r(n+1).ThermRedPowerSum  + r(k).ThermRedPowerSum;
    r(n+1).nThermRecords     = r(n+1).nThermRecords     + r(k).nThermRecords;
end
return % ShSim_sum
