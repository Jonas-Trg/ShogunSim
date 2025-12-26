function [F_target,rampStart,rampEnd,mode]=ShSim_TargetBE(TE_old,Acc_old,speed,v_target,tvehicle,F_slope,F_res,dT)
% ShSim_TargetBE
% 0.1.1.0
% This function calculates the target brake force to be used when reducing
% speed. It is defined by minimum of:
%   1. retard_max - the target retardation defined by the vehicle file and its
%                   method BrKinMethod.
%   2. retard_min - combined with retard_min_speed, which is used if "Optimise
%                   DynBr" option is selected
%   3. Total brake curve
%
F_target=ShSim_maxBE(tvehicle,speed,F_slope,F_res);
%Include the jerk rate
if speed>=v_target
    if 1
        a_target=-(F_res+F_slope)/tvehicle.totMass; %target acceleration in new point
        speed=v_target-a_target*dT;
    else
        speed=v_target;
        a_target=0;
    end
    rampStart=0;
    rampEnd=0;
    mode=0;
else %speed<v_target
    
    t_rampDown=-Acc_old/tvehicle.jerkrateB; %the time it takes to ramp down current acceleration
    if isinf(tvehicle.jerkrateB)
        v_rampDown=v_target;
    else
        v_rampDown=v_target-tvehicle.jerkrateB*t_rampDown^2/2;    %Speed where ramp down shall start
    end
    if sign(Acc_old)<0
        t_2rampDown=min(max((speed-v_rampDown)/Acc_old,0),dT);  %time to where ramp down shall start
    else
        t_2rampDown=inf;
    end
    if t_2rampDown<dT %it will happen within this sample
        mode=-1;    %Mode = -1 indicates jerk rate ramp down
        rampStart=t_2rampDown/dT;
        if t_rampDown<(dT-t_2rampDown)
            rampEnd=min((t_2rampDown+t_rampDown)/dT,1);
            a_target=0;
        else
            rampEnd=1;
            a_target=-(sqrt((v_target-speed)*2/tvehicle.jerkrateB)-dT)*tvehicle.jerkrateB;
            %                 a_target=Acc_old+(dT-t_2rampDown)*tvehicle.jerkrateB;
        end
    else
        mode=0;
        rampStart=0;
        rampEnd=0;
        a_target=(speed-v_target)/dT;
    end
end
F_target = max([a_target*tvehicle.totMass+F_res+F_slope F_target]);
F_target = min([F_target 0]);
if rampEnd==1 && rampStart==0
    
elseif F_target>TE_old+tvehicle.dDBE_jerkLimit
    F_target=TE_old+tvehicle.dDBE_jerkLimit;
    mode=-1;
    rampStart=0;
    rampEnd=1;
elseif F_target<TE_old-tvehicle.dDBE_jerkLimit
    F_target=TE_old-tvehicle.dDBE_jerkLimit;
    mode=1;
    rampStart=0;
    rampEnd=1;
end
F_target = min([0 F_target tvehicle.totMass*(v_target-speed)/dT + F_res + F_slope]);
if F_target<TE_old %ramping up
    mode=1;
    rampStart=0;
    rampEnd=(TE_old-F_target)/tvehicle.dDBE_jerkLimit;
end
return %ShSim_TargetBE
