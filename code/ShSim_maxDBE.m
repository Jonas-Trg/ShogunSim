function TE = ShSim_maxDBE(TE_demand, vehicle, speed, limitation)
% ShSim_maxDBE
% 0.1.2.1
%
%
% Ver       Date        Sign    Descr
% 0.1.2.1   2025-12-18  JR      Fix for vehicle.PwrDBE_p == inf.
%

if any(vehicle.PwrDBE_p == inf)
    TE = TE_demand;
else % below only with valid user defined brake curve
    if speed <= vehicle.PwrDBE_v(1)         % if speed is less than the first speed point (=lower fading speed)
        TE = -vehicle.PwrDBE_f(1);           % No Power
        if speed > vehicle.dynBrkFadeStart
            TE = -vehicle.maxDBE;
        end
    elseif speed < vehicle.PwrDBE_v(2)      %if speed is less than 2nd speed point
        d = speed / vehicle.PwrDBE_v(2);    %linear effort
        TE = -d * vehicle.PwrDBE_f(1) + (d - 1) * vehicle.PwrDBE_f(2);
    else
        x = find(vehicle.PwrDBE_v > speed, 1, 'first') - 1;
        if isempty(x)
            x = length(vehicle.PwrDBE_v) - 1;
        end
        d = (speed - vehicle.PwrDBE_v(x)) / (vehicle.PwrDBE_v(x + 1) - vehicle.PwrDBE_v(x));
        P = d * vehicle.PwrDBE_p(x + 1) + (1 - d) * vehicle.PwrDBE_p(x);
        TE = -P / speed;
    end
end
if speed < vehicle.dynBrkFadeStart              %If the speed is less than higher fading speed
    d = (speed - vehicle.dynBrkFadeEnd) / (vehicle.dynBrkFadeStart - vehicle.dynBrkFadeEnd);    %linear ramp down to lower fading speed
    TE = max([TE, d * -vehicle.maxDBE]);
end
TE_demand = max([TE_demand, -vehicle.maxDBE]);
if TE < TE_demand %ensure that we don't request more than demand.
    TE = TE_demand;
end
if TE > 0
    TE = 0;
end

if ~isempty(char(vehicle.customTERef))
    if exist(char(vehicle.customTERef), 'file')
        DBE_user = eval([char(vehicle.customTERef) '(TE,vehicle,speed)']);
        if DBE_user > TE_demand % if user DBE is greater (less negative) then use this.
            TE = DBE_user;
        else
            TE = TE_demand;
        end
        if TE > 0
            TE = 0;
        end
    end
end
if speed > 0
    p = min([vehicle.maxPwrDBEbase * vehicle.vBaseBrk2 / speed, vehicle.maxPwrDBEbase]);
    TE=max([TE, -p * limitation.thermEffort * limitation.thermPower / speed]);
end
return %ShSim_maxDBE
