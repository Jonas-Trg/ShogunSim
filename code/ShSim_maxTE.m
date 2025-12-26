function TE = ShSim_maxTE(TE_demand, vehicle, speed, limitation)
% ShSim_maxTE
% 0.1.2.0
%   Description
%       Gives the maximum effort at a given speed limited to TE_demand
%   Input:
%       TE_demand - reference for rail effort in kN
%       vehicle, struct of vehicle definition
%       speed, in m/s
%       limitation, struct of
%           thermEffort
%           thermPower
%   Output:
%       TE - maximum output effort on rail including gear loss
%
if speed <= vehicle.PwrTE_v(1)
    TE = vehicle.PwrTE_f(1);
elseif speed < vehicle.PwrTE_v(2)
    d = speed / vehicle.PwrTE_v(2);
    TE = d * vehicle.PwrTE_f(2) + (1 - d) * vehicle.PwrTE_f(1);
else
    x = find(vehicle.PwrTE_v > speed, 1, 'first')-1;
    if isempty(x)
        x = length(vehicle.PwrTE_v) - 1;
    end
    d = (speed - vehicle.PwrTE_v(x)) / (vehicle.PwrTE_v(x + 1) - vehicle.PwrTE_v(x));
    P = d * vehicle.PwrTE_p(x + 1) + (1 - d) * vehicle.PwrTE_p(x);
    TE = P / speed;
end
TE = TE * limitation.thermEffort * limitation.thermPower;
if speed > 0
    p = min([vehicle.maxPwrTEbase * vehicle.vBasePwr2 / speed, vehicle.maxPwrTEbase]);
    TE = min([p * limitation.thermPower / speed, vehicle.maxTE * limitation.thermEffort, TE, TE_demand]);
else
    TE = min([vehicle.maxTE, TE, TE_demand]);
end

if ~isempty(char(vehicle.customTERef))
    if exist(char(vehicle.customTERef), 'file')
        TE_user = eval([char(vehicle.customTERef) '(TE,vehicle,speed)']);
        if TE_user < TE
            TE = TE_user;
        end
    end
end
return %ShSim_maxTE
