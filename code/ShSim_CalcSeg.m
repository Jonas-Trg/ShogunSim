function [r, stop_time, thermal_out, v_cut, d_cut] = ...
    ShSim_CalcSeg(trackdata, vehicle, sim_config, start_time, thermal_in, routedata, slopedata, data, hwb)
% ShSim_CalcSeg
% 0.1.1.1
%
%   Works with all physics on one unit.
% Input signals
%   trackdata
%   vehicle
%   simulation
%   start_time
%   thermal_in
%   routedata
%   slopedata
%   do_echo
% Output signals
%   r           Result data record
%   stop_time   Time of last record (starting time for next section
%   thermal_out 
%   v_cut
%   d_cut
% Description
%   A section from station to station is calculated here. All Efforts and 
%`  masses are per train unit.
%timedata record
% 01 Time
% 02 Dist
% 03 Speed
% 04 TE/DBE
% 05 FricBE
% 06 ULine
% 07 ILine
% 08 BRPow
% 09 BRTemp
% 10 iMotor
% 11 mSlip
% 12 Switching Frquency
% 13 ULink
% 14 ESPow - Energy Saver Power
% 15 Main Transformer prim winding RMS
% 16 Main Transformer 2nd winding RMS
% 17 Line filter RMS
% 18 INV phase current RMS
% 19 CONV phase current RMS
% 20 Traction motor RMS
% 21 Power Reduction factor
% 22 Effort Reduction factor
% 23 Running Resistance
% 24 Traction loss
% 25 Slope force
% 26 Average speed
% 27 Average TE/DBE
% 28 Average FricBE
% 29 Average line current
% 30 Tunnel (effective C-factor of tunnel)
% 31 Notch
% 32 Curve resistance
% 33 Tunnel resistance
% 34 INV phase current
% 35 CONV phase current
% 
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-26  JR      New coding design (spaces)

do_echo = data.echo;
if do_echo
    fprintf(['Calculating part: ', trackdata.from_station, '-', trackdata.to_station, '\n']);
    if ishandle(hwb)
        waitbar(0, hwb, [trackdata.from_station, '-', trackdata.to_station]);
    end
end
tic;
result.fName = '';
result.runtime = 0;
result.depart = 0;
result.arrival = 0;

section_index = trackdata.from_track;
dT = sim_config.delta_T;
[tvehicle, tRoute, error, op_mode] = check_parameters(routedata ,routedata.track(section_index, 9:end), vehicle,dT, do_echo, vehicle.class);
mdone = error;
    
dist = routedata.track(section_index, 1) * 1000;
v_cut = 0;
d_cut = 0;
V_new = 0;
index = 1;
if isnan(start_time)
    start_time = 0;
end
time = start_time;

limitation.thermPower = 1;
limitation.thermEffort = 1;
limitation.ES_ChargeLevel = 0;
limitation.IL_red = 1;
limitation.I2nd_red = 1;
limitation.ILC_red = 1;
limitation.IDCF_ref = 1;
limitation.ITM_stat = 1;
limitation.IMC_stat = 1;

[TE_average, ULine, ULink, IL_new, BRPow, iMotor, IS, mSlip, fMode, ESPow, PropLoss] = ...
    ShSim_Efforts(0, tvehicle, V_new, tRoute, sim_config, limitation); 
if tvehicle.ES_Op_Mode > 0
    limitation.ES_ChargeLevel = 1;
end
timedata=zeros(1000, 50);
timedata(index, 1)  = time;            % Time
timedata(index, 2)  = dist;            % Dist
timedata(index, 3)  = V_new;           % Speed
timedata(index, 4)  = TE_average;      % TE/DBE
timedata(index, 5)  = 0;               % FBE
timedata(index, 6)  = ULine;           % ULine
timedata(index, 7)  = IL_new;          % ILine
timedata(index, 8)  = BRPow;           % BR Power
BRTemp = thermal_in.BRTemp;
timedata(index, 9)  = BRTemp;          % BR temp
timedata(index, 10) = iMotor;          % Motor stator current
timedata(index, 11) = mSlip;           % Motor slip (Hz)
timedata(index, 12) = fMode;           % Modulation index
timedata(index, 13) = ULink;           % DC-link voltage
timedata(index, 14) = ESPow;           % Energy saver power
if tvehicle.uFreq > 0
    timedata(index, 15) = thermal_in.iPrimTrafoRMS^2; % Main Transformer prim winding RMS
    timedata(index, 16) = thermal_in.i2ndTrafoRMS^2;  % Main Transformer 2nd winding RMS
    timedata(index, 17) = 0;                          % Line filter RMS
    timedata(index, 19) = thermal_in.iPhLCM_RMS^2;    % CONV phase current RMS
else
    timedata(index, 15) = 0;                          % Main Transformer prim winding RMS
    timedata(index, 16) = 0;                          % Main Transformer 2nd winding RMS
    timedata(index, 17) = thermal_in.iLineFiltRMS^2;  % Line filter RMS
    timedata(index, 19) = 0;                          % CONV phase current RMS
end
timedata(index, 18) = thermal_in.iPhMCM_RMS^2;        % INV phase current RMS
timedata(index, 20) = thermal_in.iStatorTM_RMS^2;     % Traction motor RMS
timedata(index, 21) = 1;                              % Power Reduction factor
timedata(index, 22) = 1;                              % Effort Reduction factor
timedata(index, 23) = 0;                              % Train resistance
timedata(index, 24) = PropLoss;                       % Prop loss, power
timedata(index, 25) = get_slope(dist, slopedata)*tvehicle.statmassLoad * 9.81/1000;
timedata(index, 26:49) = 0;                           % Allocate for future use
timedata(index, 29) = IL_new;                         % Average line current
if tvehicle.uFreq > 0
    timedata(index, 35) = IS;                         % Average line current
end
timedata(index, 50) = op_mode;

endSection = routedata.track(section_index + 1, 1) * 1000;
last_section_index = trackdata.to_track; % size(routedata.track, 1);
% fprintf(['Calculating section 1 of ' num2str(last_section_index) '\n']);
TE_new = 0;
BE_new = 0;
Acc_new = 0;
start_clock = datetime;
tic;
time = time + dT;
doCoast = 0;
notch = 0;
if tvehicle.notchMethod > 0
    lastNotch = time - tvehicle.notchrate - dT;
end

while (~mdone) && (time < sim_config.maxT + start_time)
    if toc > 0.2
        if do_echo
            if ishandle(hwb)
                waitbar((dist / 1000 - routedata.track(1, 1)) / (routedata.track(end, 1) - routedata.track(1, 1)), hwb);
            else % user has closed the waitbar figure
                break;
            end
        end
        if seconds(datetime - start_clock) > sim_config.time_out
            dlg = warndlg('Calculation interrupted by time-out', 'Time-out');
            uiwait(dlg);
            break
        end
        tic;
    end
    % ================= New Methodology ==================
    % Step1 - speed, dist are actuals from last calculated point
    % Get resistances from last point
    V_old = V_new;
    BE_old = BE_new;
    TE_old = TE_new;
    Acc_old = Acc_new;
    [force_resistance, force_curve, force_tunnel] = getResistances(tvehicle, V_old, dist, routedata.track(section_index,:), tRoute.Track_Guage, tRoute.headWind);
    % Get slope force
    F_slope = get_slope(dist, slopedata) * tvehicle.statmassLoad * 9.81 / 1000;
    % Target speed
    v_target = min([routedata.track(section_index, 2), tvehicle.vMax, tRoute.ATO_SpdMax / 3.6]) - tRoute.ATO_SpdMgn / 3.6;
    % Target force in new point
    [F_target, T_start, T_end, mode] = GetTargetTractionForce(TE_old, Acc_old, V_old, v_target, tvehicle, F_slope, force_resistance, dT, limitation);
    if tvehicle.notchMethod > 0
        % Notch control
        [F_target, notch, lastNotch] = Notching(time, notch, lastNotch, TE_old,V_old, v_target, Acc_old, tvehicle, mode, F_slope, force_resistance, dT, limitation);
    else
        % ATO control
        [doCoast, F_target] = ATO_control(doCoast, V_old, F_target, tRoute, dist, tvehicle, F_slope, do_echo, dT);
        notch = sign(F_target);
    end
    % Calc traction
    [force_traction0, ULine, ULink, ILine0, BRPow, iMotor, IS, mSlip, fMode, ESPow, PropLoss] = ...
        ShSim_Efforts(F_target, tvehicle, V_old, tRoute, sim_config, limitation); %#ok<ASGLU>
    % Calc friction brake - only if F_target>force_traction (can't hold speed down hill)
    BE_new = max([force_traction0 - F_target, 0]);
    % First average calc
    TE_average = (T_start + T_end) / 2 * TE_old + (1 - (T_start + T_end) / 2) * force_traction0;
    % test calculation average acceleration (AveAcc1)
    acc0 = (TE_average - force_resistance - BE_new - F_slope) / tvehicle.totMass;
    tSpeed = max([V_old + acc0 * dT, 0]); % Speed1
    % Re-calculate New point
    [TE_new, ~, ~, IL_new] = ShSim_Efforts(F_target, tvehicle, tSpeed, tRoute, sim_config, limitation);
    % Calc friction brake - only if F_target>force_traction (can't hold speed down hill)
    BE_new = max([force_traction0 - F_target, 0]);
    BE_average = (T_start + T_end) / 2 * BE_old + (1 - (T_start + T_end) / 2) * BE_new;
    % Second average calc
    TE_average = (T_start + T_end) / 2 * TE_old + (1 - (T_start + T_end) / 2) * TE_new;
    Acc_new = (TE_new - force_resistance - BE_new - F_slope) / tvehicle.totMass;
    A_average = (TE_average - force_resistance - BE_new - F_slope) / tvehicle.totMass;
    V_new = max([V_old + A_average * dT, 0]); % Speed1
    % Check this - not entirely correct, but maybe good enough
    V_average = (V_old + V_new) / 2;
%     IL_average=(IL_old + IL_new)/2;
    [TE_average, ULine, ULink, IL_average, BRPow, iMotor, IS, mSlip, fMode, ESPow, PropLoss] = ...
        ShSim_Efforts(TE_average, tvehicle, V_average, tRoute, sim_config, limitation);

%     IL_average=(T_start+T_end)/2*IL_old+(1-(T_start+T_end)/2)*IL_new;
    dist = dist + V_average * dT;

    % Check if the train is still moving
    if V_new == 0
        if TE_old == TE_new && (vehicle.notchMethod == 0 || notch == vehicle.nNotchT)
            if do_echo
                fprintf(newline)
                cprintf('*err','Warning: ')
                cprintf('err',['Train is stuck at ' num2str(dist/1000) ' km']);
                fprintf(['Acceleration:              ' num2str(Acc_new,5) ' m/s/s' newline])
                fprintf(['TE/BE:                     ' num2str(TE_new,5) ' kN' newline])
                fprintf(['Friction brake:            ' num2str(BE_new,5) ' kN' newline])
                fprintf(['Slope resistance:          ' num2str(F_slope,5) ' kN' newline])
                fprintf(['Train resistance (Davies): ' num2str(force_resistance-force_curve,5) ' kN' newline])
                fprintf(['Curve resistance:          ' num2str(force_curve,5) ' kN' newline])
                fprintf(['Tunnel resistance:         ' num2str(force_tunnel,5) ' kN' newline newline])
                dlg = errordlg(['Train is stuck at ' num2str(dist/1000) ' km'],'Warning');
                uiwait(dlg);
            end
            mdone = 1;
        end
    end
    
    
    if ~isreal(TE_average)
        dbstop;
    end
    if tvehicle.ES_Op_Mode > 0
        limitation.ES_ChargeLevel = limitation.ES_ChargeLevel - ESPow * dT / (tvehicle.ES_capacity * 3600 * tvehicle.nHVSystems);
    end
    
    index = index + 1;
    timedata = allocate(index, timedata);
    timedata(index, 1) = time; % s
    timedata(index, 2) = dist; % m
    timedata(index, 3) = V_new; % m/s
    timedata(index, 4) = TE_new;
    timedata(index, 5) = BE_new;
    timedata(index, 6) = ULine;
    timedata(index, 7) = IL_new;
    timedata(index, 8) = BRPow;
    
    if V_new >= tvehicle.BRFanOffSpeed
        if BRTemp > tvehicle.BRFanHighSpeedTemp
            TC = tvehicle.BRThermTC_FanOnHigh;
        elseif BRTemp > tvehicle.BRFanOnTemp
            TC = tvehicle.BRThermTC_FanOnLow;
        else
            TC = tvehicle.BRThermTC_FanOff;
        end
    else
        if BRTemp > tvehicle.BRFanHighSpeedTemp
            TC = tvehicle.BRThermTC_FanOnLow;
        else
            TC = tvehicle.BRThermTC_FanOff;
        end
    end
    
    tcLow = 2 * TC / (1 + tvehicle.BRThermTC_Tempfactor); % Time constant at ambient temperature
    tcHigh = tcLow * tvehicle.BRThermTC_Tempfactor;
    TC = (BRTemp - tvehicle.AmbientTemp) / (tvehicle.BRTemp_High - tvehicle.AmbientTemp) * (tcHigh - tcLow) + tcLow;
    if ~isnan(TC)
        BRTemp = BRTemp + (tvehicle.AmbientTemp - BRTemp) * (1 - exp(-dT / TC));
    end

    timedata(index, 9) = BRTemp;
    
    timedata(index, 10) = iMotor;         % Motor stator current
    timedata(index, 11) = mSlip;          % Motor slip (Hz)
    timedata(index, 12) = fMode;          % Modulation index
    timedata(index, 13) = ULink;          % DC-link voltage
    timedata(index, 14) = ESPow;          % Energy saver power
    [timedata(index-1:index, :), limitation] = doThermal(timedata(index - 1:index, :), IL_new, V_new, tvehicle, dT, limitation, 0);
    timedata(index, 23) = force_resistance;           % Running resistance
    timedata(index, 24) = PropLoss;       % Traction loss
    timedata(index, 25) = F_slope;        % Slope effective slope force
    timedata(index, 26) = V_average;      % Average speed
    timedata(index, 27) = TE_average;     % Average TE/DBE
    timedata(index, 28) = BE_average;     % Average FricBE
    timedata(index, 29) = IL_average;     % Average line current
    timedata(index, 30) = routedata.track(section_index, 8) * tvehicle.davies_c2(tvehicle.nUnits) * (V_new * 3.6)^2;     % Tunnel resistance
    timedata(index, 31) = notch;
    timedata(index, 32) = force_curve;
    timedata(index, 33) = force_tunnel;
    timedata(index, 34) = iMotor * tvehicle.nMotors;
    if tvehicle.uFreq > 0
        timedata(index, 35) = IS;
    else
        timedata(index, 35) = 0;
    end
    timedata(index, 36:49) = 0;         % allocate
    timedata(index,50) = op_mode;       % 

% 30 Tunnel (effective C-factor of tunnel)
    if dist > endSection % We have reached the end of the section.
        section_index = section_index + 1;
        if section_index < last_section_index % We have not yet reached the next station. From this point we need to go backwards.
            endSection = routedata.track(section_index + 1, 1) * 1000;
            min1 = min([routedata.track(section_index, 2) / 3.6, tvehicle.vMax]) - tRoute.ATO_SpdMgn / 3.6;
            min2 = min([routedata.track(section_index - 1, 2) / 3.6, tvehicle.vMax]) - tRoute.ATO_SpdMgn / 3.6;
            if (V_new > min1) && (min1 < min2) % only if the speed limit in next section is lower we have already exceeded this, we need to go backwards
                index = index - 1;
                % Indicate that we are dealing with the speed reduction
                if do_echo
                    fprintf(['Speed reduction at ', num2str(routedata.track(section_index, 1)), ' km\n'])
                end
                [tempresult, index, time, ttime, tindex, limitation] = backtrack(timedata, section_index, index, routedata, trackdata, slopedata, ...
                    vehicle, tvehicle.class, dT, do_echo, sim_config, min1, hwb, time, limitation, start_clock);
                if ~isempty(tempresult)
                    [timedata, limitation, v_cut, d_cut] = mergeReverseData(timedata(1:index, :), tempresult, v_target, dT, ...
                        tvehicle, limitation, tRoute, sim_config, routedata.track(1,1), routedata.track(end,1), do_echo, hwb);
                    index = index + tindex;
                end
                clear tempresult;
                V_new = min([routedata.track(section_index, 2) / 3.6, tvehicle.vMax]) - tRoute.ATO_SpdMgn / 3.6;
                time = time - ttime - dT;
                
                dist = routedata.track(section_index, 1) * 1000;
                mdone = mdone || v_cut > 0;
                Acc_new = 0;
            end % if (V_new > min1) && (min1 < min2)
            [tvehicle, tRoute, error, op_mode] = check_parameters(routedata, routedata.track(section_index, 9:end), vehicle, dT, do_echo, tvehicle.class);
            mdone = mdone || error;
        else % We have reach next station and need to go backwards until we meet
            if do_echo
                fprintf(['Approach station at ', num2str(routedata.track(section_index, 1)) ' km\n'])
            end
            [tempresult, index, time, ttime, tindex, limitation] = backtrack(timedata, section_index, index - 1, routedata, trackdata, slopedata, ...
                vehicle, tvehicle.class, dT, do_echo, sim_config, 0, hwb, time, limitation, start_clock);
            if ~isempty(tempresult)
                [timedata, limitation, v_cut, d_cut] = mergeReverseData(timedata(1:index, :), tempresult, v_target, dT, ...
                    tvehicle, limitation, tRoute, sim_config, routedata.track(1, 1), routedata.track(end, 1), do_echo, hwb);
                index = index + tindex;
            end
            
            % Station stop
            dist = routedata.track(section_index, 1) * 1000;
            V_new = 0;
            time = time - ttime - dT;
            [TE_average, ULine, ULink, IL_new, BRPow, iMotor, IS, mSlip, fMode, ESPow, PropLoss] = ...
                ShSim_Efforts(0, tvehicle, V_new, tRoute, sim_config, limitation); %#ok<ASGLU>
            
            result.depart = start_time;
            result.arrival = time;
            result.runtime = time - start_time;
            
            BRTemp = timedata(index, 9);
            oldtime = time;
            if trackdata.arrival - result.arrival > 0 % we arrived too early
                n = round((trackdata.dwellTime + trackdata.arrival - result.arrival) / dT);
                time = time + trackdata.dwellTime + trackdata.arrival - result.arrival;
            else
                n = round(trackdata.dwellTime / dT);
                time = time + trackdata.dwellTime;
            end
            timedata(index + 1:index + n, 1)   = ((1:n) * dT + oldtime)';
            timedata(index + 1:index + n, 2)   = routedata.track(section_index, 1) * 1000;
            timedata(index + 1:index + n, 3:5) = 0;
            timedata(index + 1:index + n, 6)   = ULine;
            timedata(index + 1:index + n, 7)   = IL_new;
            timedata(index + 1:index + n, 8)   = 0;   % Brake resistor power
            timedata(index + 1:index + n, 10:12) = 0; % stator current, modulation index
            timedata(index + 1:index + n, 13)  = ULink;
            timedata(index + 1:index + n, 14)  = 0;   % ES power
            timedata(index + 1:index + n, 24)  = PropLoss;
            timedata(index + 1:index + n, 25)  = F_slope;
            timedata(index + 1:index + n, 26)  = 0;      % Average speed
            timedata(index + 1:index + n, 27)  = 0;      % Average TE/DBE
            timedata(index + 1:index + n, 28)  = abs(F_slope) * 1.2;      % Average FricBE
            timedata(index + 1:index + n, 29)  = IL_new;      % Average line current
            timedata(index + 1:index + n, 30)  = routedata.track(section_index, 8);     % Tunnel;
            timedata(index + 1:index + n, 31)  = -1;
            timedata(index + 1:index + n, 32:49) = 0;
            timedata(index + 1:index + n, 50)  = op_mode;
            notch = 0;
            for i1 = index + 1:index + n
                if BRTemp > tvehicle.BRFanHighSpeedTemp
                    TC = tvehicle.BRThermTC_FanOnLow;
                else
                    TC = tvehicle.BRThermTC_FanOff;
                end
                tcLow = 2 * TC / (1 + tvehicle.BRThermTC_Tempfactor); % Time constant at ambient temperature
                tcHigh = tcLow * tvehicle.BRThermTC_Tempfactor;
                TC = (BRTemp - tvehicle.AmbientTemp) / (tvehicle.BRTemp_High - tvehicle.AmbientTemp) * (tcHigh - tcLow) + tcLow;
                BRTemp = BRTemp + (tvehicle.AmbientTemp - BRTemp) * (1 - exp(-dT / TC));
                timedata(i1, 9) = BRTemp; % Brake resistor temperature
                
                [timedata(i1-1:i1, :), limitation] = doThermal(timedata(i1-1:i1, :), IL_new, 0, tvehicle, dT, limitation, 1);
            end
            
            mdone = 1;
            
        end % dist > endSection
        if do_echo
            if ishandle(hwb)
                waitbar((dist / 1000 - routedata.track(1, 1)) / (routedata.track(end, 1) - routedata.track(1, 1)), hwb, [trackdata.from_station, '-', trackdata.to_station]);
            else
                break
            end
        end
    end % if dist > endSection
    time = time + dT;
end %(~mdone) && (time < simulation.maxT)
if time - start_time > sim_config.maxT
    dlg = errordlg(['Maximum running time of ' num2str(sim_config.maxT) ' seconds exceed']);
    uiwait(dlg);
end
if ~exist('result','dir')
    mkdir('result');
end

r = ShSim_Summary(result, timedata, sim_config, trackdata, routedata, tvehicle, data);
timedata(:, 15:20) = sqrt(timedata(:, 15:20));
% for ii=15:20
%     timedata(:,ii)=sqrt(timedata(:,ii));
% end
r.timedata = timedata;
stop_time = timedata(end, 1);

thermal_out.BRTemp        = timedata(end, 9);
thermal_out.iPrimTrafoRMS = timedata(end, 15);    % Main Transformer prim winding RMS
thermal_out.i2ndTrafoRMS  = timedata(end, 16);    % Main Transformer 2nd winding RMS
thermal_out.iLineFiltRMS  = timedata(end, 17);    % Line filter RMS
thermal_out.iPhMCM_RMS    = timedata(end, 18);    % INV phase current RMS
thermal_out.iPhLCM_RMS    = timedata(end, 19);    % CONV phase current RMS
thermal_out.iStatorTM_RMS = timedata(end, 20);    % Traction motor RMS

if do_echo
    fprintf(['Calculation time: ', num2str(round(seconds(datetime - start_clock) * 10) / 10), ' s\n']);
end
return % ShSim_CalcSeg

function [tempresult, index, time, ttime, tindex, limitation] = backtrack(timedata, section_index, index, routedata, ~, slopedata, vehicle, class, dT, do_echo, sim_config, vStart, hwb, time, limitation, start_clock)
% Go backwards. backtrackSection is the section we are reversing on
backtrackSection = section_index - 1;
% Check for parameters
[tvehicle, tRoute, error, op_mode] = check_parameters(routedata, routedata.track(backtrackSection, 9:end), vehicle, dT, do_echo, class); %#ok<ASGLU>
% We start at the speed limitation for the next section = min1
V0_new = vStart;
% Distance from which the speed limit is applied
dist1 = routedata.track(section_index, 1) * 1000;
% Calculate the slope to estimate the effort we have when approaching
F_slope = get_slope(dist1, slopedata) * tvehicle.statmassLoad * 9.81 / 1000;
force_resistance = getResistances(tvehicle, V0_new, dist1, routedata.track(section_index, :), tRoute.Track_Guage, tRoute.headWind);
% Assume zero acceleration when leaving speed restriction
Acc_new = 0;
% Assume last TE (which can be positive)
if vStart > 0
    TE_new = F_slope + force_resistance;
else
    TE_new = -tvehicle.residual_retard * tvehicle.totMass + F_slope + force_resistance;
end
% Initiate time and index
ttime = 0;
tindex = 0;
% Memory pre-allocation
tempresult = zeros(1000, 50);
% Start backwards integration. We run backwards until we reach the speed limit
if do_echo
    if ishandle(hwb)
        waitbar((dist1 / 1000 - routedata.track(1, 1)) / (routedata.track(end, 1) - routedata.track(1, 1)), hwb, 'Back tracking...');
    else
        return;
    end
end
notch = -tvehicle.nNotchB;
while V0_new < min([timedata(index, 3), tvehicle.vMax])
    % Update waitbar
    if toc > 0.2
        if do_echo
            if ishandle(hwb)
                waitbar((dist1 / 1000 - routedata.track(1, 1)) / (routedata.track(end, 1) - routedata.track(1, 1)), hwb);
            else
                tempresult = [];
                return
            end
        end
        % if elapsed time exceeds time-out limit -> escape the loop
        if seconds(datetime - start_clock) > sim_config.time_out
            tempresult = [];
            return
        end
        tic;
    end
    % Step 1 - Calculate first point
    V0_old = V0_new;
    TE_old = TE_new;
    Acc_old = Acc_new;
    dist0 = dist1;
    % Calculate the resistances
    [force_resistance, force_curve, force_tunnel] = getResistances(tvehicle, V0_new, dist0, routedata.track(backtrackSection, :), tRoute.Track_Guage, tRoute.headWind);
    % Get the force induced by the current slope
    F_slope = get_slope(dist0, slopedata) * tvehicle.statmassLoad * 9.81 / 1000;
    % Calculate the acceleration (retardation)
    v_target = max([tRoute.ATO_InitBrkSpd / 3.6, V0_new]) + tRoute.ATO_SpdMgn / 3.6;
    v_target = min([routedata.track(backtrackSection, 2), tvehicle.vMax, v_target, tRoute.ATO_SpdMax / 3.6]) - tRoute.ATO_SpdMgn / 3.6;
    % v_target=min([routedata.track(backtrackSection,2)/3.6 tvehicle.vMax tRoute.ATO_SpdMax/3.6])-tRoute.ATO_SpdMgn/3.6;
    
    [F_target, T_start, T_end] = ShSim_TargetBE(TE_old ,Acc_old, V0_new, v_target, tvehicle, F_slope, force_resistance, dT);
    % First average calc
    TE_average = (T_start + T_end) / 2 * TE_old + (1 - (T_start + T_end) / 2) * F_target;
    % test calculation average acceleration (AveAcc1)
    ret0 = (TE_average - force_resistance - F_slope) / tvehicle.totMass;
    tSpeed = max([V0_old - ret0 * dT, 0]); % Speed1
    % Check max effort in new speed point
    TE_new = ShSim_TargetBE(TE_old, Acc_old, tSpeed, inf, tvehicle, F_slope, force_resistance, dT);
    TE_new = max([F_target, TE_new]);
    % Second average calc
    TE_average = (T_start + T_end) / 2 * TE_old + (1 - (T_start + T_end) / 2) * TE_new;
    Acc_new = (TE_new - force_resistance - F_slope) / tvehicle.totMass;
    A_average = (TE_average - force_resistance - F_slope) / tvehicle.totMass;
    V0_new = max([V0_old - A_average * dT, 0]); % Speed1
    % Check this - not entirely correct, but maybe good enough
    V_average = (V0_old + V0_new) / 2;
    
    % Step down time
    ttime = ttime - dT;
    % Cut the forward record to the current position
    while timedata(index, 2) > dist1
        index = index - 1;
        time = time - dT;
    end % while timedata(index,2)>tdist
    % New distance
    dist1 = dist0 - V_average * dT;
    if dist1 < routedata.track(backtrackSection, 1) * 1000
        backtrackSection = backtrackSection - 1;
        [tvehicle, tRoute, error, op_mode] = check_parameters(routedata, routedata.track(backtrackSection, 9:end), vehicle, dT, do_echo, tvehicle.class); %#ok<ASGLU>
    end
    % Step up index
    tindex = tindex + 1;
    % check if more memory needs to be allocated
    tempresult = allocate(tindex, tempresult);
    % Create record
    tempresult (tindex, 1) = ttime;     % time
    tempresult (tindex, 2) = dist0;     % Distance
    tempresult (tindex, 3) = V0_old;    % Speed
    tempresult (tindex, 4) = TE_old;    % Total effort target - to be adjusted when merged
    tempresult (tindex, 5:22) =nan;	    % Values to be calculated when merged
    tempresult (tindex, 23) = force_resistance;  %Running resistance
    tempresult (tindex, 24) = nan;	    % Propulsion loss
    tempresult (tindex, 25) = F_slope;
    tempresult (tindex, 26) = V_average;
    tempresult (tindex, 27) = TE_average;
    tempresult (tindex, 28) = nan;
    tempresult (tindex, 29) = nan;
    tempresult (tindex, 30) = routedata.track(backtrackSection, 8) * tvehicle.davies_c2(tvehicle.nUnits) * (V0_new*3.6)^2;     %Tunnel resistance;
    tempresult (tindex, 31) = notch;
    tempresult (tindex, 32) = force_curve;
    tempresult (tindex, 33) = force_tunnel;
    tempresult (tindex, 34:49) = 0; % allocate
    tempresult (tindex, 50) = op_mode;
end % vtemp < timedata(index, 3)
tempresult = tempresult(1:tindex, :);
if ~isempty(tempresult)
    tempresult(:, 1) = tempresult(:, 1) - tempresult(end, 1);
    tempresult(:, 1) = tempresult(:, 1) + timedata(index, 1)+dT;
    tempresult = flipud(tempresult(1:tindex, :));
end

return % backtrack

function slope = get_slope(dist, slopedata)
p = find(slopedata(:, 1) <= dist, 1, 'last');
slope = slopedata(p, 2) + (dist - slopedata(p, 1)) * slopedata(p, 3);
return % get_slope

function [F_target, rampStart, rampEnd, mode] = GetTargetTractionForce(TE_old, Acc_old, speed, v_target, tvehicle, F_slope, F_res, dT, limitation)
%
FAccLimit = tvehicle.max_acc * tvehicle.totMass;
switch tvehicle.AccMethod
    case 1
        FAccLimit = FAccLimit + tvehicle.davies_a;
    case 2
        FAccLimit = FAccLimit + F_res;
    case 3
        FAccLimit = FAccLimit + F_res + F_slope;
end % switch
FAccLimit = ShSim_maxTE(FAccLimit, tvehicle, speed, limitation);
if speed >= v_target
    a_target = (v_target - speed) / dT; % target acceleration in new point
    rampStart = 0;
    rampEnd = 0;
    mode = 0;
else % speed<v_target
    t_rampDown = Acc_old / tvehicle.jerkrateT; % the time it takes to ramp down current acceleration
    if isinf(tvehicle.jerkrateT)
        v_rampDown = v_target;
    else
        v_rampDown = v_target - tvehicle.jerkrateT * t_rampDown^2 / 2;    % Speed where ramp down shall start
    end
    t_2rampDown = min(max((v_rampDown - speed) / abs(Acc_old), 0), dT);  % time to where ramp down shall start
    if t_2rampDown < dT % it will happen within this sample
        mode = -1;    % Mode = -1 indicates jerk rate ramp down
        rampStart = t_2rampDown / dT;
        if t_rampDown < (dT - t_2rampDown)
            rampEnd = min((t_2rampDown + t_rampDown) / dT, 1);
            a_target = 0;
        else
            rampEnd = 1;
            a_target = (sqrt((v_target - speed) * 2 / tvehicle.jerkrateT) - dT) * tvehicle.jerkrateT; % Acc_old-(dT-t_2rampDown)*tvehicle.jerkrateT;
        end
    else
        mode = 1;
        rampStart = 0;
        rampEnd = 1;
        a_target = (v_target - speed) / dT;
    end
end
F_target = a_target*tvehicle.totMass + F_res + F_slope;

F_target = min([F_target, FAccLimit]);

maxStep = min([tvehicle.dTE_jerkLimit, tvehicle.maxRatePwrTE / speed * dT]);
if F_target > TE_old + maxStep
    F_target = TE_old + maxStep;
    mode = 1;
    rampStart = 0;
    rampEnd = 1;
elseif F_target<TE_old-tvehicle.dTE_jerkLimit && mode==0
    F_target = TE_old - tvehicle.dTE_jerkLimit;
    mode = -1;
    rampStart = 0;
    rampEnd = 1;
end
F_target = min([tvehicle.maxTE, F_target, tvehicle.totMass * (v_target - speed) / dT + F_res + F_slope]);
if F_target > TE_old % ramping up
    mode = 1;
    rampStart = 0;
    rampEnd = (F_target - TE_old) / maxStep;
end
return % GetTargetTractionForce

function [F_target, notch, notchTime] = Notching(time, oldNotch, lastNotch, TE_old, speed, v_target, Acc_old, tvehicle, mode, F_slope, F_res, dT, limitation)
% lastNotch: (1)=time since the last notch down, (2)=time since the last notch up
notchTime = lastNotch;
if tvehicle.notchMethod == 1 % Time based notching
    if oldNotch >= 0
        F_Old = oldNotch / tvehicle.nNotchT * min([tvehicle.maxTE, tvehicle.maxPwrTEbase / speed]);
    else
        F_Old = oldNotch / tvehicle.nNotchB * min([tvehicle.maxDBE, tvehicle.maxPwrDBEbase / speed]);
    end
    tNotch = min(max(oldNotch + mode, -tvehicle.nNotchB), tvehicle.nNotchT);
    a_target = (sqrt((v_target - speed) * 2 / tvehicle.jerkrateT) - dT) * tvehicle.jerkrateT; % Acc_old-(dT-t_2rampDown)*tvehicle.jerkrateT;
    if Acc_old > a_target
        tNotch = min(0, tNotch);
    end
    if tNotch > oldNotch
        if all(time - lastNotch >= tvehicle.notchrate) || (tNotch <= 0 && time - lastNotch(2) >= tvehicle.notchrate(2))
            notch = tNotch;
            notchTime(2) = time;
        else
            notch = oldNotch;
        end
    elseif tNotch < oldNotch && F_Old >= TE_old
        notch = oldNotch - 1;
        notchTime(1:2) = time;
    else
        notch = oldNotch;
    end
    
    if notch >= 0
        F_target1 = notch * ShSim_maxTE(inf, tvehicle, speed, limitation) / tvehicle.nNotchT;
    else
        F_target1 = notch / tvehicle.nNotchB * min([tvehicle.maxDBE, tvehicle.maxPwrDBEbase / speed]);
    end
    
    if F_target1 >= 0
        if F_target1 > TE_old
            maxStep = min([tvehicle.dTE_jerkLimit, tvehicle.maxRatePwrTE / speed * dT]);
        else
            maxStep = tvehicle.dTE_jerkLimit;
        end
    else
        maxStep = tvehicle.dDBE_jerkLimit;
    end
    F_target = max([min([F_target1, TE_old + maxStep]), TE_old - maxStep]);
elseif tvehicle.notchMethod == 2 % Speed based
    notch = oldNotch;
    F_Max = ShSim_maxTE(inf, tvehicle, speed, limitation);
    F_Min = ShSim_maxBE(tvehicle, speed, F_slope, F_res);
    if oldNotch >= 0
        dSpeed = tvehicle.notchDeltaV(1);
    else
        dSpeed = tvehicle.notchDeltaV(2);
    end
    if speed < v_target - dSpeed
        if all(time > lastNotch + tvehicle.notchrate) && speed < v_target - dSpeed + Acc_old * dT % Full speed
            notch = min([tvehicle.nNotchT, oldNotch+1]);
        elseif oldNotch < tvehicle.nNotchT
            notch = min(ceil((F_slope + F_res) / F_Max * tvehicle.nNotchT), tvehicle.nNotchT);
            notch = max([notch, oldNotch]); % To ensure we don't step down during acceleration
        end
        
    elseif mode == -1 || speed >= v_target
        notch = max(floor((F_slope + F_res) / F_Max * tvehicle.nNotchT), -tvehicle.nNotchB);
    end
    
    a_target = (v_target - speed - [0, dSpeed]) / dT * 1.1;
    F_target = a_target * tvehicle.totMass + F_slope + F_res;
    if notch < oldNotch
        notchTime(1:2) = time;
    elseif notch > oldNotch
        notchTime(2) = time;
    end

    if notch >= 0
        F_target1 = min([notch / tvehicle.nNotchT * F_Max, F_target(1)]);
    else
        F_target1 = max([-notch / tvehicle.nNotchB * F_Min F_target(2)]);
    end

    if F_target1 >= 0
        if F_target1 > TE_old
            maxStep = min([tvehicle.dTE_jerkLimit, tvehicle.maxRatePwrTE / speed * dT]);
        else
            maxStep = tvehicle.dTE_jerkLimit;
        end
    else
        maxStep = tvehicle.dDBE_jerkLimit;
    end
    F_target = max([min([F_target1, TE_old + maxStep]), TE_old - maxStep]);
elseif tvehicle.notchMethod == 3 % controller
    a_target = (v_target - speed) / tvehicle.notchDeltaVAcc * 3.6 * tvehicle.notchMaxAcc;
    F_target = a_target * tvehicle.totMass + F_slope + F_res;
    F_Max = ShSim_maxTE(inf, tvehicle, speed, limitation);
    F_Min = ShSim_maxBE(tvehicle, speed, F_slope, F_res);
    if F_target >= 0
        tnotch=floor(F_target / F_Max * tvehicle.nNotchT);
    else
        tnotch = -ceil(F_target / F_Min * tvehicle.nNotchB);
    end
    notch = oldNotch;
    if time > lastNotch(2) + tvehicle.notchrate(2)
        if tnotch > oldNotch && oldNotch < tvehicle.nNotchT
            notch = oldNotch + 1;
            notchTime(2) = time;
        elseif tnotch < oldNotch && oldNotch > -tvehicle.nNotchB
            notch = oldNotch - 1;
            notchTime(2) = time;
        end
    end
    if notch >= 0
        F_target1 = notch / tvehicle.nNotchT * F_Max;
    else
        F_target1 = -notch / tvehicle.nNotchB * F_Min;
    end
    if F_target1 >= 0
        if F_target1 > TE_old
            maxStep = min([tvehicle.dTE_jerkLimit, tvehicle.maxRatePwrTE / speed * dT]);
        else
            maxStep = tvehicle.dTE_jerkLimit;
        end
    else
        maxStep = tvehicle.dDBE_jerkLimit;
    end
    F_target = max([min([F_target1, TE_old + maxStep]), TE_old - maxStep]);
end
return %Notching

function [tvehicle, troutedata, error, op_mode]=check_parameters(routedata, values, vehicle, delta_T, do_echo, class)
tvehicle = vehicle;
troutedata = routedata;
error = 0;
op_mode = 0;
if ~isempty(troutedata.parameters)
    for ii = 1:size(troutedata.parameters, 2)
        [pName, idx] = stripPar(troutedata.parameters(ii).name);
        if isfield(tvehicle, pName)
            if size(values, 2) >= ii
                if ~isnan(values(ii))
                    if idx > 0
                        tvehicle.(char(pName))(idx) = values(ii);
                    else
                        tvehicle.(char(pName)) = values(ii);
                    end
                end
            end
        end
        tvehicle = ShSim_CalcMassRes(tvehicle);
        if isfield(troutedata, troutedata.parameters(ii).name)
            if size(values, 2) >= ii
                if ~isnan(values(ii))
                    troutedata.(char(troutedata.parameters(ii).name)) = values(ii);
                end
            end
        end
        if strcmpi(char(troutedata.parameters(ii).name), 'OpMode')
            op_mode = values(ii);
        end
        if strcmpi(char(troutedata.parameters(ii).name), 'vehicle')
            if values(ii) > 0
                if exist(fullfile(vehicle.path, char(troutedata.parameters(ii).vehicles(values(ii)))), 'file')
                    tvehicle = ShSim_Read_Vehicle(vehicle.path, char(troutedata.parameters(ii).vehicles(values(ii))), 0);
                    if do_echo && ~strcmp(class, tvehicle.class)
                        fprintf(['Changing to vehicle ''',  char(troutedata.parameters(ii).vehicles(values(ii))), '''\n']);
                        if abs((tvehicle.unitLength - vehicle.unitLength) / vehicle.unitLength) > 0.05
                            warning('Vehicle length differ. Speed limits will not be valid.');
                        end
                    end
                else
                    dlg = errordlg(['Failed to load vehicle file:' newline char(troutedata.parameters(ii).vehicles(values(ii))) newline 'Check file name and path'], 'File not found');
                    uiwait(dlg);
                    tvehicle = vehicle;
                    error = 1;
                end
            end
        end
    end
end
if any(tvehicle.LoadCondition == 1:5) 
    tvehicle.dTE_jerkLimit = tvehicle.jerkrateT * tvehicle.totMass * delta_T;
    if tvehicle.dTE_jerkLimit == 0
        tvehicle.dTE_jerkLimit = inf;
        tvehicle.jerkrateT = inf;
    end
    tvehicle.dDBE_jerkLimit = tvehicle.jerkrateB * tvehicle.totMass * delta_T;
    if tvehicle.dDBE_jerkLimit == 0
        tvehicle.dDBE_jerkLimit = inf;
        tvehicle.jerkrateT = inf;
    end
    tvehicle = ShSim_CalcMassRes(tvehicle);
else
    dlg = errordlg('Parameter ''LoadCondition'' must be an integer between 1..5', 'Illegal value');
    uiwait(dlg);
    error = 1;
end
return % check_parameters

function [pName, idx] = stripPar(inName)
p = strfind(char(inName), '(');
if isempty(p)
    idx = 0;
    pName = inName;
else
    inName = char(inName);
    idx = num2str(inName(p + 1:end - 1));
    pName = cellstr(inName(1:p - 1));
end
return % stripPar

function [r, limitation_out, v_cut, d_cut] = mergeReverseData(timedata, reversedata, ~, dt, vehicle, limitation, routedata, simulation, fromdist, todist, do_echo, hwb)
lnRD = size(timedata, 1);
r = [timedata; reversedata];
BRTemp = timedata(end, 9); % Incoming brake resistor temperature
limitation.thermPower = 1; 
limitation.thermEffort = 1;
red = (vehicle.BRTemp_Max - BRTemp) / (vehicle.BRTemp_Max - vehicle.BRTemp_High);
if ~isnan(red)
    limitation.thermPower = min([red, limitation.thermPower]);
end
speed = r(lnRD, 3);
TE_ref = r(lnRD, 4);

if lnRD > 1 && size(reversedata, 1) > 1
    a_end = min([(r(lnRD, 3) - r(lnRD - 1, 3)) / dt, TE_ref / vehicle.totMass]);
    a_brk = min([(reversedata(1, 3) - reversedata(2, 3)) / dt, -reversedata(1, 4) / vehicle.totMass]);
    if (a_end > vehicle.jerkrateT * dt || a_brk > vehicle.jerkrateB * dt) && (a_end + a_brk > (vehicle.jerkrateT + vehicle.jerkrateB) * dt) && (a_end > 0)
        x = lnRD;
        b = lnRD;
        d_cut = [0, inf];
        while r(b, 2) < d_cut(2)
            a1 = (r(x, 3) - r(x - 1, 3)) / dt;
            jr = min([vehicle.maxRatePwrTE / vehicle.totMass / timedata(end, 3), vehicle.jerkrateT]);
            t1 = a1 / jr * 1.1;
            v_cut = r(x, 3) + a1 * t1 - jr * t1^2 / 2;
            dd1 = r(x, 3) * t1 + a1 * t1^2 / 2 - jr * t1^3 / 8;
            a2 = -inf;
            while a2 < a_brk
                a2 = a_brk;
                t2 = a2 / vehicle.jerkrateB * 1.1;
                v_end = v_cut - a2 * t2 + vehicle.jerkrateB * t2^2 / 2;
                dd2 = v_cut * t2 - a2 * t2^2 / 2 + vehicle.jerkrateB * t2^3 / 8;
                d_cut = r(x, 2) + [-dd1, dd1 + dd2];
                b = find(r(lnRD:end, 3) < v_end, 1, 'first') + lnRD - 1;
                if isempty(b) || b == size(r, 1)
                    break
                else
                    a_brk = (r(b, 3) - r(b + 1, 3)) / dt;
                end
            end
            if r(b, 2) < d_cut(2)
                x = x - 1;
            end
        end
        if abs(1 - v_cut / r(lnRD, 3)) < 0.001
            v_cut = 0;
            d_cut = 0;
        end
    else
        v_cut = 0;
        d_cut = 0;
    end
else
    v_cut = 0;
    d_cut = 0;
end

[DBE_effort, ULine, ULink, ILine, BRPow, iMotor, IS, mSlip, fMode, ESPow, PropLoss]=...
    ShSim_Efforts(TE_ref, vehicle, speed, routedata, simulation, limitation); %#ok<ASGLU>

for ii = lnRD + 1:size(r, 1)
    if toc > 0.1
        if do_echo
            if ishandle(hwb)
                waitbar((r(ii, 2) / 1000 - fromdist) / (todist - fromdist), hwb);
            end
        end
        tic;
    end
    speed = r(ii, 3);

    TE_ref = r(ii, 4);
    Brk_average = r(ii, 27);
    V_average = r(ii, 26);
    [DBE_effort, ULine, ULink, ILine, BRPow, iMotor, IS, mSlip, fMode, ESPow, PropLoss]=...
        ShSim_Efforts(TE_ref, vehicle, speed, routedata, simulation, limitation);  %#ok<ASGLU>
    [DBE_Average, ULine, ULink, IL_Average, BRPow, iMotor, IS, mSlip, fMode, ESPow, PropLoss]=...
        ShSim_Efforts(Brk_average, vehicle, V_average, routedata, simulation, limitation); 
    
    r(ii, 4) = DBE_effort;
    r(ii, 5) = DBE_effort-TE_ref;
    r(ii, 6) = ULine;
    r(ii, 7) = ILine;
    r(ii, 8) = BRPow;
    WBr = BRPow * dt;
    if speed >= vehicle.BRFanOffSpeed
        if BRTemp > vehicle.BRFanHighSpeedTemp
            TC = vehicle.BRThermTC_FanOnHigh;
        elseif BRTemp > vehicle.BRFanOnTemp
            TC = vehicle.BRThermTC_FanOnLow;
        else
            TC = vehicle.BRThermTC_FanOff;
        end
    else
        if BRTemp > vehicle.BRFanHighSpeedTemp
            TC = vehicle.BRThermTC_FanOnLow;
        else
            TC = vehicle.BRThermTC_FanOff;
        end
    end
    tcLow = 2 * TC / (1 + vehicle.BRThermTC_Tempfactor); % Time constant at ambient temperature
    tcHigh = tcLow * vehicle.BRThermTC_Tempfactor;
    TC = (BRTemp - vehicle.AmbientTemp) / (vehicle.BRTemp_High - vehicle.AmbientTemp) * (tcHigh - tcLow) + tcLow;
    BRThermC = vehicle.BRThermC * vehicle.nINV * vehicle.nHVSystems;
    BRTemp = BRTemp + WBr / BRThermC + (vehicle.AmbientTemp - BRTemp) * (1 - exp(-dt / TC)); % Brake resistor temperature
    red = (vehicle.BRTemp_Max - BRTemp) / (vehicle.BRTemp_Max - vehicle.BRTemp_High);
    limitation.thermPower = min([red, limitation.thermPower]);
    
    r(ii, 9)  = BRTemp;
    r(ii, 10) = iMotor;         % Motor stator currentp
    r(ii, 11) = mSlip;          % Motor slip (Hz)
    r(ii, 12) = fMode;          % Modulation index
    r(ii, 13) = ULink;          % DC-link voltage
    r(ii, 14) = ESPow;          % Energy saver power

    [r(ii-1:ii, :), limitation] = doThermal(r(ii-1:ii, :), ILine, speed, vehicle, dt, limitation, 1);
    
    r(ii, 24) = (PropLoss + r(ii - 1, 24)) / 2;
    r(ii, 27) = DBE_Average;   % 27 Average TE/DBE
    r(ii, 28) = DBE_Average-Brk_average;   %28 Average FricBE
    r(ii, 29) = IL_Average;
    r(ii, 34) = iMotor * vehicle.nMotors;
    if vehicle.uFreq > 0
        r(ii, 35) = IS;
    else
        r(ii, 35) = 0;
    end
end
limitation_out = limitation;
return % mergeReverseData

function [r, limitation] = doThermal(r_in, ILine, speed, vehicle, dt, limitation_in, inBrk)
r = r_in;
limitation = limitation_in;
if vehicle.uFreq > 0
% Primary line current
    r(2, 15) = r(1, 15) + (ILine^2 - r(1, 15)) * (1 - exp(-dt / vehicle.IPrimTrafoContTC));
    if r(2, 15) > vehicle.IPrimTrafoCont^2
        if abs(ILine) > vehicle.IPrimTrafoCont * (1 - dt / 10)
            limitation.IL_red = limitation.IL_red * (1 - dt / 10);
        end
    else
        limitation.IL_red = min([limitation.IL_red * (1 + dt / 10), 1]);
    end
    limitation.thermPower = min([limitation.IL_red, 1]);
% Main Transformer 2nd winding RMS
    ISec = ILine * vehicle.MT_ratio / (vehicle.nHVSystems * vehicle.nCONV * 2);
    r(2, 16) = r(2 - 1, 16) + (ISec^2 - r(1, 16)) * (1 - exp(-dt / vehicle.I2ndTrafoContTC));
    if r(2, 16) > vehicle.I2ndTrafoCont^2
        if abs(ISec) > vehicle.I2ndTrafoCont * (1 - dt / 10)
            limitation.I2nd_red = limitation.I2nd_red * (1 - dt / 10);
        end
    else
        limitation.I2nd_red = min([limitation.IL_red * (1 + dt / 10), 1]);
    end
    limitation.thermPower = min([limitation.thermPower, limitation.I2nd_red]);
% Line filter RMS (n/a in AC mode)
    r(2, 17) = 0;
% CONV phase current RMS
    r(2, 19) = r(2 - 1, 19) + (ISec^2 - r(1, 19)) * (1 - exp(-dt / vehicle.IphLCMContTC));
    if r(2, 19) > vehicle.IphLCMCont^2
        if abs(ISec) > vehicle.IphLCMCont * (1 - dt / 10)
            limitation.ILC_red = limitation.ILC_red * (1 - dt / 10);
        end
    else
        limitation.ILC_red = min([limitation.ILC_red * (1 + dt / 10), 1]);
    end
    limitation.thermPower = min([limitation.thermPower, limitation.ILC_red]);
else % DC mode
% Main Transformer prim winding RMS
    r(2, 15) = 0;
% Main Transformer 2nd winding RMS
    r(2, 16) = 0;
% Line filter RMS
    if isinf(vehicle.IDCFilterCont)
        r(2, 17) = inf;
    else
        r(2, 17) = r(1, 17) + (ILine^2 - r(1, 17)) * (1 - exp(-dt / vehicle.IDCFilterContTC));
        red = max([(vehicle.IDCFilterCont^2 / r(2, 17) - 1) / 0.04, 0]);
        limitation.thermPower = min([limitation.thermPower, red]);
    end
% CONV phase current RMS
    r(2, 19) = 0;
end
% INV phase current RMS
%     MotorRPM=speed*tvehicle.gearRatio/(tvehicle.wheelsize*pi)*60;
iMotor = r_in(2, 10);
if isinf(vehicle.IphMCMContHex)
    r(2, 18) = inf;
else
    if speed < vehicle.Speed0(1 + inBrk) % Sinus modulation
        r(2, 18) = r(1, 18) + ((iMotor * vehicle.nMotors * vehicle.IphMCMContHex / vehicle.IphMCMContSin)^2 ...
            - r(1, 18)) * (1 - exp(-dt / vehicle.IphMCMContTC));
    else %Hex
        r(2, 18) = r(1, 18) + ((iMotor * vehicle.nMotors)^2 - r(1, 18)) * (1 - exp(-dt / vehicle.IphMCMContTC));
    end
    if ~isinf(vehicle.IphMCMContHex)
        red = max([(vehicle.IphMCMContHex^2 / r(2, 18) - 1) / 0.04, 0]);
    else
        red = 1;
    end
    limitation.thermEffort = min([limitation.thermEffort, red]);
end

if vehicle.TMventType == 1
    mFactor = sqrt((vehicle.Speed0(1 + inBrk) + 1) / (r(2, 3) + 1));
else
    mFactor = 1;
end

% Traction motor RMS
if isinf(vehicle.IstatTMCont)
    r(2, 20) = inf;
else
    r(2, 20) = r(1, 20) + ((iMotor)^2 * mFactor - r(1, 20)) * (1 - exp(-dt / mFactor / vehicle.IstatTMContTC));
    red = (vehicle.IstatTMCont^2 / r(2, 20) - 1) / 0.04;
    limitation.thermEffort = min([limitation.thermEffort, red]);
end
r(2, 21) = limitation.thermPower;      %Power Reduction factor
r(2, 22) = limitation.thermEffort;     %Effort Reduction factor
return % doThermal

function [Coast, F_target] = ATO_control(Coasting, speed, F_target_in, troutedata, dist, tvehicle, F_slope, do_echo, dT)
if Coasting
    F_target = min(0, F_target_in);
    if speed <= troutedata.ATO_vCoastEnd / 3.6
        Coast = 0;
    elseif dist > troutedata.ATO_dCoastEnd * 1000
        Coast = 0;
    else
        Coast = 1;
    end
else
    Coast = 0;
    if speed >= troutedata.ATO_vCoastStart / 3.6
        Coast = 1;
        if do_echo
            fprintf(['Coasting at: ', num2str(dist / 1000, '%4.1f'), ' km\n'])
        end
    elseif dist < troutedata.ATO_dCoastStart * 1000 + speed * dT * 2
        if dist > troutedata.ATO_dCoastStart * 1000
            Coast = 1;
            if do_echo
                fprintf(['Coasting at: ', num2str(dist / 1000, '%4.1f'), ' km\n'])
            end
        end
    end
    
    if speed > troutedata.ATO_V1 / 3.6
        ATO_force = max([0, (1 - (speed * 3.6 - troutedata.ATO_V1) / (troutedata.ATO_V2 - troutedata.ATO_V1)) * tvehicle.maxPwrTEbase / speed + F_slope * troutedata.ATO_cGrad]);
        F_target = min([F_target_in, ATO_force]); % + F_slope*troutedata.ATO_cGrad
    else
        F_target = F_target_in;
    end
    F_target = F_target * troutedata.ATO_TBC;
end
return % ATO_control

function [force_resistance, curve_resistance, tunnel_resistance] = getResistances(tvehicle, speed, dist, track, gauge, wind)
% getResistances
% input
%   tvehicle    Vehicle structure
%   speed       Train speed in m/s
%   dist        distance in meter on the line
%   track       1x8 vector with 
%                   1 - dist
%                   2 - xx
%                   3 - xx
%                   4 - base curvature rad/m
%                   5 - curvature gradient rad/m/m
%                   6 - base cant mm
%                   7 - cant gradient mm/m
%                   8 - Tunnel
% output
%   force_resistance    Total train resistance in this point
%   curve_resistance    curve resistance part of above
%   tunnel_resistance   additional resistance force caused by the tunnel
%
curvature = abs((track(4) + (dist / 1000 - track(1)) * track(5)) / 1000); % rad/m
if curvature == 0
    curve_resistance = 0;
else
    cant = track(6) + (dist / 1000 - track(1)) * track(7); % mm
%     speed0=sqrt(cant*9.81*curvature/(gauge*1000))*1000/curvature;
    curve_resistance = (tvehicle.curve_cr2 * tvehicle.statmassLoad + tvehicle.haul_curve_cr2 * tvehicle.haul_mass) * (9.81 * cant / gauge - curvature * speed^2);
    curve_resistance = curve_resistance + tvehicle.curve_cr0 * tvehicle.statmassLoad * curvature / (1 - tvehicle.curve_cr1 * curvature);
    curve_resistance = curve_resistance + tvehicle.haul_curve_cr0 * tvehicle.haul_mass * curvature / (1 - tvehicle.haul_curve_cr1 * curvature);
end

if track(8) %tunnel
    force_davies = (tvehicle.davies_a + tvehicle.davies_b * speed * 3.6 + track(8) * tvehicle.davies_c2(tvehicle.nUnits) * (speed * 3.6)^2);
    tunnel_resistance = track(8) * tvehicle.davies_c2 * (speed*3.6)^2 - tvehicle.davies_c1 * (speed*3.6)^2;
else %not tunnel
    force_davies = (tvehicle.davies_a + tvehicle.davies_b * (wind + speed) * 3.6 + tvehicle.davies_c1(tvehicle.nUnits) * ((wind + speed) * 3.6)^2);
    tunnel_resistance = 0;
end
force_resistance = force_davies + curve_resistance;
return % getResistances

function r = allocate(idx, r)
if idx > size(r, 1)
    r(idx + 1000, 1) = 0;
end % allocate