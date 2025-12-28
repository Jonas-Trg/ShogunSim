function locResult = ShSim_Calc(data, sim_config)
% ShSim_Calc
% 0.1.1.2
% ShSim_Calc
% The Shogun Tracion Performance Calculation core
% Autor: Jonas Rosengren
%
% Format:
%   r = ShSim_Calc(data)
%
% Description
%   This function is managing the calculation of the route.
%
% Input data
%   struct data (input data)
%       vehicle_file
%       track_file
%       doOptimise
%       parameters (each field name, value)
%   sim_config
%
% Output data
%   struct r (output data)
%            fName: 'Station A-Station B_Cl375 1x4DC750A60R25_047.mat'
%     from_station: 'Station A'
%       to_station: 'Station B'
%          vehicle: 'Cl375 1x4DC750A60R25'
%          runtime: 351.3000
%           depart: 0
%          arrival: 351.3000
%       t_powerMax: 0
%       t_powerMin: 0
%      t_effortMax: 0
%      t_effortMin: 0
%         t_cruise: 0
%          t_coast: 0
%       w_traction: 6.4725e+004
%          w_regen: 7.4400e+003
%       w_friction: 3.2596e+004
%      i_rms_sqsum: 4.1929e+009
%        i_rms_num: 4115
%     i_rms2_sqsum: 4.1929e+009
%       i_rms2_num: 3513
%       parameters: same as input
%          timedata: run details
%
% Calculation parameters
%   doOptimise
%   Parameter list
%
% Output data
%   Result of calcuation
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-19  JR      Cleaning code.
%
warning off MATLAB:singularMatrix
locResult = [];
if ~isfield(data, 'echo')
    data.echo = 1;
end
if ~isfield(data, 'skipTimeTable')
    data.skipTimeTable = 0;
end
if ~isfield(sim_config, 'ThermReduction')
    sim_config.ThermReduction = 0;
end
if ~isfield(sim_config, 'useSubStations')
    sim_config.useSubStations = 1;
end
if ~isfield(sim_config, 'plotRef')
    sim_config.plotRef = 2;
end
global logfileref %#ok<GVMIS> 

if exist(data.vehicle_file, 'file')
    locVehicle = ShSim_Read_Vehicle(data.vehicle_file);
    if isempty(locVehicle)
        return
    end
else
    return
end
% Replace vehicle parameters
if isfield(data, 'parameter')
    for ii = 1:size(data.parameter, 2)
        if isfield(locVehicle, data.parameter(ii).name)
            locVehicle.(data.parameter(ii).name) = data.parameter(ii).value;
        end
    end
    locVehicle = ShSim_CalcMassRes(locVehicle);
end

if exist(data.track_file, 'file')
    [trackdata, routedata, routeProfile] = ShSim_Read_trackdata(data.track_file, locVehicle, sim_config);
    if isempty(trackdata)
        return
    end
else
    trackdata = [];
end
if isempty(trackdata)
    return
end
fileattrib(data.vehicle_file, '-w'); % Set file read only
fileattrib(data.track_file, '-w'); % Set file read only

% Replace route parameters
if isfield(data, 'parameter')
    for ii = 1:size(data.parameter, 2)
        if isfield(routedata, data.parameter(ii).name)
            routedata.(data.parameter(ii).name) = data.parameter(ii).value;
        end
    end
end
gradientdata = ShSim_CalcGradients(locVehicle, routeProfile, data.echo);
plot_fig = data.echo;

number_of_sections = size(trackdata, 2);
if number_of_sections < 1
    return
end

start_time_stamp = datetime;
thermal_in.BRTemp          = locVehicle.AmbientTemp;
stoptime = trackdata(end).dwellTime;
thermal_in.iPrimTrafoRMS   = locVehicle.IPrimTrafoCont * exp(-stoptime / locVehicle.IPrimTrafoContTC) * 0.98;
thermal_in.i2ndTrafoRMS    = locVehicle.I2ndTrafoCont  * exp(-stoptime / locVehicle.I2ndTrafoContTC)  * 0.98;
thermal_in.iLineFiltRMS    = locVehicle.IDCFilterCont  * exp(-stoptime / locVehicle.IDCFilterContTC)  * 0.98;
thermal_in.iPhMCM_RMS      = locVehicle.IphMCMContHex  * exp(-stoptime / locVehicle.IphMCMContTC)     * 0.98;
thermal_in.iPhLCM_RMS      = locVehicle.IphLCMCont     * exp(-stoptime / locVehicle.IphLCMContTC)     * 0.98;
thermal_in.iStatorTM_RMS   = locVehicle.IstatTMCont    * exp(-stoptime / locVehicle.IstatTMContTC)    * 0.98;
if ~isfield(data, 'sections')
    data.sections = [];
end
if isempty(data.sections)
    data.sections = 1:number_of_sections;
end
if data.skipTimeTable
    time_dep = 0;
else
    time_dep = trackdata(data.sections(1)).departure;
end

if ~isfield(data, 'doOptimise')
    data.doOptimise = 0;
end
if nargin > 2
    p = data.sections <= number_of_sections;
    if ~any(data.sections == 0)
        data.sections=data.sections(p);
    end
end
if ~isempty(logfileref)
    fprintf(logfileref, [char(datetime, 'HH:mm:ss') ': Route calculated with:\r\n']);
    fprintf(logfileref, [pad('', 10), pad('Track:', 13), strrep(data.track_file, '\', '\\'), newline]);
    fprintf(logfileref, [pad('', 10), pad('Vehicle:', 13), strrep(data.vehicle_file, '\', '\\'), newline]);
    if size(data.sections, 2) > 1
        fprintf(logfileref, [pad('', 10), pad('Sections:', 13), num2str(data.sections(1)), '-', num2str(data.sections(end)), newline]);
    else
        fprintf(logfileref, [pad('', 10), pad('Sections:', 13), num2str(data.sections), newline]);
    end
    
    if size(data.doOptimise, 2) == 1
        if data.doOptimise > 0
            fprintf(logfileref, [pad('', 10), pad('Optimise:', 13), 'method: ', num2str(data.doOptimise), newline]);
        else
            fprintf(logfileref, [pad('', 10), pad('Optimise:', 13), 'None', newline]);
        end
    elseif size(data.doOptimise, 2) == 2
        fprintf(logfileref, [pad('', 10), pad('Optimise:', 13), 'method: ', num2str(data.doOptimise(1)), '\tSpeed relation: ', num2str(data.doOptimise(2)), newline]);
    elseif size(data.doOptimise, 2) == 3
        fprintf(logfileref, [pad('', 10), pad('Optimise:', 13), 'method:', num2str(data.doOptimise(1)), '\tV2/V1: ', num2str(data.doOptimise(2)), '\tV1/vBrk: ', num2str(data.doOptimise(3)), newline]);
    else
        fprintf(logfileref,[pad('', 10), pad('Optimise:', 13), 'None', newline]);
    end
    fprintf(logfileref, [pad('', 10), pad('Delta-T:', 13), num2str(sim_config.delta_T * 1000), newline]);
end
clear locResult;

index = 0;
delay = 0;
lastFig = 0;
if data.echo
    hwb = waitbar(0, '');
    ps = hwb.Position;
    ps(2) = ps(2) + 80;
    s = ['Calculating route: ', newline, trackdata(data.sections(1)).from_station, ' - ', trackdata(data.sections(end)).to_station];
    hwb0 = waitbar(0, s, 'Position', ps);
    hwb0.Position = ps;
else
    hwb = nan;
    hwb0 = nan;
end
for ii = data.sections
    prog = (ii - data.sections(1)) / (data.sections(end) - data.sections(1));
    if ishandle(hwb0)
        waitbar(prog, hwb0);
    end
    vehicle = locVehicle;
    index = index + 1;
    v_cut = inf;
    if data.skipTimeTable
        trackdata(ii).arrival = 0;
    end
    while v_cut > 0
        if trackdata(ii).from_track < trackdata(ii).to_track
            [r_org, time_arr, thermal_out, v_cut, d_cut]  = ...
                ShSim_CalcSeg(trackdata(ii), vehicle, sim_config, time_dep, thermal_in, routedata, gradientdata, data, hwb);
        else
            revTrack = trackdata(ii);
            tRoutedata = routedata;
            tRoutedata.track = routedata.revtrack;
            tRoutedata.revtrack = routedata.track;
            revTrack.from_track = size(routedata.track, 1) - revTrack.from_track;
            revTrack.to_track = size(routedata.track, 1) - revTrack.to_track;
            [r_org, time_arr, thermal_out, v_cut, d_cut]  = ...
                ShSim_CalcSeg(revTrack, vehicle, sim_config, time_dep, thermal_in, tRoutedata, gradientdata, data, hwb);
        end
        if v_cut > 0
            v_cut = v_cut + routedata.ATO_SpdMgn / 3.6;
            track = trackdata(ii).track;
            b = find(track(:, 1) < d_cut(1) / 1000, 1, 'last');
            if ~isempty(b)
                r1 = track(b, :); 
                r1(1, 1) = d_cut(1) / 1000;
                r1(1, 2) = v_cut * 3.6;
                track = [track(1:b, :); r1; track(b + 1:end, :)];
            else
                track(1, 2) = v_cut * 3.6;
                b = 0;
            end
            d_end = (d_cut(2) + vehicle.unitLength * vehicle.nUnits) / 1000;
            c = find(track(:, 1) < d_end, 1, 'last');
            r2 = track(c, :);
            if c < size(track, 1) && track(c, 2) > v_cut * 3.6
                r2(1, 1) = d_end;
                r2(1, 2) = v_cut * 3.6;
                r2(1, 4) = track(c, 4) + track(c, 5) * (r2(1, 1) - track(c, 1));
                r2(1, 6) = track(c, 6) + track(c, 5) * (r2(1, 1) - track(c, 1));
                track = [track(1:c, :); r2; track(c + 1:end, :)];
            end
            fprintf(['Speed: ', num2str(v_cut * 3.6), ' km/h between ', num2str(d_cut(1) / 1000), ' - ', num2str(d_end), newline])
            track(b + 1:c, 2) = min(track(b + 1:c, 2), v_cut * 3.6);
            trackdata(ii).track = track;
        end
    end
    if data.doOptimise(1) % 0 = no optimisation, 1 = ATO_SpdMax, 2 = ATO_InitBrkSpd, 3 = ATO_TBC
        switch data.doOptimise(1)
            case 1
                par{1} = 'ATO_SpdMax';
                orgval = GetValues(par, routedata);
                if routedata.(par{1}) == inf
                    routedata.(par{1}) = max(r_org.timedata(:, 3)) * 3.6; % vehicle.vMax*3.6;
                end
            case 2
                par{1} = 'ATO_InitBrkSpd';
                orgval = GetValues(par, routedata);
                if routedata.(par{1}) == inf
                    routedata.(par{1}) = max(r_org.timedata(:, 3)) * 3.6; % vehicle.vMax*3.6;
                end
            case 3
                par{1} = 'ATO_TBC';
                orgval = GetValues(par, routedata);
            case 4 % ATO_V1 & ATO_V2
                par{1} = 'ATO_V1';
                par{2} = 'ATO_V2';
                orgval = GetValues(par, routedata);
                routedata.ATO_V1 = max(r_org.timedata(:, 3)) * 3.6;
                routedata.ATO_V2 = data.doOptimise(2) * routedata.ATO_V1;
            case 5 % ATO_SpdMax & ATO_InitBrkSpd
                par{1} = 'ATO_SpdMax';
                par{2} = 'ATO_InitBrkSpd';
                orgval = GetValues(par, routedata);
                routedata.ATO_InitBrkSpd = max(r_org.timedata(:, 3)) * 3.6;
                routedata.ATO_SpdMax = routedata.ATO_InitBrkSpd * data.doOptimise(2);
            case 6 % ATO_V1 & ATO_V2 & ATO_InitBrkSpd
                par{1} = 'ATO_V1';
                par{2} = 'ATO_V2';
                par{3} = 'ATO_InitBrkSpd';
                orgval = GetValues(par, routedata);
                routedata.ATO_InitBrkSpd = max(r_org.timedata(:, 3)) * 3.6;
                routedata.ATO_V1 = routedata.ATO_InitBrkSpd * data.doOptimise(3);
                routedata.ATO_V2 = routedata.ATO_V1 * data.doOptimise(2);
            case 7 % User defined
        end
        if plot_fig
            if ishandle(199 + ii)
                close(199 + ii)
            end
            if lastFig && ishandle(lastFig)
                close(lastFig)
            end
            stdfig([6, 16], 199 + ii, ['ECO Optimisation for ', r_org.from_station, ' - ', r_org.to_station]); 
            hold on; grid on;
            lastFig = 199 + ii;
            k_plot = plot(r_org.timedata(:, 2) / 1000, r_org.timedata(:, 3) * 3.6, 'k', 'LineWidth', 3);
            title(['ECO Optimisation for ', r_org.from_station, ' - ', r_org.to_station]);
            xlabel('Distance [km]');
            ylabel('Speed [km/h]');
            drawnow;
        end
        r_temp = r_org;
        target_runtime = trackdata(ii).arrival - trackdata(ii).departure;
        val = GetValues(par, routedata);
        tdata = data;
        tdata.echo = 0;
        if target_runtime > r_temp.runtime+delay
            if data.echo
                fprintf([' => Target: ', num2str(target_runtime - delay), ' s\n'])
                fprintf([' => Initial time: ', num2str(r_temp.runtime + delay), ' s\n'])
            end
            delta = 0.5;
            
            oldmin = 0;
            min_found = 0;
            oldmax = inf;
            max_found = 0;
            b_plot = [];
            r_plot = [];
            while abs(target_runtime - (r_temp.runtime + delay)) >= sim_config.delta_T % exact time not found
                if target_runtime > r_temp.runtime + delay % Too fast
                    if min_found && max_found
                        break;
                    end
                    routedata = OffsetValues(par, -val * delta, routedata);

                    % extend running time
                    if plot_fig && ishandle(199 + ii)
                        set(groot, 'CurrentFigure', 199 + ii);
                        b_plot = plot(r_temp.timedata(:, 2) / 1000, r_temp.timedata(:, 3) * 3.6, 'b');
                        drawnow;
                    end
                    echoPars(par, routedata, data.echo);
                    [r_temp, time_arr, thermal_out]  = ShSim_CalcSeg(trackdata(ii), vehicle, sim_config, time_dep, thermal_in, routedata, gradientdata, tdata, 0);
                    if data.echo
                        fprintf([' => Time: ', num2str(r_temp.runtime), ' (diff=', num2str(target_runtime - (r_temp.runtime + delay)), ' s)\n'])
                    end
                    if oldmin == r_temp.runtime && r_temp.runtime > 0
                        min_found = 1; %#ok<NASGU>
                        break;
                    elseif r_temp.runtime > 0
                        oldmin = r_temp.runtime;
                    else
                        r_temp.runtime = inf;
                    end
                else % too slow
                    if plot_fig && ishandle(199 + ii)
                        set(groot, 'CurrentFigure', 199 + ii);
                        r_plot = plot(r_temp.timedata(:, 2) / 1000, r_temp.timedata(:, 3) * 3.6, 'r');
                        drawnow;
                    end
                    routedata = OffsetValues(par, val * delta, routedata);
                    echoPars(par, routedata, data.echo);
                    [r_temp, time_arr, thermal_out]  = ShSim_CalcSeg(trackdata(ii),vehicle, sim_config, time_dep, thermal_in, routedata, gradientdata, tdata, 0);
                    if data.echo
                        fprintf([' => Time: ', num2str(r_temp.runtime), ' (diff=', num2str(target_runtime - (r_temp.runtime + delay)), ' s)\n'])
                    end
                    if oldmax == r_temp.runtime
                        max_found = 1;
                    elseif r_temp.runtime > 0
                        oldmax = r_temp.runtime;
                    else
                        r_temp.runtime = inf;
                    end
                end
                delta = delta / 2;
            end
            
        else
            r_temp = r_org;
        end
        val = GetValues(par, routedata);
        for jj = 1:size(par, 2)
            r_temp.(par{jj}) = val(jj);
        end
        slask = 1;
        if slask
            if data.echo
                fprintf(['Time: ', num2str(r_temp.runtime), ' s, Energy: ', num2str((r_temp.w_traction + r_temp.w_regen) / 1000, '%5.2f'), ' kWh', newline])
            end
            for jj = 1:size(par, 2) + 2
                locVehicle = vehicle;
                troutedata = routedata;
                switch jj
                    case 1 % Effort
                        locVehicle.maxTE = locVehicle.maxTE * .9;
                    case 2 % Power
                        locVehicle.maxPwrTEbase = locVehicle.maxPwrTEbase * .9;
                    otherwise
                        troutedata.(par{jj - 2}) = val(jj - 2) * .99;
                end
                [r_temp2, ~, ~]  = ShSim_CalcSeg(trackdata(ii), locVehicle, sim_config, time_dep, thermal_in, troutedata, gradientdata, tdata, 0);
                if data.echo
                    fprintf(['Time: ', num2str(r_temp2.runtime), ' s, Energy: ', num2str((r_temp2.w_traction + r_temp2.w_regen) / 1000, '%5.2f'), ' kWh, relative saving: '])
                end
                sTime = r_temp.runtime - r_temp2.runtime;
                sEnergy = r_temp2.w_traction + r_temp2.w_regen - r_temp.w_traction - r_temp.w_regen;
                dw = sEnergy / sTime * 1000;
                if data.echo
                    fprintf([num2str(dw), ' Wh/s', newline])
                end
                if jj == 1
                    r_temp.('fMax_dw') = dw;
                elseif jj == 2
                    r_temp.('pMax_dw') = dw;
                else
                    r_temp.([par{jj - 2}, '_dw']) = dw;
                end
            end
        end
        routedata = SetValues(par, orgval, routedata);
        delay = r_temp.runtime - target_runtime + delay;
        if data.echo
            fprintf(['Current delay is: ', num2str(delay,'%5.2f'), ' seconds\n'])
        end
        if plot_fig && ishandle(199 + ii)
            set(groot, 'CurrentFigure', 199 + ii);
            g_plot = plot(r_temp.timedata(:, 2) / 1000, r_temp.timedata(:, 3) * 3.6, 'g', 'LineWidth', 2);
            if exist('b_plot', 'var') || exist('r_plot', 'var')
                if ~isempty(b_plot) && ~isempty(r_plot)
                    legend([k_plot, b_plot, r_plot, g_plot], 'Original run', 'Too fast', 'Too slow', 'Final optimisation', 'Location', 'best');
                elseif ~isempty(b_plot)
                    legend([k_plot, b_plot, g_plot], 'Original run', 'Too slow', 'Final optimisation', 'Location', 'best');
                elseif ~isempty(r_plot)
                    legend([k_plot, r_plot, g_plot], 'Original run', 'Too fast', 'Final optimisation', 'Location', 'best');
                else
                    legend([k_plot, g_plot], 'Original run', 'Final optimisation', 'Location', 'best');
                end
            else
                legend([k_plot, g_plot], 'Original run', 'Final optimisation', 'Location', 'best');
            end
            drawnow;
        end
    else
        r_temp = r_org;
    end
    r_temp = find_coasting(r_temp);
    locResult(index) = r_temp; %#ok<AGROW>
    time_dep = time_arr;
    thermal_in = thermal_out;
    if time_arr == 0
        break
    end
end
if lastFig && ishandle(lastFig)
    close(lastFig)
end
if ishandle(hwb)
    close(hwb)
end
if ishandle(hwb0)
    close(hwb0)
end
if data.echo
    fprintf(['Total calculation time: ', num2str(round(seconds(datetime - start_time_stamp) * 10) / 10), ' s\n']);
end
for ii = 1:size(locResult, 2)
    locResult(ii).timedata(:, 2) = locResult(ii).timedata(:, 2) - (sim_config.plotRef - 1) / 2 * vehicle.unitLength * vehicle.nUnits;
end
locResult = ShSim_sum(locResult);
return %ShSim_Calc

function val = GetValues(par, data)
val(size(par, 2)) = 0;
for ii = 1:size(par, 2)
    val(ii) = data.(par{ii});
end
return % GetValues

function data_out = OffsetValues(par, offsets, data)
for ii = 1:size(par, 2)
    data.(par{ii}) = data.(par{ii}) + offsets(ii); 
end
data_out = data;
return % OffsetValues

function data_out = SetValues(par, val, data)
for ii = 1:size(par, 2)
    data.(par{ii}) = val(ii);
end
data_out = data;
return % SetValues

function echoPars(par, data, echo)
if echo
    for ii = 1:size(par, 2)
        fprintf(['''', par{ii}, ''' = ', num2str(data.(par{ii})), '\n'])
    end
end
return % echoPars

function r = find_coasting(r)
r.ATO_dCoastStart = [];
r.ATO_dCoastEnd = [];
c = find(abs(r.timedata(:, 4)) > 0.01, 1, 'first');
b = find(abs(r.timedata(c:end, 4)) < 0.01, 1, 'first') + c - 1;
while b
    c = find(abs(r.timedata(b:end, 4)) > 0.01, 1, 'first') + b - 1;
    if r.timedata(c, 1) - r.timedata(b, 1) > 5
        r.ATO_dCoastStart = [r.ATO_dCoastStart, r.timedata(b, 2)];
        r.ATO_dCoastEnd = [r.ATO_dCoastEnd, r.timedata(c, 2)];
    end
    b = find(abs(r.timedata(c:end, 4)) < 0.01, 1, 'first') + c - 1;
end

return % find_coasting

