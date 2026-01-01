function nFigs = ShSim_Plot_Result(timedata, vehicle, input_data, type, figs, dispName, sim_config)
% ShSim_Plot_Result
% 0.1.1.1
%
%   Desciption
%       This function creates the plots deined in figs for the data
%       defined by 'timedata', and using the vehicle data
%       sim_config
%           struct 'vehicle' with 
%           nUnits
%           uNorm
%           rLineInd
%           ILDriveMax
%           ILBrakeMax
%           BRThermC
%       matrix 'timedata' columns
%           01 Time
%           02 Dist
%           03 Speed
%           04 TE/DBE
%           05 FricBE
%           06 ULine
%           07 ILine
%           08 BRPow
%           09 BRTemp
%           10 iMotor
%           11 mSlip
%           12 Modulation index
%           13 ULink
%           14 ESPow - Energy Saver Power
%           15 Main Transformer prim winding RMS
%           16 Main Transformer 2nd winding RMS
%           17 Line filter RMS
%           18 MCM phase current RMS
%           19 LCM phase current RMS
%           20 Traction motor RMS
%           21 Spare
%           22 Spare
%           23 Spare
%           24 Spare
%           25 Mode of operation
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-29  JR      Adjusted to new ShogunSim environment
%
%
nFigs = [];
if type == 0 || isempty(timedata)
    return
end
switch type
    case 1 % Foreground
        colb = [0, 0, 1]; % blue
        colr = [1, 0, 0]; % red
        colg = [0, 1, 0]; % green
        colk = [0, 0, 0]; % black
        lSize = 2;
    case 2 % Alternativ
        colb = [0.75, 0.75, 1];    % blue
        colr = [1, 0.75, 0.75];    % red
        colg = [0.75, 1, 0.75];    % green
        colk = [0.25, 0.25, 0.25]; % black
        lSize = 3;
    case 3 % Background
        colb = [0.75, 0.75, 0.75]; % blue
        colr = [0.75, 0.75, 0.75]; % red
        colg = [0.75, 0.75, 0.75]; % green
        colk = [0.75, 0.75, 0.75]; % black
        lSize = 3;
    case 4 % Importdata
        colb = [1, 0, 1];       % magnenta
        colr = [0, 1, 1];       % cyan
        colg = [1, 1, 0];       % yellow
        colk = [0.5, 0.5, 0.5]; % black
        lSize = 2;
    case 5 % Maxdata
        colb = [1, 0.5, 1];     % magnenta
        colr = [0.5, 1, 1];     % cyan
        colg = [1, 1, 0.5];     % yellow
        colk = [0.75, 0, 0.75]; % black
        lSize = 2;
end
% scrn= get(0, 'ScreenSize');
[sections, route, profile] = ShSim_Read_Route(input_data.track_file, vehicle, sim_config); %#ok<ASGLU>
offset = [0, 0];
if figs(1) % Speed vs Distance
    nFigs = [nFigs, 100];
    stdfig([0.5, 0.5], 100, 'Speed vs Distance', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2); 
    offset(2) = mod(offset(2) + 1, 2);  
    hold on;
    grid on;
    plot3(timedata(:, 2) / 1000, timedata(:, 3) * 3.6, timedata(:, 1), 'x-', 'Color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
    title('Speed vs distance');
    xlabel('Distance [km]');
    ylabel('Speed [km/h]');
    zlabel('Time [s]');
    v = axis;
    ylim([0, max([1.1 * max(timedata(:, 3) * 3.6), v(4)])]);
    zoom on;
end

if figs(2) % Speed vs Time
    nFigs = [nFigs, 101];
    stdfig([0.5, 0.5], 101, 'Speed vs Time', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    hold on;
    grid on;
    if timedata(1, 1) > 0
        plot3(timedata(:, 1) / 3600, timedata(:, 3) * 3.6, timedata(:, 2) / 1000, 'x-', 'Color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot3(timedata(:, 1) / 3600, timedata(:, 26) * 3.6, timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average');
        end
        if input_data.detorial % Stupid word, but result is not a merge of several simulations
            for ii = 1:size(sections, 2)
                if sections(ii).departure > 0 && type == 1
                    h = plot(sections(ii).departure / 3600 * ones(1, 5), (0:4) * vehicle.vMax * 3.6 / 5, '-gx', 'LineWidth', 3, 'DisplayName', 'Departure time');
                    if ii > 1
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
                end
                if sections(ii).arrival > 0 && type == 1
                    h = plot(sections(ii).arrival / 3600 * ones(1, 5), (0:4) * vehicle.vMax * 3.6 / 5, '-rx', 'LineWidth', 3, 'DisplayName', 'Arrival time');
                    if ii > 1
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
                end
            end
        end
        xlabel('Time [hrs]');
    else
        plot3(timedata(:, 1), timedata(:, 3) * 3.6, timedata(:, 2) / 1000, 'x-', 'Color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot3(timedata(:, 1), timedata(:, 26) * 3.6, timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average');
        end
        for ii=1:size(sections, 2)
            if sections(ii).departure<timedata(end, 1)
                if sections(ii).departure > 0 && type == 1
                    h = plot(sections(ii).departure * ones(1, 5), (0:4) * vehicle.vMax * 3.6 / 5, '-gx', 'LineWidth', 1, 'DisplayName', 'Departure time');
                    if ii < size(sections, 2)
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
                end
                if sections(ii).arrival > 0 && type == 1
                    h = plot(sections(ii).arrival * ones(1, 5), (0:4) * vehicle.vMax * 3.6 / 5, '-rx', 'LineWidth', 1, 'DisplayName', 'Arrival time');
                    if ii > 1
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
                end
            end
        end
        xlabel('Time [s]');
    end
    title('Speed vs Time');
    v = axis;
    ylim([0, max([1.1 * max(timedata(:, 3) * 3.6) v(4)])])
    ylabel('Speed [km/h]');
    zlabel('Dist [km]');
    zoom on;
end
if figs(3) % TE/DBE/FBrE vs distance
    nFigs = [nFigs, 102];
    stdfig([0.5, 0.5], 102, 'TE/DBE/FBrE vs distance', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);  
    hold on;
    grid on;
    plot3(timedata(:, 2) / 1000, -timedata(:, 5), timedata(:, 1), 'x-', 'Color', colr, 'LineWidth', lSize, 'DisplayName', [dispName, ' - Friction']);
    plot3(timedata(:, 2) / 1000, timedata(:, 4), timedata(:, 1), 'x-', 'Color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Electric']);
    title('TE/DBE/FBrE vs distance');
    xlabel('Distance [km]');
    ylabel('Effort [kN]');
    zlabel('Time [s]');
    zoom on;
end
if figs(4) % TE/DBE/FBrE vs Time
    nFigs = [nFigs, 103];
    stdfig([0.5, 0.5], 103, 'TE/DBE/FBrE vs Time', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);  
    hold on;
    grid on;
    if timedata(1, 1) > 0
        plot3(timedata(:, 1) / 3600, -timedata(:, 5), timedata(:, 2) / 1000, 'x-', 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Friction']);
        plot3(timedata(:, 1) / 3600, timedata(:, 4), timedata(:, 2) / 1000, 'x-', 'Color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Electric']);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot3(timedata(:, 1) / 3600, timedata(:, 27), timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average');
        end
        xlabel('Time [hrs]');
    else
        plot3(timedata(:, 1), -timedata(:, 5), timedata(:, 2) / 1000, 'x-', 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Friction']);
        plot3(timedata(:, 1), timedata(:, 4), timedata(:, 2) / 1000, 'x-', 'Color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Electric']);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot3(timedata(:, 1), timedata(:, 27), timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average');
        end
        xlabel('Time [s]');
    end
    title('TE/DBE/FBrE vs Time');
    ylabel('Effort [kN]');
    zlabel('Distance [km]');
    zoom on;
end
if figs(5) % Line voltage and current vs distance
    nFigs = [nFigs, 104];
    ax1 = stdfig([0.5, 0.5], 104, 'Line voltage and current vs distance', [2, 1, 1], offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    hold on;
    grid on;
    title('Line voltage vs distance');
    if vehicle.U_Line_Nom >= 1000 % bigger than 1000V - scale in kV
        plot([timedata(1, 2), timedata(end, 2)] / 1000, [vehicle.U_Line_Nom, vehicle.U_Line_Nom] / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Nominal');
        if vehicle.rLineInd > 0
            plot(timedata(:, 2) / 1000, timedata(:, 6) / 1000 - timedata(:, 7) * vehicle.rLineInd / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        end
        plot(timedata(:, 2) / 1000, timedata(:, 6) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        ylabel('Line voltage [kV]');
    else
        plot([timedata(1, 2), timedata(end, 2)] / 1000, [vehicle.U_Line_Nom, vehicle.U_Line_Nom], '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', dispName);
        if vehicle.rLineInd > 0
            plot(timedata(:, 2) / 1000, timedata(:, 6) - timedata(:, 7) * vehicle.rLineInd, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        end
        plot(timedata(:, 2) / 1000, timedata(:, 6), 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        ylabel('DC-link voltage [V]');
    end
    zoom on;
    ax2 = stdfig([0.5, 0.5], 104, 'Line voltage and current vs distance', [2, 1, 2], offset);
    offset(1)=mod(offset(1)+offset(2), 2); offset(2)=mod(offset(2)+1, 2); 
%     ax2=subplot(2, 1, 2);
    hold on;
    grid on;
    
    title('Line current vs distance');
    if ~isinf(vehicle.ILDriveMax)
        plot([timedata(1, 2) timedata(end, 2)]/1000, [vehicle.ILDriveMax vehicle.ILDriveMax], '--', 'color', colb, 'LineWidth', lSize, 'DisplayName', 'Max limit');
    end
    if ~isinf(vehicle.ILBrakeMax)
        plot([timedata(1, 2) timedata(end, 2)]/1000, [-vehicle.ILBrakeMax -vehicle.ILBrakeMax], '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Min limit');
    end
    plot(timedata(:, 2)/1000, timedata(:, 7), 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
    ylabel('Line current [A]');
    xlabel('Distance [km]');
    zoom on;
    linkaxes([ax1 ax2], 'x');
end
if figs(6) % Line voltage and current vs Time
    nFigs = [nFigs, 105];
    ax1 = stdfig([0.5, 0.5], 105, 'Line voltage and current vs Time', [2, 1, 1], offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 

    hold on;
    grid on;
    title('Line voltage vs Time');
    if vehicle.U_Line_Nom >= 1000
        if timedata(1, 1) > 0
            plot([timedata(1, 1), timedata(end, 1)] / 3600, [vehicle.U_Line_Nom, vehicle.U_Line_Nom] / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Nominal');
            if vehicle.rLineInd>0
                plot(timedata(:, 1) / 3600, timedata(:, 6) / 1000 - timedata(:, 7) * vehicle.rLineInd / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
            end
            plot(timedata(:, 1) / 3600, timedata(:, 6) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        else
            plot([timedata(1, 1), timedata(end, 1)], [vehicle.U_Line_Nom, vehicle.U_Line_Nom] / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Nominal');
            if vehicle.rLineInd > 0
                plot(timedata(:, 1), timedata(:, 6) / 1000 - timedata(:, 7) * vehicle.rLineInd / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
            end
            plot(timedata(:, 1), timedata(:, 6) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        end
        ylabel('Line voltage [kV]');
    else
        if timedata(1, 1) > 0
            plot([timedata(1, 1), timedata(end, 1)] / 3600, [vehicle.U_Line_Nom, vehicle.U_Line_Nom], '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', dispName);
            if vehicle.rLineInd > 0
                plot(timedata(:, 1) / 3600, timedata(:, 6) - timedata(:, 7) * vehicle.rLineInd, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
            end
            plot(timedata(:, 1) / 3600, timedata(:, 6), 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        else
            plot([timedata(1, 1), timedata(end, 1)], [vehicle.U_Line_Nom, vehicle.U_Line_Nom], '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', dispName);
            if vehicle.rLineInd > 0
                plot(timedata(:, 1), timedata(:, 6) - timedata(:, 7) * vehicle.rLineInd, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
            end
            plot(timedata(:, 1), timedata(:, 6), 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        end
        ylabel('DC-link voltage [V]');
    end
    zoom on;
    ax2 = stdfig([0.5, 0.5], 105, 'Line voltage and current vs Time', [2, 1, 2], offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    title('Line current vs Time');
    if ~isinf(vehicle.ILDriveMax)
        if timedata(1, 1) > 0
            plot([timedata(1, 1), timedata(end, 1)] / 3600, [vehicle.ILDriveMax, vehicle.ILDriveMax], '--', 'color', colb, 'LineWidth', lSize, 'DisplayName', 'Max limit');
        else
            plot([timedata(1, 1), timedata(end, 1)], [vehicle.ILDriveMax, vehicle.ILDriveMax], '--', 'color', colb, 'LineWidth', lSize, 'DisplayName', 'Max limit');
        end
    end
    if ~isinf(vehicle.ILBrakeMax)
        if timedata(1, 1) > 0
            plot([timedata(1, 1), timedata(end, 1)] / 3600, [-vehicle.ILBrakeMax, -vehicle.ILBrakeMax], '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Min limit');
        else
            plot([timedata(1, 1), timedata(end, 1)], [-vehicle.ILBrakeMax, -vehicle.ILBrakeMax], '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Min limit');
        end
    end
    if timedata(1, 1) > 0
        plot(timedata(:, 1) / 3600, timedata(:, 7), 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot(timedata(:, 1) / 3600, timedata(:, 29) / 1000, 'kx-', 'DisplayName', 'Average');
        end
        xlabel('Time [hrs]');
    else
        plot(timedata(:, 1), timedata(:, 7), 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot(timedata(:, 1), timedata(:, 29), 'kx-', 'DisplayName', 'Average');
        end
        xlabel('Time [s]');
    end
    ylabel('Line current [A]');
    zoom on;
    linkaxes([ax1, ax2], 'x');
end
if figs(7) % Line- and Rail power vs distance
    nFigs = [nFigs, 106];
    stdfig([0.5, 0.5], 106, 'Line- and Rail power vs distance', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    hold on;
    grid on;
    plot3(timedata(:, 2) / 1000, timedata(:, 6) .* timedata(:, 7) / 1000, timedata(:, 1), 'color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Line']);
    plot3(timedata(:, 2) / 1000, timedata(:, 3) .* timedata(:, 4), timedata(:, 1), 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Rail']);
    title('Line & Rail power vs distance');
    ylabel('Line/Rail Power [kW]');
    xlabel('Distance [km]');
    zlabel('Time [s]');
    zoom on;
end
if figs(8) % Line- and Rail power vs Time
    nFigs = [nFigs, 107];
    stdfig([0.5, 0.5], 107, 'Line- and Rail power vs Time', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    hold on;
    grid on;
    if timedata(1, 1) > 0
        plot3(timedata(:, 1) / 3600, timedata(:, 6) .* timedata(:, 7) / 1000, timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Line']);
        plot3(timedata(:, 1) / 3600, timedata(:, 3) .* timedata(:, 4), timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Rail']);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot3(timedata(:, 1) / 3600, timedata(:, 26) .* timedata(:, 27), timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average - Rail');
            plot3(timedata(:, 1) / 3600, timedata(:, 6) .* timedata(:, 29) / 1000, timedata(:, 2) / 1000, 'ko-', 'DisplayName', 'Average - Line');
        end
        xlabel('Time [hrs]');
    else
        plot3(timedata(:, 1), timedata(:, 6) .* timedata(:, 7) / 1000, timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Line']);
        plot3(timedata(:, 1), timedata(:, 3) .* timedata(:, 4), timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Rail']);
        if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
            plot3(timedata(:, 1), timedata(:, 6) .* timedata(:, 29) / 1000, timedata(:, 2) / 1000, 'ko-', 'DisplayName', ' - Line');
            plot3(timedata(:, 1), timedata(:, 26) .* timedata(:, 27), timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average - Line');
        end
        xlabel('Time [s]');
    end
    title('Line & Rail power power vs time');
    zlabel('Distance [km]');
    ylabel('Line/Rail Power [kW]');
    zoom on;
end
if figs(9) || figs(10)
    dTemp = diff([0; timedata(:, 9)]);  % Temp
    dTime=diff([0;timedata(:, 1)]);     % Time
    if isfield(vehicle, 'BRThermC')
        cPow=timedata(:, 8) - dTemp * vehicle.BRThermC ./ dTime * vehicle.nINV * vehicle.nHVSystems;
    end
end
if figs(9)
    nFigs = [nFigs, 108];
    stdfig([0.5, 0.5], 108, 'Brake resistor Power & Temperature vs distance', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    ax1 = subplot(2, 1, 1);
    hold on;
    grid on;
    title('Brake resistor power vs distance');
    plot(timedata(:, 2) / 1000, timedata(:, 8) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - BRePow']);
    if exist('cPow', 'var')
        plot(timedata(:, 2) / 1000, cPow / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', 'Cooling power');
    end
    ylabel('Power [kW]');
    zoom on;

    ax2 = subplot(2, 1, 2);
    hold on;
    grid on;
    title('Brake resistor temperature vs distance');

    plot(timedata(:, 2)/1000, timedata(:, 9), 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
    ylabel('BR temperature [C]');
    xlabel('Distance [km]');
    zoom on;
    linkaxes([ax1, ax2], 'x');
end
if figs(10)
    nFigs = [nFigs, 109];
    stdfig([0.5, 0.5], 109, 'Brake resistor Power & Temperature vs Time', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    ax1 = subplot(2, 1, 1);
    hold on;
    grid on;
    title('Brake resistor power vs time');
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 8) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - BRePow']);
    if exist('cPow', 'var') 
        plot(timedata(:, 1) - timedata(1, 1), cPow / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', 'Cooling power');
    end
    ylabel('Power [kW]');
    zoom on;
    ax2 = subplot(2, 1, 2);
    hold on;
    grid on;
    title('Brake resistor temperature vs Time');
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 9), 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
    ylabel('BR temperature [\circC]');
    xlabel('Time [s]');
    zoom on;
    linkaxes([ax1 ax2], 'x');
end
if figs(11)
    nFigs = [nFigs, 110];
    stdfig([0.5, 0.5], 110, 'Distance vs Time', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    hold on;
    grid on;
    title('Distance vs time');
    if timedata(1, 1) > 0
        plot(timedata(:, 1) / 3600, timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
    else
        plot(timedata(:, 1) - timedata(1, 1), timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
    end
    ylabel('Distance [km]');
    xlabel('Time [s]');
    zoom on;
end
if figs(12)
    nFigs = [nFigs, 111];
    stdfig([0.5, 0.5], 111, 'Effort vs Speed', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    title('Effort vs speed');
    plot3(timedata(:, 3) * 3.6, timedata(:, 4), timedata(:, 2) / 1000, 'x-', 'color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Electric']);
    plot3(timedata(:, 3) * 3.6, -timedata(:, 5), timedata(:, 2) / 1000, 'x-', 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Friction']);
    if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
        plot3(timedata(:, 26) * 3.6, timedata(:, 27), timedata(:, 2) / 1000, 'ko-', 'DisplayName', 'Average - Electric');
        plot3(timedata(:, 26) * 3.6, -timedata(:, 28), timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average - Friction');
    end
    ylabel('Effort [kN]');
    xlabel('Speed [km/h]');
    zlabel('Distance [km]');
    zoom on;
end
if figs(13)
    nFigs = [nFigs, 112];
    stdfig([0.5, 0.5], 112, 'Acceleration vs Distance', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    t = timedata(:, 1);
    a=[0; diff(timedata(:, 3)) ./ diff(t)];
    d = timedata(:, 2) / 1000;
    p = find(abs(a) == inf);
    while ~isempty(p)
        p = p(1);
        t = [t(1:p - 1); t(p + 1:end)];        
        a = [a(1:p - 1); a(p + 1:end)];
        d = [d(1:p - 1); d(p + 1:end)];
        p = find(abs(a) == inf);
    end
    
    plot(d, a, 'x-', 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
    title('Acceleration vs Distacce');
    xlabel('Distance [km]');
    ylabel('Acceleration [m/s^2]');
    zoom on;
end
if figs(14)
    nFigs = [nFigs, 113];
    stdfig([0.5, 0.5], 113, 'Acceleration vs Speed', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    hold on;
    grid on;
    t = timedata(:, 1);
    d = timedata(:, 2) / 1000;
    a = [0; diff(timedata(:, 3)) ./ diff(t)];
    v = timedata(:, 3);
    p = find(abs(a) == inf);
    while ~isempty(p)
        a(p) = a(p-1);
        p = find(abs(a) == inf);
    end
    plot3(v * 3.6, a, d, 'x-', 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
    x = xlim;
    xlim([0, x(2)]);
    title('Acceleration vs Speed');
    xlabel('Speed [km/h]');
    ylabel('Acceleration [m/s^2]');
    zlabel('Distance [km]');
    zoom on;
end
if figs(15)
    nFigs = [nFigs, 114];
    v = timedata(:, 3) * 3.6;
    vmax = max(v);
    ax1 = stdfig([0.5, 0.5], 114, 'Line voltage and current vs Speed', [2, 1, 1], offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    hold on;
    grid on;
    title('Line voltage and current vs Speed');
    if vehicle.U_Line_Nom >= 1000
        plot([0, vmax], [vehicle.U_Line_Nom, vehicle.U_Line_Nom] / 1000, '--', 'color', colg, 'LineWidth', lSize+1, 'DisplayName', 'Nominal');
        if vehicle.rLineInd > 0
            plot3(v, timedata(:, 6) / 1000 - timedata(:, 7) * vehicle.rLineInd / 1000, timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        end
        plot3(v, timedata(:, 6) / 1000, timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        ylabel('Voltage [kV]');
    else
        plot([0, vmax], [vehicle.U_Line_Nom, vehicle.U_Line_Nom], 'color', colg, 'LineWidth', lSize, 'DisplayName', dispName);
        if vehicle.rLineInd > 0
            plot3(v, timedata(:, 6) - timedata(:, 7) * vehicle.rLineInd, timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
        end
        plot3(v, timedata(:, 6), timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
        ylabel('Voltage [V]');
    end
    x = xlim;
    xlim([0, x(2)]);
    zlabel('Dist [km]');
    xlabel('Speed [km/h]');
    zoom on;

    ax2 = stdfig([0.5, 0.5], 114, 'Line voltage and current vs Speed', [2, 1, 2], offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    if vehicle.ILDriveMax ~= 0 && type == 1
        plot([0, vmax], [vehicle.ILDriveMax, vehicle.ILDriveMax], '--', 'color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Max Limit']);
    end
    if vehicle.ILBrakeMax ~= 0 && type == 1
        plot([0, vmax], [-vehicle.ILBrakeMax, -vehicle.ILBrakeMax], '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', [dispName ' - Min Limit']);
    end
    plot3(v, timedata(:, 7), timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', dispName);
    if size(timedata, 2) >= 24 && vehicle.ES_Op_Mode && size(timedata, 2) > 14
        plot3(v, timedata(:, 14) ./ timedata(:, 6) * vehicle.nMotors * vehicle.nHVSystems * vehicle.nINV, timedata(:, 2) / 1000, 'b', 'DisplayName', '???');
    end
    x = xlim;
    xlim([0, x(2)]);
    ylabel('Line current [A]');
    xlabel('Speed [km/h]');
    zlabel('Dist [km]');
    zoom on;
    linkaxes([ax1, ax2], 'x');
end
if figs(16)
    nFigs = [nFigs, 115];
    v = timedata(:, 3) * 3.6;
    stdfig([0.5, 0.5], 115, 'Line- and Rail power vs Speed', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    title('Line/Rail power vs Speed');
    plot3(v, timedata(:, 6) .* timedata(:, 7) / 1000, timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Line']);
    plot3(v, timedata(:, 3) .* timedata(:, 4), timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Rail']);
    if size(timedata, 2) > 30 && isfield(sim_config, 'plotAverageVFI') && sim_config.plotAverageVFI
        plot3(v, timedata(:, 29) .* timedata(:, 7) / 1000, timedata(:, 2) / 1000, 'ko-', 'DisplayName', 'Average - Line');
        plot3(v, timedata(:, 26) .* timedata(:, 27), timedata(:, 2) / 1000, 'kx-', 'DisplayName', 'Average - Rail');
    end
    x = xlim;
    xlim([0, x(2)]);
    ylabel('Power [kW]');
    xlabel('Speed [km/h]');
    zlabel('Dist [km]');
    zoom on;
end
if figs(17)
    nFigs = [nFigs, 116];
    v = timedata(:, 3) * 3.6;
    ax1 = stdfig([0.5, 0.5], 116, 'Brake resistor Power & Temperature vs Speed', [2, 1, 1], offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    title('Brake resistor power vs Speed');
    plot3(v, timedata(:, 8) / 1000, timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - BR Power']);
    dTemp = diff([0; timedata(:, 9)]);
    dTime = diff([0; timedata(:, 1)]);
    if isfield(vehicle, 'BRThermC')
        cPow = timedata(:, 8) / 1000 - dTemp * vehicle.BRThermC ./ dTime;
        plot3(v, cPow, timedata(:, 2) / 1000, 'color', colb, 'LineWidth', lSize, 'DisplayName', [dispName ' - Cooling power']);
    end
    zlabel('Dist [km]');
    ylabel('Power [kW]');
    xlabel('Speed [km/h]');
    x=xlim;
    xlim([0, x(2)]);
    zoom on;

    ax2 = stdfig([0.5, 0.5], 116, 'Brake resistor Power & Temperature vs Speed', [2, 1, 2], offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    title('Brake resistor temperature vs distance');
    plot3(v, timedata(:, 9), timedata(:, 2) / 1000, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName ' - Temp']);
    x=xlim;
    xlim([0, x(2)]);
    ylabel('BR temperature [C]');
    xlabel('Speed [km/h]');
    zlabel('Dist [km]');
    zoom on;
    linkaxes([ax1, ax2], 'x');
end
if figs(18) % floating average speed
    v_avg = [0; (timedata(2:end, 2) - timedata(1, 2)) ./ (timedata(2:end, 1) - timedata(1, 1)) * 3.6];
    if figs(1) %speed vs distance
        figure(100);
        plot3(timedata(:, 2) / 1000, v_avg, timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Floating average speed');
    end
    if figs(2) % Speed vs time
        figure(101);
        if timedata(1, 1) > 0
            plot3(timedata(:, 1) / 3600, v_avg, timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Floating average speed');
        else
            plot3(timedata(:, 1), v_avg, timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Floating average speed');
        end
    end
    if ~(figs(1) || figs(2))
        nFigs = [nFigs, 117];
        stdfig([0.5, 0.5], 117, 'Floating average speed vs Time', 1, offset);
        offset(1) = mod(offset(1) + offset(2), 2);
        offset(2) = mod(offset(2) + 1, 2);
        grid on;
        hold on;
        title('Floating average vs Time');
        plot3(timedata(:, 1) - timedata(1, 1), v_avg, timedata(:, 2) / 1000, 'color', colk, 'LineWidth', lSize, 'DisplayName', dispName);
        ylabel('Speed [km/h]');
        xlabel('Time [s]');
        zlabel('Dist [km]');
    end
end
if figs(19) || figs(20)
    if size(timedata, 2) > 10
        speed = timedata(:, 26);
    else
        speed = timedata(:, 3);
    end
    MotorRPM = speed * vehicle.gearRatio / (vehicle.wheelsize(vehicle.ws_selection) * pi) * 60; % rpm
end
if figs(19) % Torque vs. Motor speed
    TE = timedata(:, 4);
    gearloss = vehicle.gearLoss0 * vehicle.gearRatio * vehicle.nHVSystems * vehicle.nINV * 2 / vehicle.wheelsize(vehicle.ws_selection) / 1000 + vehicle.gearLoss1 * TE / 100 + MotorRPM * vehicle.gearLoss2;
    TE = TE + gearloss;
    nMotors = vehicle.nMotors * vehicle.nINV * vehicle.nHVSystems;
    tq = TE * vehicle.wheelsize(vehicle.ws_selection) * 1000 / 2 / vehicle.gearRatio / nMotors;
    nFigs = [nFigs, 118];
    stdfig([0.5, 0.5], 118, 'Motor torque vs. motor speed', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    title('Torque vs. Motor speed');
    plot3(MotorRPM, tq, timedata(:, 2) / 1000, 'x-', 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
    ylabel('Motor Torque [Nm]');
    xlabel('Motor Speed [rpm]');
end
if figs(20) && size(timedata, 2) > 10 % Inverter current vs. Motor speed
    nFigs = [nFigs, 119];
    stdfig([0.5, 0.5], 119, 'Inverter current vs. motor speed', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2);
    title('Torque vs. Motor speed');
    hold on;
    grid on;
    iStator = timedata(:, 10);
    plot3(MotorRPM, iStator, timedata(:, 2) / 1000, 'x-', 'color', colb, 'LineWidth', lSize, 'DisplayName', dispName);
    ylabel('Inverter Current [A]');
    xlabel('Motor Speed [rpm]');
end

if figs(21) && size(timedata, 2)>10 %Thermal limits Line voltage system
    nFigs = [nFigs, 120];
    stdfig([0.5, 0.5], 120, 'Thermal Limits Line Voltage System vs. Time', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    ax1 = subplot(2, 1, 1);
    hold on;
    grid on;
    if vehicle.uFreq > 0 % AC system
        plot(timedata(:, 1) - timedata(1, 1), timedata(:, 15) / vehicle.IPrimTrafoCont * 100, 'color', colb, 'LineWidth', lSize, 'DisplayName', 'IPrim');
        plot(timedata(:, 1) - timedata(1, 1), timedata(:, 16) / vehicle.I2ndTrafoCont * 100, 'color', colr, 'LineWidth', lSize, 'DisplayName', 'ISec');
        plot(timedata(:, 1) - timedata(1, 1), timedata(:, 19) / vehicle.IphLCMCont * 100, 'color', colg, 'LineWidth', lSize, 'DisplayName', 'IPhLCM');
    else
        plot(timedata(:, 1) - timedata(1, 1), timedata(:, 17) / vehicle.IDCFilterCont * 100, 'color', colb, 'LineWidth', lSize, 'DisplayName', 'DC Line filter');
    end
    title('Thermal Limits Line Voltage System vs. Time');
    ylabel('Load level [%]');
    xlabel('Time [s]');
    zoom on;
    ax2 = subplot(2, 1, 2);
    hold on;
    grid on;
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 21) * 100, 'color', colb, 'LineWidth', lSize+1, 'DisplayName', 'Power limitation');
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 22) * 100, 'color', colr, 'LineWidth', lSize, 'DisplayName', 'Effort limitation');
    title('Thermal utilisation vs. Time');
    ylabel('Utilisation [%]');
    xlabel('Time [s]');
    zoom on;
    linkaxes([ax1, ax2], 'x');
end
if figs(22) && size(timedata, 2) > 10 % Thermal limits Drive system
    nFigs = [nFigs, 121];
    stdfig([0.5, 0.5], 121, 'Thermal Limits Drive System vs. Time', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    ax1 = subplot(2, 1, 1);
    hold on;
    grid on;
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 18) / vehicle.IphMCMContHex * 100, 'color', colb, 'LineWidth', lSize+1, 'DisplayName', 'MCM phase current');
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 20) / vehicle.IstatTMCont * 100, 'color', colr, 'LineWidth', lSize, 'DisplayName', 'Traction motor');
    title('Thermal Limits Drive System vs. Time');
    ylabel('Load level [%]');
    xlabel('Time [s]');
    zoom on;
    ax2 = subplot(2, 1, 2);
    hold on;
    grid on;
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 21) * 100, 'color', colb, 'LineWidth', lSize+1, 'DisplayName', 'Power limitation');
    plot(timedata(:, 1) - timedata(1, 1), timedata(:, 22) * 100, 'color', colr, 'LineWidth', lSize, 'DisplayName', 'Effort limitation');
    title('Thermal utilisation vs. Time');
    ylabel('Utilisation [%]');
    xlabel('Time [s]');
    zoom on;
    linkaxes([ax1, ax2], 'x');
end
if (figs(23) || figs(24)) && size(timedata, 2) > 30
    dt = diff(timedata(:, 1));
    % Line power => energy
    w1 = [0; cumsum(timedata(2:end, 29) .* timedata(2:end, 6) / 3600000 .* dt)];
    % Rail power => energy
    w2 = [0; cumsum(timedata(2:end, 27) .* timedata(2:end, 26) / 3600 .* dt)];
    % Potential energy
    w3 = (vehicle.totMass * timedata(:, 3).^2) / 7200;
    if figs(23) % Energy vs Distance
        nFigs = [nFigs, 122];
        stdfig([0.5, 0.5], 122, 'Energy concumption vs. Distance', 1, offset);
        offset(1) = mod(offset(1) + offset(2), 2);
        offset(2) = mod(offset(2) + 1, 2);
        hold on;
        grid on;
        plot(timedata(:, 2) / 1000, w1, 'color', colb, 'LineWidth', lSize + 1, 'DisplayName', [dispName, ' - Line']);
        plot(timedata(:, 2) / 1000, w2, 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName, ' - Rail']);
        title('Energy consumption vs. Distance');
        ylabel('Energy [kWh]');
        xlabel('Distance [km]');
        zoom on;
    end
    if figs(24) % Energy vs Time
        nFigs = [nFigs, 123];
        stdfig([0.5, 0.5], 123, 'Energy concumption vs. Time', 1, offset);
        offset(1) = mod(offset(1) + offset(2), 2);
        offset(2) = mod(offset(2) + 1, 2); 
        hold on;
        grid on;
        plot(timedata(:, 1) - timedata(1, 1), w1, 'x-', 'color', colb, 'LineWidth', lSize+1, 'DisplayName', [dispName, ' - Line']);
        plot(timedata(:, 1) - timedata(1, 1), w2, 'o-', 'color', colr, 'LineWidth', lSize, 'DisplayName', [dispName, ' - Rail']);
        plot(timedata(:, 1) - timedata(1, 1), w3, 'ko--', 'LineWidth', lSize, 'DisplayName', 'Potential energy');
        title('Energy consumption vs. Time');
        ylabel('Energy [kWh]');
        xlabel('Time [s]');
        zoom on;
    end
end
if figs(25) && size(timedata, 2) > 10
    nFigs = [nFigs, 124];
    stdfig([0.5, 0.5], 124, 'Efficiency vs. Speed', 1, offset);
    offset(1) = mod(offset(1) + offset(2), 2);
    offset(2) = mod(offset(2) + 1, 2); 
    hold on;
    grid on;
    railPow = timedata(:, 26) .* timedata(:, 27);       % kW, average speed times average effort
    linePow = timedata(:, 6) .* timedata(:, 29) / 1000; % kW, line voltage times rms current
    b = find((railPow - vehicle.AuxPower) > 0.1 * max(railPow));
    eff = (railPow(b) + vehicle.AuxPower) ./ linePow(b) * 100;
    plot(timedata(b, 3) * 3.6, eff, '.', 'color', colb, 'LineWidth', lSize+1, 'DisplayName', [dispName ' - Traction']);
    b = find((railPow - vehicle.AuxPower) < 0.1 * min(railPow));
    eff = linePow(b) ./ (railPow(b) + vehicle.AuxPower) * 100;
    plot(timedata(b, 3) * 3.6, eff, '.', 'color', colr, 'LineWidth', lSize+1, 'DisplayName', [dispName ' - Braking']);
    title('Traction Efficiency vs. Speed');
    ylabel('Efficiency [%]');
    xlabel('Speed [km/h]');
    ylim([0, 100]);
    zoom on;
end
if figs(26) && size(timedata, 2) > 10
    if figs(3) % Effort vs distance, add notch
        figure(102);
        yyaxis right
        plot3(timedata(:, 2) / 1000, timedata(:, 31), timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Notch');
        ylabel('Notch [-]');
        ylim([-vehicle.nNotchB - 1, vehicle.nNotchT + 1]);
        yyaxis left
    end
    if figs(4) % Effort vs time, add notch
        figure(103);
        yyaxis right
        if timedata(1, 1) > 0
            plot3(timedata(:, 1) / 3600, timedata(:, 31), timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Notch');
        else
            plot3(timedata(:, 1), timedata(:, 31), timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Notch');
        end
        ylabel('Notch [-]');
        ylim([-vehicle.nNotchB - 1, vehicle.nNotchT + 1]);
        yyaxis left
    end
    if ~(figs(3) || figs(4))
        nFigs = [nFigs, 125];
        stdfig([0.5, 0.5], 125, 'Notches vs. Time', 1, offset);
        offset(1) = mod(offset(1) + offset(2), 2);
        offset(2) = mod(offset(2) + 1, 2); %#ok<NASGU>
        hold on;
        grid on;
        if timedata(1, 1) > 0
            plot3(timedata(:, 1) / 3600, timedata(:, 31), timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Notch');
        else
            plot3(timedata(:, 1), timedata(:, 31), timedata(:, 2) / 1000, '--', 'color', colg, 'LineWidth', lSize, 'DisplayName', 'Notch');
        end
        title('Notch selection vs. Time');
        ylabel('Notch [-]');
        xlabel('Time [s]');
        ylim([-vehicle.nNotchB - 1, vehicle.nNotchT + 1]);
        zoom on;
    end

end
return % ShSim_Plot_Result