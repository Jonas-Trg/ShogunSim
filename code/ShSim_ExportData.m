function xFile = ShSim_ExportData(input_data, r, dsmpl, prefix, sim_config, textformat, fInName, echo)
% ShSim_ExportData
% 0.1.1.1
%
%   Description
%       This function exports the result to either a mat file or a text
%       file.
%
%   Input signals
%       input_data - a struct with:
%           track_file - used to suggest the default name
%       r - result record struct
%       dsmpl - down sampling
%       prefix - string to amend the proposed file name with
%       sim_config - struct with
%           workingDirectory - the starting working directory
%       textformat - 1 is to export to text, otherwise to mat
%       fInName - default export name - if defined no uiputfile is needed
%       echo - in this context - show the progress bar
%
%   Output signals
%       xFile - the name of the exported file
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-30  JR      Code cleaning

global logfileref %#ok<GVMIS> 
if textformat
    fExt = '.txt';
else
    fExt = '.mat';
end
oldwd = cd(fullfile(sim_config.workingDirectory, 'Export'));
[~, fName] = fileparts(input_data.track_file);
if isempty(prefix)
    fName = [fName, '_', char(datetime, 'yyyyMMdd'), 'T', char(datetime, 'HHmmss'), fExt];
else
    fName = [prefix, '_', fName, '_', char(datetime, 'yyyyMMdd'), 'T', char(datetime, 'HHmmss'), fExt];
end
if isempty(fInName)
    [fName, p] = uiputfile(['*' fExt], 'Save export file as...', fName);
else
    p = pwd;
    fName = [fInName, fExt];
end
cd(fileparts(mfilename('fullpath')));
if ischar(fName) && textformat
    data = [];
    for ii = 1:size(r, 2) - 1
        data=[data; r(ii).timedata(1:end - 1, :)]; %#ok<AGROW>
    end
    data = [data; r(ii).timedata(end, :)];
    data(:, 2) = data(:, 2) / 1000;  % export in km
    data(:, 3) = data(:, 3) * 3.6;   % convert m/s to kph
    data(:, 26) = data(:, 26) * 3.6; % convert m/s to kph
    data(:, 5) = -data(:, 5);        % change sign of friction brake
    data(:, 28) = -data(:, 28);      % change sign of friction brake
    data(:, 8) = data(:, 8) / 1000;  % export brake resistor power in kW
    if ~isempty(logfileref)
        fprintf(logfileref, [char(datetime, 'HH:mm:ss'),': Route exported as: ', fName, newline]);
    end
    if echo
        hwb = waitbar(0, ['Exporting file: ', strrep(fName, '_', '\_')]);
    else
        hwb = 0;
    end
    xFile = fullfile(p, fName);
    f1 = fopen(xFile, 'w');
    if f1 > 0
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
        fprintf(f1, 'Time [s]\tDist [km]\tSpeed [km/h]\tEffort [kN]\tFriction [kN]\tULine\tILine\tBR power [kW]\tBR temp [C]\tIInverter [A]\tSlip'); % 1-11
        % 12 Modulation index
        % 13 ULink
        % 14 ESPow - Energy Saver Power
        % 15 Main Transformer prim winding RMS
        % 16 Main Transformer 2nd winding RMS
        % 17 Line filter RMS
        % 18 MCM phase current RMS
        % 19 LCM phase current RMS
        % 20 Traction motor RMS
        % 21 Power Reduction factor
        % 22 Effort Reduction factor
        % 23 Running Resistance
        % 24 Traction loss
        fprintf(f1, '\tMod idx\tUDC [V]\tES_Pow[kW]\tMTPrWi\tMTSeWi\tILFiltRMS\tMCM_IPH_RMS\tLCM_IPH_RMS\tTM_RMS\tThPower_red\tThEffort_red\tTrain Resistance\tProp Loss [kW]'); % 12-24
        % 25 Slope force
        % 26 Average speed
        % 27 Average TE/DBE
        % 28 Average FricBE
        % 29 Average line current
        % 30 Tunnel (effective C-factor of tunnel)
        % 31 Notch
        % 32 Curve resistance
        % 33 Tunnel extra
        fprintf(f1, ['\tF_Slope [kN]\tAverage speed [km/h]\tAverage TE/DBE [kN]\tAverage FricBE [kN]\tAverage IL [A]\tTunnel\tNotch\tF_Curve [kN]\tF_Tunnel [kN]', newline]);
        
        for ii = 1:size(data, 1)
            if toc > 0.1
                if echo
                    if ishandle(hwb)
                        waitbar(ii / size(data, 1), hwb);
                    else
                        break; % Break the loop
                    end
                end
                tic;
            end
            if (mod(ii - 1, dsmpl) == 0) || (ii == size(data, 1))
                fprintf(f1, '%1.3f\t%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f', data(ii, 1:24));
                fprintf(f1, '\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.0f\t%1.3f\t%1.3f\r\n', data(ii, 25:33));
            end
        end
        fclose(f1);
    else
        dlg = errordlg(['Can''''t write to file name: ', fName], 'Write to file error');
        uiwait(dlg);
        xFile = '';
    end
    if echo && ishandle(hwb)
        close(hwb);
    end
elseif ischar(fName)
    if ~isempty(logfileref)
        fprintf(logfileref, [char(datetime, 'HH:mm:ss'), ': Route exported as: ', fName, newline]);
    end
    data = [];
    for ii = 1:size(r, 2) - 1
        data = [data; r(ii).timedata(1:end - 1, :)]; %#ok<AGROW>
    end
    data = [data; r(ii).timedata(end, :)];
    xFile = fullfile(p, fName);
    if dsmpl > 1
        data2 = [];
        for ii = 1:size(data, 1)
            if (mod(ii - 1, dsmpl) == 0) || (ii == size(data, 1))
                data2(end + 1, :) = data(ii, :); %#ok<AGROW>
            end
        end
        data = data2;
    end
    save(xFile, 'data');
end
cd(oldwd)
return % ShSim_ExportData
