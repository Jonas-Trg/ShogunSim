function [sections, route, profile, curves] = ShSim_Read_Route(fileName, vehicle, sim_config)
% function ShSim_Read_Route
% 0.1.1.1
%
% Description
%   This function reads the excel file (or its image file stored as mat)
%   and extracts the raw data into the route definition used for the
%   calculations.
%
% Input data
%   pName - char
%       path to the route data excel file
%   fName - char
%       file name of the route data excel file
%   vehicle - struct of (least contain):
%       length      length of unit
%       nUnits      number of units
% Output data
%   sections
%
% Check if mat-file exists
%
% 0.1.1.1   2025-12-11  JR      New coding design standard
%
if exist(fileName,'file')
    load(fileName, 'route_data');
    [sections, route, profile, curves] = read_route(route_data, vehicle, sim_config);
end
return % ShSim_Read_Route

function [sections, route, track_profile, curves] = read_route(route_data, theVehicle, sim_config)
%This to ensure the function doesn't crash at 'return'
sections = [];
route = [];
track_profile = [];
curves = [];
if ~isempty(theVehicle)
    TrainLength = theVehicle.unitLength * theVehicle.nUnits / 1000;
else
    TrainLength = 0;
end

if isfield(route_data, 'stations')
    stations = route_data.stations;
    for i1 = 1:size(stations, 1)
        stations{i1, 1} = stations{i1, 1};
        stations{i1, 3} = stations{i1, 3};
        if isempty(stations{i1, 4})
            stations{i1, 4} = 0.5;
        else
            stations{i1, 4} = stations{i1, 4};
        end
    end
else
    dlg = errordlg(['The trackdata file ''', route_data.FileName, ''' has no field ''Stations'''], mfilename);
    uiwait(dlg);
    return
end

switch route_data.ATOSpeedUnitSel
    case 1 %
        spd_scale = 3.6;
    case 2
        spd_scale = 2.2369;
    otherwise
        spd_scale = 1;
end

if isfield(route_data, 'speedLimits')
    n = size(route_data.speedLimits, 1);
    speeds(n, 2) = 0;
    for i1 = 1:n
        speeds(i1, 1) = route_data.speedLimits{i1, 1};
        speeds(i1, 2) = route_data.speedLimits{i1, 2} / spd_scale;
    end
else
    dlg = errordlg(['The trackdata file ''', route_data.FileName, ''' has no field ''Speedlimits'''], mfilename);
    uiwait(dlg);
    return
end

if isfield(route_data, 'gradients')
    n = size(route_data.gradients, 1);
    gradients(n, 2) = 0;
    for i1 = 1:n
        gradients(i1, 1) = route_data.gradients{i1, 1};
        gradients(i1, 2) = route_data.gradients{i1, 2};
    end
else
    gradients = [0, 0];
end

if isfield(route_data, 'curves')
    n = size(route_data.curves, 1);
    curves(n, 3) = 0;
    curveNames = route_data.curves(:, 5);
    for i1 = 1:n
        curves(i1, 1) = route_data.curves{i1, 1};
        if strcmp(route_data.curves{i1, 2}, '==')
            curves(i1, 2) = 0;
        else
            % Convert to curvature 1000 / radius
            curves(i1, 2) = 1000 / route_data.curves{i1, 2};
        end
    end
    
else
    curves = [0, 0, 0];
end

if isfield(route_data, 'tunnels') && size(route_data.tunnels, 1) > 0
    n = size(route_data.tunnels, 1);
    tunnels(n, 3) = 0;
    tunnelNames = route_data.tunnels(:, 3);
    for i1 = 1:n
        tunnels(n, 1) = route_data.tunnels{i1, 1};
        tunnels(n, 2) = route_data.tunnels{i1, 2};
        tunnels(n, 3) = route_data.tunnels{i1, 3};
    end
else
    tunnels = [];
    tunnelNames = {};
end
if isfield(route_data,'parameters')
    pars = route_data.parameters;
    parsT = route_data.parameterNames(2:end);
else
    pars = [];
end

nStations = size(stations, 1);
original_track(nStations, 8) = 0;
%matrix containing 9 columns + number of parameters
%   dist            section start in km 
%   speed limit     in km/h
%   slope           in per mille
%   curvature       in rad/km
%   curvature'      
%
parameter.name = {};
parameter = [];
departureTime(nStations) = 0;
dwellTime(nStations) = 0;
for i1 = 1:nStations
    stations{i1, 1} = stations{i1, 1} + stations{i1, 4} * TrainLength;
    original_track(i1, 1) = stations{i1, 1};
    dwellTime(i1) = stations{i1, 3};
    T_dur = seconds(duration(stations{i1, 5}));
    if ~isnan(T_dur)
        departureTime(i1) = T_dur;
    end
end

%     [pars,vStrs,vIdx]=check4vehicle(pars,parsT);


dist = -inf;
for i1 = 1:size(speeds, 1)
    if speeds(i1, 1) > dist
        dist = speeds(i1, 1);
    else
        dlg = errordlg('Error in speed distances','Reading Track data');
        uiwait(dlg);
        return
    end
end

i1 = 1;
while i1 < size(speeds, 1)
    i1 = i1 + 1;
    if speeds(i1, 2) >= speeds(i1 - 1, 2) % if new speed is higher...
        speeds(i1, 1) = speeds(i1, 1) + TrainLength; % train length will be added
        if i1 < size(speeds, 1)
            if speeds(i1, 1) > speeds(i1 + 1, 1) && speeds(i1, 2) > speeds(i1 + 1, 2)
                if speeds(i1 + 1, 2) > speeds(i1 - 1, 2)
                    speeds(i1 + 1, 1) = speeds(i1, 1) - TrainLength;
                end
                speeds = speeds([1:i1 - 1, i1 + 1:end], :);
                i1 = i1 - 1;
            end
        end
    elseif speeds(i1, 2) < speeds(i1 - 1, 2)
        if speeds(i1, 1) < speeds(i1 - 1, 1)
            speeds(i1 - 1, 2) = speeds(i1, 2);
            speeds = [speeds(1:i1 - 1, :); speeds(i1 + 1:end, :)];
            i1 = i1 - 1;
        end
    end
end
for i1 = 1:size(speeds, 1)
    b = find(original_track(:, 1) >= speeds(i1, 1), 1, 'first');
    if ~isempty(b)
        if original_track(b, 1) > speeds(i1, 1) && b > 1
            original_track = [original_track(1:b - 1, :); [speeds(i1, 1:2), original_track(b - 1, 3:end)]; original_track(b:end, :)];
        else
            original_track(b, 2) = speeds(i1, 2);
        end
        original_track(b + 1:end, 2) = speeds(i1, 2);
    end
end


% orgtrack: dist, speed, slope, curve, curve transition, cant, cant transition, tunnel coeff
track_profile = gradients(:, 1:2); % only interested in the first two columns

dist = -inf;
for i1 = 1:size(track_profile, 1)
    if track_profile(i1, 1) > dist
        dist = track_profile(i1, 1);
    else
        dlg = errordlg('Error in profile distances', 'Reading Track data');
        uiwait(dlg);
        return
    end
end

track_profile(:, 1) = round(track_profile(:, 1) * 1000); % rounded to meters
nGradients = size(gradients, 1);
for i1 = 1:nGradients
    idx = find(original_track(:,1) >= gradients(i1, 1), 1, 'first');
    if ~isempty(idx)
        if original_track(idx, 1) > gradients(i1, 1) && idx > 1
            original_track = [original_track(1:idx-1,:); [gradients(i1,1), original_track(idx - 1, 2), gradients(i1, 2), original_track(idx - 1, 4:end)]; original_track(idx:end, :)];
        else
            original_track(idx, 3) = gradients(i1, 2);
        end
        original_track(idx + 1:end, 3) = gradients(i1, 2);
    end
end

% Curve radius is translated to "curvature" eq to rad/km described with a
% complex number curvature/transition (rad/km & rad/km2)
dist = -inf;
for i1 = 1:size(curves, 1)
    if curves(i1, 1) > dist
        dist = curves(i1, 1);
    else
        dlg = errordlg(['Error in curve distances. Row ' num2str(i1)], 'Reading Track data');
        uiwait(dlg);
        return
    end
end

nCurves = size(curves, 1);
b=find(isnan(curves(:, 2)));
for i1=b' % Create transitions (replace nan's)
    curves(i1, 2) = curves(i1 - 1, 2) + 1i * (curves(i1 + 1, 2) - curves(i1 - 1, 2)) / (curves(i1 + 1, 1) - curves(i1, 1));
    curves(i1, 3) = curves(i1 - 1, 3) + 1i * (curves(i1 + 1, 3) - curves(i1 - 1, 3)) / (curves(i1 + 1, 1) - curves(i1, 1));
end
lastIdx = 0;
for i1 = 1:nCurves
    idx = find(original_track(:, 1) >= curves(i1, 1), 1, 'first');
    if ~isempty(idx)
        if original_track(idx, 1) > curves(i1, 1) && idx > 1
%             orgtrack = [orgtrack(1:idx-1,:);[curves(ii,1) orgtrack(idx-1,2:3) real(curves(ii,2)) imag(curves(ii,2)) real(curves(ii,3)) imag(curves(ii,3)) orgtrack(idx-1,8:end)];orgtrack(idx:end,:)];
            original_track = [original_track(1:idx-1, :); original_track(idx-1:end, :)]; % duplicate line idx
            original_track(idx, 1) = curves(i1, 1);
        end
        if lastIdx > 0
            for i2 = lastIdx + 1:idx - 1
                original_track(i2, 4) = original_track(lastIdx, 4) + original_track(lastIdx, 5) * (original_track(i2, 1) - original_track(lastIdx, 1));
                original_track(i2, 5) = original_track(lastIdx, 5);
                original_track(i2, 6) = original_track(lastIdx, 6) + original_track(lastIdx, 7) * (original_track(i2, 1) - original_track(lastIdx, 1));
                original_track(i2, 7) = original_track(lastIdx, 7);
            end
        end
        original_track(idx, 4:7) = [real(curves(i1, 2)), imag(curves(i1, 2)), real(curves(i1, 3)), imag(curves(i1, 3))];
        original_track(idx + 1:end, 4) = real(curves(i1, 2));
        original_track(idx + 1:end, 6) = real(curves(i1, 3));
        lastIdx = idx;
    end
end

dist = -inf;
for i1 = 1:size(tunnels, 1)
    if tunnels(i1, 1) - dist > -0.001
        dist = tunnels(i1, 1) + tunnels(i1, 2);
    else
        dlg = errordlg('Error in tunnel distances','Reading Track data');
        uiwait(dlg);
        return
    end
end

nTunnels = size(tunnels, 1);
for i1 = 1:nTunnels
    idx = find(original_track(:, 1) >= tunnels(i1, 1), 1, 'first');
    if ~isempty(idx)
        if original_track(idx, 1) > tunnels(i1, 1) && idx > 1
            original_track = [original_track(1:idx - 1, :); [tunnels(i1, 1), original_track(idx - 1, 2:7), tunnels(i1, 3), original_track(idx - 1, 9:end)]; original_track(idx:end, :)];
            original_track(idx, 4) = original_track(idx - 1, 4) + original_track(idx - 1, 5) * (original_track(idx, 1) - original_track(idx - 1, 1));
            original_track(idx, 6) = original_track(idx - 1, 6) + original_track(idx - 1, 7) * (original_track(idx, 1) - original_track(idx - 1, 1));
        else
            original_track(idx, 8) = tunnels(i1, 3);
        end
        original_track(idx + 1:end, 8) = tunnels(i1, 3);
        tunnelEnd = tunnels(i1, 1) + tunnels(i1, 2);
        idx = find(original_track(:, 1) >= tunnelEnd, 1, 'first');
        if ~isempty(idx) && idx > 1
            original_track = [original_track(1:idx-1, :); [tunnelEnd, original_track(idx - 1, 2:7), 0, original_track(idx - 1, 9:end)]; original_track(idx:end, :)];
            original_track(idx + 1:end, 8) = 0;
        else
            original_track(idx:end, 8) = 0;
        end
    end
end
if ~isempty(pars) % && size(pars,2)==size(parsT,2)
    [pars, vStrs, vIdx] = check4vehicle(pars, parsT);
    nParsT = size(parsT, 2) - 2;
    
    for i1=1:nParsT
        parameter(end + 1).name = parsT(1, 2 + i1); %#ok<AGROW>
        original_track(:, end + 1) = NaN; %#ok<AGROW>
    end
    for i1 = 1:nParsT
        parameter(end - (nParsT - i1)).dist = pars(:, 1)'; %#ok<AGROW>
        parameter(end - (nParsT - i1)).value = pars(:, 2 + i1)'; %#ok<AGROW>
    end
    if ~isempty(vStrs)
        parameter(end + vIdx).vehicles = vStrs;
    end
    nParRows = size(pars, 1);
    pars(:, 1) = pars(:, 1) + TrainLength * (sim_config.plotRef - 1) / 2;
    for i1 = 1:nParRows
        idx = find(original_track(:, 1) >= pars(i1, 1), 1, 'first');
        if ~isempty(idx)
            if original_track(idx, 1) > pars(i1, 1) && idx > 1
                original_track = [original_track(1:idx - 1, :); [pars(i1, 1), original_track(idx - 1, 2:end - nParsT), pars(i1, 3:2 + nParsT)]; original_track(idx:end, :)];
                original_track(idx, 4) = original_track(idx - 1, 4) + original_track(idx - 1, 5) * (original_track(idx, 1) - original_track(idx - 1, 1));
                original_track(idx, 6) = original_track(idx - 1, 6) + original_track(idx - 1, 7) * (original_track(idx, 1) - original_track(idx - 1, 1));
            else
                original_track(idx, end - nParsT + 1:end) = pars(i1, 3:nParsT + 2);
            end
            for i2 = 1:nParsT
                original_track(idx + 1:end, end - nParsT + i2) = pars(i1, 2 + i2);
            end
        end
    end
end

i1=1;
while i1 < size(stations, 1)
    if (dwellTime(i1 + 1) == 0 || isnan(dwellTime(i1 + 1))) && i1 < size(stations, 1) - 1 % remove this station from the list
        stations = [stations(1:i1, :); stations(i1 + 2:end, :)];
        dwellTime = [dwellTime(1:i1), dwellTime(i1 + 2:end)];
        departureTime = [departureTime(1:i1), departureTime(i1 + 2:end)];
    else
        trackdata(i1).from_station = stations{i1, 2}; %#ok<AGROW>
        trackdata(i1).to_station = stations{i1 + 1, 2}; %#ok<AGROW>
        trackdata(i1).dwellTime = dwellTime(i1 + 1); %#ok<AGROW>
        trackdata(i1).departure = departureTime(i1) * 24 * 3600; %#ok<AGROW>
        trackdata(i1).arrival = departureTime(i1 + 1) * 24 * 3600 - dwellTime(i1 + 1); %#ok<AGROW>
        trackdata(i1).stopPos = [stations{i1, 5}, stations{i1 + 1, 5}]; %#ok<AGROW>
        b = find(original_track(:, 1) == stations{i1, 1}, 1, 'first');
        c = find(original_track(:, 1) == stations{i1 + 1, 1}, 1, 'first');
        if ~isempty(b & c)
            trackdata(i1).from_track = b; %#ok<AGROW>
            trackdata(i1).to_track = c; %#ok<AGROW>
        else
            fprintf('Error...!\n');
        end
        i1 = i1 + 1;
    end
end

if exist('trackdata','var')
    sections = trackdata;
else
    dlg = errordlg('Track data file is non valid');
    uiwait(dlg);
    return
end
% trackdata - matrix of struct with:
%     The matrix of trackdata contains all sections of the route. Each
%     section is described with the following fields
%         from_station - char
%             Station name of departure
%         to_station - char
%             Station name of arrival
%         dwellTime - num
%             dwell time at arrival station
%         departure - num
%             time for departure at "from_station"
%         arrival - num
%             time for arrival at "to_station"
%         track - matrix of sub-sections describing
%             1-dist
%             2-max speed
%             3-slope
route.name = route_data.name;
route.purpose = route_data.purpose;
if route_data.VoltageSel >= 3
    uScale = 1000;
else
    uScale = 1;
end
route.USupply = replaceNAN(route_data.USupply * uScale, 1500);
route.UMotor = replaceNAN(route_data.UL_Motoring * uScale, route.USupply);
route.URegen = replaceNAN(route_data.UL_Regen * uScale, route.USupply);
route.Regen_percent = replaceNAN(route_data.Regen_percent, 1);
route.Regen_current = replaceNAN(route_data.Regen_current, inf);
route.Start_direction = replaceNAN(route_data.Start_direction, 90) * 2 * pi / 360;
route.start_lat = replaceNAN(route_data.start_lat, 0);
route.start_lon = replaceNAN(route_data.start_lon, 0);
route.ambTemp = replaceNAN(route_data.ambTemp, 15);
route.relHumidity = replaceNAN(route_data.relHumidity, 0);
route.headWind = replaceNAN(route_data.headWind, 0);
route.Track_Guage = replaceNAN(route_data.Track_Guage, 1435);

route.ATO_SpdMgn = replaceNAN(route_data.ATO_SpdMgn, 0) / spd_scale;
route.ATO_SpdMax = replaceNAN(route_data.ATO_SpdMax, inf) / spd_scale;
route.ATO_InitBrkSpd = replaceNAN(route_data.ATO_InitBrkSpd, inf) / spd_scale;
route.ATO_TBC = replaceNAN(route_data.ATO_TBC, 100) / 100;
route.ATO_dCoastStart = replaceNAN(route_data.ATO_dCoastStart, inf);
route.ATO_dCoastEnd = replaceNAN(route_data.ATO_dCoastEnd, inf);
route.ATO_vCoastStart = replaceNAN(route_data.ATO_vCoastStart, inf) / spd_scale;
route.ATO_vCoastEnd = replaceNAN(route_data.ATO_vCoastEnd, 0) / spd_scale;
route.ATO_V1 = replaceNAN(route_data.ATO_V1, inf) / spd_scale;
route.ATO_V2 = replaceNAN(route_data.ATO_V2, inf) / spd_scale;
route.ATO_cGrad = replaceNAN(route_data.ATO_cGrad, 0);

route.tunnelNames = tunnelNames;
route.curveNames = curveNames;
route.parameters = parameter;
route.track = original_track;
route.revtrack = reverseTrack(original_track);
return % read_track

function revTrack = reverseTrack(track)
revTrack = flipud(track);
revTrack(:, 1) = track(end, 1) - revTrack(:, 1);
revTrack(1:end - 1, 2) = revTrack(2:end, 2);
revTrack(1:end - 1, 3) = -revTrack(2:end, 3);
% Curves, including transients
revTrack(:, 4) = -revTrack(:, 4);
revTrack(1:end - 1, 5) = revTrack(2:end, 5);
revTrack(:, 6) = revTrack(:, 6);
revTrack(1:end - 1, 7) = -revTrack(2:end, 7);
% tunnels
revTrack(1:end - 1, 8) = revTrack(2:end, 8);
% Parameters
return % reverseTrack

function retVal = replaceNAN(val, replacementValue)
if isempty(val) || isnan(val)
    retVal = replacementValue;
else
    retVal = val;
end
return % replaceNAN

function [pars, vStrs, idx] = check4vehicle(iPars, iParsT)
% Check iParsT (text array from excel) for any 'Vehicle', creates an cell
% text array of the vehicles in 'vStrs' and sets the parameter value to the
% index in this cell text array.
% Input
%    ipars -     numerical array from excel read
%    iparsT      cell text array from excel read
% Output
%    pars        
%    vStrs
%    idx
%
pars=iPars;
vStrs={};n=0;idx=2;
for ii=3:size(iParsT,2)
    if strcmpi(char(iParsT(1,ii)),'vehicle')
        for jj=2:size(iParsT,1)
            if ~isempty(char(iParsT(jj,ii)))
                n=n+1;
                if ~contains(char(iParsT(jj,ii)),'.xls')
                    iParsT(jj,ii)=cellstr([char(iParsT(jj,ii)) '.xlsx']);
                end
                vStrs(n)=iParsT(jj,ii); %#ok<AGROW>
                idx=ii-size(iParsT,2);
                pars(jj-1,ii)=n;
            end
        end
    end
end
return %check4vehicle

function mTrackFile=check4mat(pathname, filename)
p=find(filename=='.',1,'last');
mTrackFile=[filename(1:p-1) '.mat'];
if exist(fullfile(pathname, mTrackFile),'file')
    xlsinfo=dir(fullfile(pathname, filename));
    matinfo=dir(fullfile(pathname, mTrackFile));
    if xlsinfo.datenum>matinfo.datenum
        fprintf(['XL file ''' filename ''' is newer - re-load track data\n']);
        delete(fullfile(pathname, mTrackFile));
    end
else
    fprintf(['No mat-file found, ''' mTrackFile ''' will be generated\n']);
end

return %check4mat
