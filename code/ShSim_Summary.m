function r_output = ShSim_Summary(r_input, timedata, simulation, trackdata, routedata, vehicle, data)
% ShSim_Summary
% 0.1.1.1
%
% Format
%   function r_output=ShSim_Summary(r_input,timedata,simulation,trackdata,vehicle)
%
% Description ShSim_Summary
%   This function creates the summary from the simulation of the section.
%   
% Input data
%   struct r_input
%       fName - file name?
%   matrix timedata
%   struct simulation
%
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-26  JR      New coding design (spaces)
%
%
% Generic data
result.from_station = trackdata.from_station;
result.to_station = trackdata.to_station;
result.vehicle_file = data.vehicle_file;
result.track_file = data.track_file;
result.depart = r_input.depart;
result.arrival = r_input.arrival;

result.runtime = r_input.runtime;
result.distance = routedata.track(trackdata.to_track, 1) - routedata.track(trackdata.from_track, 1);

% Set the order of the result record
result.w_traction   = 0;
result.w_regen      = 0;
result.w_Aux        = 0;
result.w_railPos    = 0;
result.w_railNeg    = 0;
result.w_resistance = 0;
result.w_friction   = 0;

b = find(timedata(2:end, 29) >= 0) + 1; % part of run with positive line current
result.w_traction   = sum(timedata(b, 29) .* timedata(b, 6) * simulation.delta_T / 3.6e6); % ULine x ILine x time kWh
b = find(timedata(2:end, 29) < 0) + 1; % part of run with regenerative line current
result.w_regen      = sum(timedata(b, 29) .* timedata(b,6) * simulation.delta_T / 3.6e6); % ULine x ILine x time kWh
b = find(timedata(2:end, 27) >= 0) + 1; % part of run with positive effort
result.w_railPos    = sum(timedata(b,26) .* timedata(b,27) * simulation.delta_T/3600); %Speed x Effort x time kWh
b=find(timedata(2:end,27)<0)+1;
result.w_railNeg    = abs(sum(timedata(b,26) .* timedata(b,27) * simulation.delta_T/3600)); %Speed x Effort x time kWh
result.w_friction   = sum(timedata(:,3) .* timedata(:,5) * simulation.delta_T/3600); %Speed x Effort x time kWh
result.w_friction2  = sum(timedata(:,26) .* timedata(:,28) * simulation.delta_T/3600); %Speed x Effort x time kWh
result.w_resistance = sum(timedata(:,3) .* timedata(:,23) * simulation.delta_T/3600); %Speed x force x time kWh
result.w_PropLoss   = sum(timedata(:,24)*simulation.delta_T/3.6); %power x time kWh
result.w_Aux        = (timedata(end,1)-timedata(1,1))*vehicle.AuxPower/vehicle.APS_efficiency/3600; %Power kWh
result.i_rms_sqsum  = sum(timedata(:,29).^2); %sqrt(sum(timedata(:,7).^2)/size(timedata,1));
result.i_rms_num    = size(timedata,1);

result.i_rms2_sqsum = sum((timedata(:,7).*timedata(:,3)>0).^2);
b=find(timedata(:,3)>0);
result.i_rms2_sqsum = result.i_rms_sqsum;
result.i_rms2_num   = size(b,1);
b=find(timedata(:,27)+timedata(:,28)==0 & timedata(:,26)>0); %costing
dTime=[0;diff(timedata(:,1))];
result.t_coast=sum(dTime(b));
result.d_coast=sum(dTime(b).*timedata(b,26));

% result.v_average=result.distance*3600/r_input.runtime;
% result.n_average=result.v_average/3.6/vehicle.wheelsize/pi*vehicle.gearRatio*60;
if simulation.ThermReduction
    result.ThermRedEffortSum=sum(timedata(:,22));
    result.ThermRedPowerSum=sum(timedata(:,21));
%     sum(timedata(:,22).*timedata(:,3))/sum(timedata(:,3));
%     result.ThermRedPower=sum(timedata(:,21))/size(timedata(:,21),1);
else
    result.ThermRedEffortSum=size(timedata,1);
    result.ThermRedPowerSum=size(timedata,1);
end
result.nThermRecords=size(timedata,1);
result.ATO_SpdMax=[];
result.ATO_InitBrkSpd=[];
result.ATO_TBC=[];
result.ATO_V1=[];
result.ATO_V2=[];
result.ATO_dCoastStart=[];
result.ATO_dCoastEnd=[];
r_output=result;
return % ShSim_Summary
