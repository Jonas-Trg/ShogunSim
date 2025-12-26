function ShSim_Plot_Vehicle(s, data_main, imported_data, sim_config)
% ShSim_Plot_Vehicle
% 0.1.1.3
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-09  JR      New coding design (spaces) plus minor corrections
% 0.1.1.2   2025-12-19  JR      field 'vehicle_file' contains full path and file name of the vehicle file
%
if exist(fullfile(data_main.vehicle_file), 'file')
    if isempty(s)
        fprintf('\nFigures to plot:\n');
        fprintf('a. Fig 200 - Effort vs Speed\n');
        fprintf('b. Fig 201 - Acceleration vs Speed\n');
        fprintf('c. Fig 202 - Power vs Speed\n');
        fprintf('d. Fig 203 - Motor torque vs Speed\n');
        fprintf('e. Fig 204 - System efficiency vs Speed\n');
        fprintf('f. Fig 205 - System losses vs. torque/speed\n');
        fprintf('g. Fig 206 - System losses vs. power/speed\n');
        s=input('Select figure[s] to be plotted: ','s');
    end
    if ~isempty(s)
        vehicle = ShSim_Read_Vehicle(data_main.vehicle_file, 0);
        if exist(data_main.vehicle2_file, 'file')
            vehicle_alt = ShSim_Read_Vehicle(data_main.vehicle2_file, 1);
        else
            vehicle_alt = [];
        end
        s = lower(s);
        figs = zeros(1, 30);
        figs(s(s >= 'a' & s <= 'z') - 'a' + 1) = 1; % decide which figure based on command "s"
        do_plot_vehicle(vehicle,vehicle_alt, imported_data, sim_config, figs);
    end
else
    dlg = warndlg('Vehicle not found!');
    uiwait(dlg);
end
return % ShSim_Plot_Vehicle

function do_plot_vehicle(vehicle, vehicle_alt, imported_data, sim_config, figs)
data = execute_vehicle(vehicle, sim_config); % create sample data
if ~isempty(vehicle_alt)
    data2 = execute_vehicle(vehicle_alt, sim_config);
end
if figs(1)
    stdfig([6, 10] * 80, 200, 'Vehicle effort vs Speed');
    hold on; grid on;
    plot(data(:, 1, 1) * 3.6, data(:, 2, 1), 'b', 'LineWidth', 2, 'DisplayName', 'Tractive effort');
    plot(data(:, 1, 2) * 3.6, data(:, 2, 2), 'g', 'LineWidth', 2, 'DisplayName', 'Dynamic brake effort');
    plot(data(:, 1, 2) * 3.6, data(:, 3, 2), 'r', 'LineWidth', 2, 'DisplayName', 'Friction brake effort');
    plot(data(:, 1, 2) * 3.6, data(:, 11, 1), '--k', 'DisplayName', 'Resistance - Open');
    plot(data(:, 1, 2) * 3.6, data(:, 12, 1), ':k', 'DisplayName', 'Resistance - Tunnel');
    if ~isempty(vehicle_alt)
        plot(data2(:, 1, 1) * 3.6, data2(:, 2, 1), 'c', 'DisplayName', 'Alt: Tractive effort');
        plot(data2(:, 1, 2) * 3.6, data2(:, 2, 2), 'y', 'DisplayName', 'Alt: Dynamic brake effort');
        plot(data2(:, 1, 2) * 3.6, data2(:, 3, 2), 'm', 'DisplayName', 'Alt: Friction brake effort');
        plot(data2(:, 1, 2) * 3.6, data2(:, 11, 1), '--r', 'DisplayName', 'Resistance - Open');
        plot(data2(:, 1, 2) * 3.6, data2(:, 12, 1), ':r', 'DisplayName', 'Resistance - Tunnel');
    end
    if ~isempty(imported_data)
        plot(imported_data(:, 3) * 3.6, imported_data(:,4), 'Color', [1 0 1], 'Linewidth', 1, 'DisplayName', 'Imported data')
    end
    if ~isempty(vehicle_alt)
        xlim([0, max([vehicle.vMax, vehicle_alt.vMax]) * 3.6 * 1.05])
    else
        xlim([0, vehicle.vMax * 3.6 * 1.05])
    end
    ylim([min(data(:, 2, 2)) max(data(:, 2, 1))] * 1.05)
    title('Effort vs Speed');
    ylabel('Effort [kN]');
    xlabel('Speed [km/h]');
    legend('Show', 'Location', 'best')
end
if figs(2)
    stdfig([6, 10] * 80, 201, 'Vehicle acceleration vs Speed');
    hold on; grid on;
    acc = (data(:, 2, 1) - data(:, 11, 1)) / vehicle.totMass(vehicle.LoadCondition);
    plot(data(:, 1, 1) * 3.6, acc, 'b', 'LineWidth', 2, 'DisplayName', 'Open');
    acc = (data(:, 2, 1) - data(:, 12, 1)) / vehicle.totMass(vehicle.LoadCondition);
    plot(data(:, 1, 1) * 3.6, acc, 'g', 'LineWidth', 2, 'DisplayName', 'Tunnel');
    ret = (data(:, 2, 2) + data(:, 3, 2) - data(:, 12, 2)) / vehicle.totMass(vehicle.LoadCondition);
    plot(data(:, 1, 2) * 3.6, ret, 'r', 'LineWidth', 2, 'DisplayName', 'Braking');
    if ~isempty(vehicle_alt)
        acc = (data2(:, 2, 1) - data2(:, 11, 1)) / vehicle_alt.totMass;
        plot(data2(:, 1, 1) * 3.6, acc, 'c', 'LineWidth', 2,'DisplayName', 'Alt: Open');
        acc = (data2(:, 2, 1) - data2(:, 12, 1)) / vehicle_alt.totMass;
        plot(data2(:, 1, 1) * 3.6, acc, 'y', 'LineWidth', 2, 'DisplayName', 'Alt: Tunnel');
        ret = (data2(:, 2, 2) + data2(:, 3, 2) - data2(:, 12, 2)) / vehicle_alt.totMass;
        plot(data2(:, 1, 2) * 3.6, ret, 'm', 'LineWidth', 2, 'DisplayName', 'Alt: Braking');
    end
    if ~isempty(imported_data)
        acc = [0; diff(imported_data(:, 3)) ./ diff(imported_data(:, 1))];
        plot(imported_data(:, 3) * 3.6, acc, 'Color', [1 0 1], 'Linewidth', 1, 'DisplayName', 'Imported data')
    end
    if ~isempty(vehicle_alt)
        xlim([0, max([vehicle.vMax, vehicle_alt.vMax]) * 3.6 * 1.05])
    else
        xlim([0, vehicle.vMax * 3.6 * 1.05])
    end
    ylim([min(ret), max(acc)] * 1.1);
    title('Vehicle acceleration vs Speed');
    ylabel('Acceleration [m/s^2]');
    xlabel('Speed [km/h]');
    legend('Show', 'Location', 'best')
end
if figs(3)
    stdfig([6, 10] * 80, 202, 'Power vs Speed');
    hold on; grid on;
    pow = data(:, 5, 1) .* data(:, 6, 1) / 1000;
    reg = data(:, 5, 2) .* data(:, 6, 2) / 1000;
    plot(data(:, 1, 1) * 3.6, pow, 'b', 'LineWidth', 2, 'DisplayName', 'Line power');
    h = plot(data(:, 1, 1) * 3.6, reg, 'b', 'LineWidth', 2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    pow = data(:, 1, 1) .* data(:, 2, 1);
    reg = data(:, 1, 2) .* data(:, 2, 2);
    plot(data(:, 1, 1) * 3.6, pow, 'g', 'LineWidth', 2, 'DisplayName', 'Rail power');
    h = plot(data(:, 1, 1) * 3.6, reg, 'g', 'LineWidth', 2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    if ~isempty(vehicle_alt)
        pow = data2(:, 5, 1) .* data2(:, 6, 1) / 1000;
        reg = data2(:, 5, 2) .* data2(:, 6, 2) / 1000;
        plot(data2(:, 1, 1) * 3.6, pow, 'c', 'LineWidth', 2, 'DisplayName', 'Line power');
        h = plot(data2(:, 1, 1) * 3.6, reg, 'c', 'LineWidth', 2);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        pow = data2(:, 1, 1) .* data2(:, 2, 1);
        reg = data2(:, 1, 2) .* data2(:, 2, 2);
        plot(data2(:, 1, 1) * 3.6, pow, 'y', 'LineWidth', 2, 'DisplayName', 'Rail power');
        h = plot(data2(:, 1, 1) * 3.6, reg, 'y', 'LineWidth', 2);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    if ~isempty(imported_data)
        plot(imported_data(:, 3) * 3.6, imported_data(:, 6) .* imported_data(:, 7) / 1000, 'm', 'DisplayName', 'Line power - Imported data')
        plot(imported_data(:, 3) * 3.6, imported_data(:, 3) .* imported_data(:, 4), 'c', 'DisplayName', 'Rail power - Imported data')
    elseif ~isempty(vehicle_alt)
        pow = data(:, 5, 1) .* data(:, 6, 1) / 1000 - data(:, 1, 1) .* data(:, 2, 1) - vehicle.AuxPower;
        reg = data(:, 5, 2) .* data(:, 6, 2) / 1000 - data(:, 1, 2) .* data(:, 2, 2) - vehicle.AuxPower;
        plot(data(:, 1, 1) * 3.6, pow, '--r', 'LineWidth', 1, 'DisplayName', 'System losses');
        h = plot(data(:, 1, 1) * 3.6, reg, '--r', 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    if ~isempty(vehicle_alt)
        xlim([0, max([vehicle.vMax, vehicle_alt.vMax]) * 3.6 * 1.05]);
    else
        xlim([0, vehicle.vMax * 3.6 * 1.05]);
    end
    title('Line/Rail vs Speed');
    ylabel('Power [kW]');
    xlabel('Speed [km/h]');
    legend('Show', 'Location', 'best')
end
if figs(4)
	stdfig([6, 10] * 80, 203, 'Motor torque vs Speed');
    hold on; grid on;
    limitation.thermEffort = 1;
    limitation.thermPower = 1;
    TE(1000, 1) = 0;
    DBE(1000, 1) = 0;
    for ii = 1:1001
        TE(ii, 1) = ShSim_maxTE(data(ii, 2, 1), vehicle,data(ii, 1, 1), limitation);
        DBE(ii, 1) = ShSim_maxTE(data(ii, 2, 2), vehicle,data(ii, 1, 2), limitation);
    end
    tq = TE * vehicle.wheelsize{vehicle.ws_selection} * 1000 / 2 / vehicle.gearRatio / vehicle.nAxDrive;
    rpm = data(:, 1, 1) * vehicle.gearRatio / (vehicle.wheelsize{vehicle.ws_selection} * pi) * 60;
    plot(rpm, tq, 'b', 'DisplayName', 'Motor', 'LineWidth', 2);
    tq = DBE * vehicle.wheelsize{vehicle.ws_selection} * 1000 / 2 / vehicle.gearRatio / vehicle.nAxDrive;
    rpm = data(:, 1, 2) * vehicle.gearRatio / (vehicle.wheelsize{vehicle.ws_selection} * pi) * 60;
    plot(rpm, tq, 'r', 'DisplayName', 'Generator', 'LineWidth', 2);
    switch vehicle.motorMethod
        case 2
            plot(vehicle.motorMap.SpeedMT(:, end), vehicle.motorMap.TorqueMT(:, end), '--k', 'Displayname', 'Motor map');
            h = plot(vehicle.motorMap.SpeedMP(:, end), vehicle.motorMap.PowerMP(:, end) ./ (vehicle.motorMap.SpeedMP(:, end) * 2 * pi / 60), '--k', 'Displayname', 'Motor map');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h = plot(vehicle.motorMap.SpeedGT(:, end), vehicle.motorMap.TorqueGT(:, end), '--k', 'Displayname', 'Motor map');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h = plot(vehicle.motorMap.SpeedGP(:, end), vehicle.motorMap.PowerGP(:, end) ./ (vehicle.motorMap.SpeedGP(:, end) * 2 * pi / 60), '--k');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h = plot([vehicle.motorMap.SpeedMP(end, end), vehicle.motorMap.SpeedGP(end, end)],[vehicle.motorMap.PowerMP(end, end), vehicle.motorMap.PowerGP(end, end)] ./ (vehicle.motorMap.SpeedGP(end, end) * 2 * pi / 60), '--k');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            plot(vehicle.motorMap.rpmBaseMotor, vehicle.motorMap.TorqueMT(end, end), 'sk', 'Displayname', 'Rated point')
        case 3
            plot(vehicle.motor.N_RPM0, vehicle.motor.Torque0, 'sk', 'Displayname', 'Rated point')
        case 4
            plot(vehicle.motor.N_RPM0, vehicle.motor.Torque0, 'sk', 'Displayname', 'Rated point')
    end
    if ~isempty(imported_data)
        tq = imported_data(:,4) * vehicle.wheelsize{vehicle.ws_selection} * 1000 / 2 / vehicle.gearRatio / vehicle.nAxDrive;
        rpm = imported_data(:,3) * vehicle.gearRatio / (vehicle.wheelsize{vehicle.ws_selection} * pi) * 60;
        b = find(tq >= 0);
        plot(rpm(b), tq(b), '.c', 'DisplayName', 'Imported data');
        b = find(tq < 0);
        plot(rpm(b), tq(b), '.m', 'DisplayName', 'Imported data');
    end
    title('Motor Torque vs Motor Speed');
    ylabel('Torque [Nm]');
    xlabel('Speed [rpm]');
    legend('Show', 'Location', 'best')
end
if figs(5) %system efficiency
    stdfig([6, 10] * 80, 204, 'Efficiency vs Speed');
    hold on; grid on;
    LP = data(:, 5, 1) .* data(:, 6, 1) / 1000 - vehicle.AuxPower / vehicle.APS_efficiency;
    RP = data(:, 1, 1) .* data(:, 2, 1);
    eff = RP ./ LP * 100;
    plot(data(:, 1, 1) * 3.6, eff, 'r', 'DisplayName', 'Traction');
    LP = data(:, 5, 2) .* data(:, 6, 2) / 1000 - vehicle.AuxPower / vehicle.APS_efficiency;
    RP = data(:, 1, 2) .* data(:, 2, 2);
    eff = LP ./ RP * 100;
    plot(data(:, 1, 2) * 3.6, eff, 'g', 'DisplayName', 'Braking');
    plot([0, vehicle.vMax*3.6],[1, 1] * vehicle.APS_efficiency * 100, 'b', 'DisplayName', 'APS');
%     if ~isempty(imported_data)
%         b=find(imported_data(:,4)>=0);
%         eff=imported_data(b,3).*imported_data(b,4)./(imported_data(b,6).*imported_data(b,7)/1000 - vehicle.AuxPower/vehicle.APS_efficiency)*100;
%         eff=imported_data(b,3).*imported_data(b,4)./(imported_data(b,6).*imported_data(b,7)/1000)*100;
%         plot(imported_data(b,3)*3.6,eff,'.m','DisplayName','Imported - traction')
%         b=find(imported_data(:,4)<0);
%         eff=(imported_data(b,6).*imported_data(b,7)/1000 - vehicle.AuxPower/vehicle.APS_efficiency)./(imported_data(b,3).*imported_data(b,4))*100;
%         plot(imported_data(b,3)*3.6,eff,'.c','DisplayName','Imported - braking')
%     end
    ylim([0, 100]);
    xlim([0, vehicle.vMax * 3.6 * 1.05])
    title('System efficiency vs Speed');
    ylabel('Efficiency [%]');
    xlabel('Speed [km/h]');
    legend('Show', 'Location', 'best')
end
if figs(6) || figs(7)
    [X, Y, Z] = create_surf(vehicle,sim_config);
    if figs(6)
    	stdfig([6, 8] * 80, 205, 'System losses vs Speed/Torque');
        hold on; grid on;
        if size(imported_data)
            h = surf(X * 3.6, Y, Z, 'LineStyle', 'none', 'FaceAlpha', 0.75, 'DisplayName', char(vehicle.class)); %#ok<NASGU>
            loss = imported_data(:, 6) .* imported_data(:, 7) / 1000 - imported_data(:, 3) .* imported_data(:, 4) - vehicle.AuxPower / vehicle.APS_efficiency;
            plot3(imported_data(:, 3) * 3.6, imported_data(:, 4), loss, '.r', 'DisplayName', 'Imported Data');
            legend('Show')
        else
            surf(X * 3.6, Y, Z);
        end
        rotate3d on
        title('System losses vs Speed/Torque');
        xlabel('Speed [km/h]');
        ylabel('Effort [kN]');
        zlabel('Loss [kW]');
    end
    if figs(7)
        stdfig([6, 8] * 80, 206, 'System losses vs Speed/Power');
        hold on; grid on;
        if size(imported_data)
            h = surf(X * 3.6, Y .* X, Z, 'LineStyle', 'none', 'FaceAlpha', 0.75, 'DisplayName', char(vehicle.class)); %#ok<NASGU>
            loss = imported_data(:, 6) .* imported_data(:, 7) / 1000 - imported_data(:, 3) .* imported_data(:, 4) - vehicle.AuxPower / vehicle.APS_efficiency;
            plot3(imported_data(:, 3) * 3.6, imported_data(:, 4) .* imported_data(:, 3), loss, '.r', 'DisplayName', 'Imported Data');
            legend('Show')
        else
            surf(X * 3.6, Y .* X, Z);
        end
        rotate3d on
        title('System losses vs Speed/Power');
        xlabel('Speed [km/h]');
        ylabel('Power [kW]');
        zlabel('Loss [kW]');
    end
end
return % do_plot_vehicle

function [X, Y, Z] = create_surf(vehicle,sim_config)
vBaseT = vehicle.maxPwrTEbase / vehicle.maxTE;
vBaseB = vehicle.maxPwrDBEbase / vehicle.maxDBE;
pp = round(vBaseT / vehicle.vMax * 100);
bp = round(vBaseB / vehicle.vMax * 100);
X(101, 41) = 0;
Y(101, 41) = 0;
Z(101, 41) = 0;
vehicle.ILBrakeMax = inf;
for ii = 0:100 % speed index
    if ii <= pp
        spdT = ii / pp * vBaseT;
    else
        spdT = (ii - pp) / (100 - pp) * (vehicle.vMax - vBaseT) + vBaseT;
    end
    if ii <= bp
        spdB = ii / bp * vBaseB;
    else
        spdB = (ii - bp) / (100 - bp) * (vehicle.vMax - vBaseB) + vBaseB;
    end
    X(ii + 1, 1:20) = spdB;
    X(ii + 1, 21:41) = spdT;
    for jj = -20:20 % torque
        if jj >= 0
            if ii <= pp
                F_target = jj * vehicle.maxTE / 20;
            else
                F_target = jj * vehicle.maxPwrTEbase / spdT / 20;
            end
            [F, UL, UDC, IL, BrP, iPh, IS, slip, mode, ESP, loss] = ShSim_Efforts(F_target, vehicle, spdT, [], sim_config); %#ok<ASGLU>
        else
            if ii <= bp
                F_target = jj * vehicle.maxDBE / 20;
            else
                F_target = jj * vehicle.maxPwrDBEbase / spdB / 20;
            end
            [F, UL, UDC, IL, BrP, iPh, IS, slip, mode, ESP, loss] = ShSim_Efforts(F_target, vehicle, spdB, [], sim_config); %#ok<ASGLU>
        end
        Y(ii + 1, jj + 21) = F;
        Z(ii + 1, jj + 21) = loss;
    end
end
return % create_surf

function data=execute_vehicle(vehicle,sim_config)
data(1001, 20, 2) = 0; % allocate memory, speed,effort
davies_a = vehicle.davies_a;
davies_b = vehicle.davies_b;
davies_c1 = vehicle.davies_c1;
davies_c2 = vehicle.davies_c2;
vehicle.dDBE_jerkLimit = inf;
vehicle.jerkrateB = 0;
for ii = 1:1001
    if ii < 1001
        speed = ii / 1000 * vehicle.vMax;
        [F, UL, UDC, IL, BrP, iPh, IS, slip, mode, ESP,loss] = ShSim_Efforts(vehicle.maxTE, vehicle, speed, [], sim_config); %#ok<ASGLU>
    else
        speed = vehicle.vMax;
        [F, UL, UDC, IL, BrP, iPh, IS, slip, mode, ESP,loss] = ShSim_Efforts(0, vehicle, speed, [], sim_config); %#ok<ASGLU>
    end
    data(ii,1,1) = speed;
    data(ii,2,1) = F;
    data(ii,3,1) = 0;
    data(ii,4,1) = UDC;
    data(ii,5,1) = UL;
    data(ii,6,1) = IL;
    data(ii,7,1) = iPh;
    data(ii,8,1) = slip;
    data(ii,9,1) = mode;
    data(ii,10,1) = loss;
    data(ii,11,1)=(davies_a + davies_b * speed * 3.6 + davies_c1 * (speed * 3.6)^2);
    data(ii,12,1)=(davies_a + davies_b * speed * 3.6 + davies_c2 * (speed * 3.6)^2);
    
    if ii < 1001
        [F_target, d] = ShSim_TargetBE(-vehicle.maxDBE, 0, speed, inf, vehicle, 0, 0, eps); %#ok<ASGLU>
    else
        F_target = 0;
    end
    [F, UL, UDC, IL, BrP, iPh, IS, slip, mode, ESP, loss] = ShSim_Efforts(F_target, vehicle, speed, [], sim_config); %#ok<ASGLU>
    data(ii,1,2) = speed;
    data(ii,2,2) = F;
    data(ii,3,2) = F_target - F;
    data(ii,4,2) = UDC;
    data(ii,5,2) = UL;
    data(ii,6,2) = IL;
    data(ii,7,2) = iPh;
    data(ii,8,2) = slip;
    data(ii,9,2) = mode;
    data(ii,10,2) = loss;
    data(ii,11,2) = data(ii,11,1);
    data(ii,12,2) = data(ii,12,1);
end
return % execute_vehicle

