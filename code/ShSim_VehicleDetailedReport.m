function ShSim_VehicleDetailedReport(ActXWord, vehicle, data_main, sim_config, level, head, head2)
% ShSim_VehicleDetailedReport
% 0.1.1.1
%
% 0.1.1.1   2025-12-22  JR      New coding design (spaces) plus removal of 'hasEfficiency'
%
WordText(ActXWord, 'General Data', head, [0, 1]);
WordIndent(ActXWord, 7.5, -7.5, 'No Spacing')
WordText(ActXWord, ['Class name:', char(9), char(vehicle.class)], [],[0, 1], [], 1);
WordText(ActXWord, ['Maximum vehicle speed:', char(9), num2str(vehicle.vMax * [3.6, 3.6 / 1.609344, 1], '%2.1f km/h (%2.1f mph, %2.2f m/s)')], [], [0, 1], [], 1);
WordText(ActXWord, ['Vehicle length:', char(9), num2str(vehicle.unitLength, '%2.2f m')], [], [0, 1], [], 1);
if level > 1
    if vehicle.uFreq > 0 % AC Line
        if vehicle.U_Line_Nom > 2000
            str = num2str([[vehicle.U_Line_Nom, vehicle.MinUL, vehicle.OVClevel] / 1000, vehicle.uFreq], '%2.1f kV (%2.1f/%2.1f kV), AC %2.1f Hz');
        else
            str = num2str([vehicle.U_Line_Nom, vehicle.MinUL, vehicle.OVClevel, vehicle.uFreq], '%2.1f V (%2.1f/%2.1f V), AC %2.1fHz');
        end
    else
        if vehicle.U_Line_Nom > 2000
            str = num2str([vehicle.U_Line_Nom, vehicle.MinUL, vehicle.OVClevel] / 1000, 'DC %2.1f kV (%2.1f/%2.1f kV)');
        else
            str = num2str([vehicle.U_Line_Nom, vehicle.MinUL, vehicle.OVClevel], 'DC %2.1f V  (%2.1f/%2.1f V)');
        end
    end
    WordText(ActXWord,['Line voltage, nominal (power/regen):', char(9), str], [], [0, 1], [], 1);

    if vehicle.APS_efficiency == 1
        WordText(ActXWord, ['Auxiliary converter output power:', char(9), num2str(vehicle.AuxPower, '%2.1f kW')], [], [0, 1], [], 1);
    else
        WordText(ActXWord, ['APS input/output power:', char(9), num2str([vehicle.AuxPower ./ [vehicle.APS_efficiency, 1], vehicle.APS_efficiency * 100], '%2.1f / %2.1f kW (%2.1f%% efficiency)')], [], [0, 1], [], 1);
    end
    if vehicle.AuxPF ~= 1
        WordText(ActXWord, ['Auxiliary converter rating:', char(9), num2str([vehicle.AuxPower / vehicle.AuxPF, vehicle.AuxPF], '%2.1f kVA (Power Factor=%2.2f)')], [], [0, 1], [], 1);
    end
    if ~isempty(vehicle.AmbientTemp)
        WordText(ActXWord, ['Ambient temperature:', char(9), num2str(vehicle.AmbientTemp, ['%2.1f', char(176), 'C'])], [], [0, 1], [], 1);
    end
    if vehicle.ILDriveMax < inf
        WordText(ActXWord, ['Line current limitation, traction:', char(9), num2str(vehicle.ILDriveMax, '%2.2f A')], [], [0, 1], [], 1);
    end
    if vehicle.ILBrakeMax < inf
        WordText(ActXWord, ['Line current limitation, braking:', char(9), num2str(vehicle.ILBrakeMax, '%2.2f A')], [], [0, 1], [], 1);
    end
end
WordText(ActXWord, '', [], [0, 1], [], 0);

% --------- Architecture ---------
if level > 1
    WordText(ActXWord, 'Architecture', head, [0, 1]);
    WordText(ActXWord, '', 'No Spacing');
    WordIndent(ActXWord, 7.5, -7.5);
    WordText(ActXWord, ['Train configuration:', char(9)]);
    if vehicle.nUnits == 1
        WordText(ActXWord, 'Single unit operation', [], [0, 1], [], 1);
    elseif vehicle.nUnits == 2
        WordText(ActXWord, 'Two units in multiple operation', [], [0, 1], [], 1);
    else
        WordText(ActXWord, [num2str(vehicle.nUnits) ' units in multiple operation'], [], [0, 1], [], 1);
    end
    WordText(ActXWord, ['Number of High Voltage systems:', char(9) num2str(vehicle.nHVSystems) ' per unit'], [], [0, 1], [], 1);
    WordText(ActXWord, ['Number of axles:', char(9) num2str(vehicle.nAxDrive+vehicle.nAxTrail) ', of which ' num2str(vehicle.nAxDrive) ' are driven'], [], [0, 1], [], 1);
    if vehicle.nHVSystems==1
        WordText(ActXWord, ['Number of Motor Inverters:', char(9) num2str(vehicle.nINV) ' per unit'], [], [0, 1], [], 1);
        WordText(ActXWord, ['Number of Traction motors:', char(9) num2str(vehicle.nMotors) ' per Motor inverter, ' num2str(vehicle.nMotors*vehicle.nINV*vehicle.nHVSystems) ' in total per unit'], [], [0, 1], [], 1);
        if vehicle.uFreq>0
            WordText(ActXWord, ['Number of Line Converters:', char(9) num2str(vehicle.nCONV) ' per unit'], [], [0, 1], [], 1);
        end
        WordText(ActXWord, ['Number of APS Inverters:', char(9) num2str(vehicle.nAPUs) ' per unit'], [], [0, 1], [], 1);
    else
        WordText(ActXWord, ['Number of Motor Inverters:', char(9) num2str(vehicle.nINV) ' per HV system, ' num2str(vehicle.nINV*vehicle.nHVSystems) ' in total per unit'], [], [0, 1], [], 1);
        WordText(ActXWord, ['Number of Traction motors:', char(9) num2str(vehicle.nMotors) ' per Motor inverter, ' num2str(vehicle.nMotors*vehicle.nINV*vehicle.nHVSystems) ' in total per unit'], [], [0, 1], [], 1);
        if vehicle.uFreq>0
            WordText(ActXWord, ['Number of Line Converters:', char(9) num2str(vehicle.nCONV) ' per HV system, ' num2str(vehicle.nCONV*vehicle.nHVSystems) ' in total per unit'], [], [0, 1], [], 1);
        end
        WordText(ActXWord, ['Number of APS Inverters:', char(9) num2str(vehicle.nAPUs) ' per HV system, ' num2str(vehicle.nAPUs*vehicle.nHVSystems) ' in total per unit'], [], [0, 1], [], 1);
    end
    WordText(ActXWord, '', [], [0, 1], [], 0);
end

% --------- Performance ---------

WordText(ActXWord, 'Traction Performance', head, [0, 1]);
WordText(ActXWord, '', 'No Spacing');
WordIndent(ActXWord, 7.5, -7.5)
fmax=min(max(vehicle.PwrTE_f), vehicle.maxTE);
WordText(ActXWord, ['Max initial Tractive Effort:', char(9), num2str(fmax, '%2.1f kN')], [], [0, 1], [], 1);
WordText(ActXWord, ['Base Speed:', char(9), num2str(vehicle.vBasePwr * [3.6, 3.6 / 1.609344, 1], '%2.1f km/h (%2.1f mph, %2.2f m/s)')], [], [0, 1], [], 1);
WordText(ActXWord, ['Max Power (at base speed):', char(9) num2str(min(max(vehicle.PwrTE_p), vehicle.maxPwrTEbase), '%2.1f kW')], [], [0, 1], [], 1);
if vehicle.vBasePwr2<vehicle.vMax
    WordText(ActXWord, ['Weak field speed:', char(9), num2str(vehicle.vBasePwr2 * [3.6, 3.6 / 1.609344,  1], '%2.1f km/h (%2.1f mph, %2.2f m/s)')], [], [0, 1], [], 1);
end
if size(vehicle.PwrTE_f, 2) > 2
    WordText(ActXWord, '', [], [0, 1], [], 0);
    WordText(ActXWord, 'Tractive Effort curve:', [], [0, 1])
    data={'Speed [km/h]', 'Effort [kN]', 'Power [kW]'};
    for ii = 1:size(vehicle.PwrTE_f, 2)
        if ~isnan(vehicle.PwrTE_f)
            data{ii+1, 1}=num2str(vehicle.PwrTE_v(1, ii) * 3.6, '%8.2f');
            data{ii+1, 2}=num2str(vehicle.PwrTE_f(1, ii), '%8.2f');
            data{ii+1, 3}=num2str(vehicle.PwrTE_p(1, ii), '%8.1f');
        end
    end
    WordCreateTable(ActXWord, data, {'Centered'}, [0, 0]);
    WordCaption(ActXWord, 'Table', 'Tractive Effort curve', 1);
%     WordText(ActXWord,': Tractive Effort curve',[],[0, 1],[],0);
    WordIndent(ActXWord, 7.5, -7.5, 'No Spacing');
end
if vehicle.jerkrateT<inf
    WordText(ActXWord, ['Jerk rate in traction mode:', char(9), num2str(vehicle.jerkrateT, '%2.2f m/s^3')], [], [0, 1], [], 1);
end
if vehicle.maxRatePwrTE<inf
    WordText(ActXWord, ['Max rail power application rate:', char(9), num2str(vehicle.maxRatePwrTE, '%2.2f kW/s (per unit)')], [], [0, 1], [], 1);
end
loadCond = vehicle.LoadCondition;
if vehicle.max_acc < inf
    if vehicle.AccMethod < 2
        Fbr=vehicle.totMass(loadCond) * vehicle.max_acc;
        if vehicle.AccMethod == 1
            Fbr = Fbr + vehicle.davies_a(loadCond);
            a=min(vehicle.max_acc, (fmax - vehicle.davies_a(loadCond)) / vehicle.totMass(loadCond));
            WordText(ActXWord, ['Max initial acceleration:', char(9), num2str(a, '%2.2f m/s^2')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Max utilised Tractive effort:', char(9), num2str(min(fmax, Fbr), '%2.1f kN')], [], [0, 1], [], 1);
        else
            acc = (min(Fbr, fmax) - vehicle.davies_a(loadCond)) / vehicle.totMass(loadCond);
            WordText(ActXWord, ['Max initial acceleration:', char(9), num2str(acc, '%2.2f m/s^2')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Max utilised Tractive effort:', char(9), num2str(min(Fbr, fmax), '%2.1f kN')], [], [0, 1], [], 1);
        end
    else
        WordText(ActXWord, ['Max acceleration:', char(9), num2str(vehicle.max_acc, '%2.2f m/s^2')], [], [0, 1], [], 1);
    end
else
    WordText(ActXWord, ['Max acceleration:' char(9), 'No limit applied'], [], [0, 1], [], 1);
end
WordText(ActXWord, ['Max Dynamic Brake Effort:', char(9), num2str(min(max(vehicle.PwrDBE_f), vehicle.maxDBE), '%2.1f kN')], [], [0, 1], [], 1);
WordText(ActXWord, ['Max Dynamic Brake Power:', char(9), num2str(min(max(vehicle.PwrDBE_p), vehicle.maxPwrDBEbase), '%2.1f kW')], [], [0, 1], [], 1);
WordText(ActXWord, ['Dynamic Brake Base Speed:', char(9), num2str(vehicle.vBaseBrk * [3.6, 3.6 / 1.609344, 1], '%2.1f km/h (%2.1f mph, %2.2f m/s)')], [], [0, 1], [], 1);
if vehicle.vBaseBrk2 < vehicle.vMax
    WordText(ActXWord, ['Dynamic Brake Weak field speed:', char(9), num2str(vehicle.vBaseBrk2 * [3.6, 3.6 / 1.609344, 1], '%2.1f km/h (%2.1f mph, %2.2f m/s)')], [], [0, 1], [], 1);
end
if size(vehicle.PwrDBE_f, 2) > 2
    WordText(ActXWord, '', [], [0, 1], [], 0);
    WordText(ActXWord, 'Dynamic Braking curve:', [], [0, 1])
    data={'Speed [km/h]', 'Effort [kN]', 'Power [kW]'};
    for ii = 1:size(vehicle.PwrDBE_v, 2)
        if ~isnan(vehicle.PwrDBE_v)
            data{ii+1, 1}=num2str(vehicle.PwrDBE_v(1, ii)*3.6, '%8.2f');
            data{ii+1, 2}=num2str(vehicle.PwrDBE_f(1, ii), '%8.2f');
            data{ii+1, 3}=num2str(vehicle.PwrDBE_p(1, ii), '%8.1f');
        end
    end
    WordCreateTable(ActXWord, data, {'Centered'}, [0, 1]);
    WordCaption(ActXWord, 'Table', 'Dynamic brake curve', 1);
%     WordText(ActXWord, ': Dynamic brake curve', [], [0, 1], [], 0);
end
WordText(ActXWord, '', [], [0, 1], [], 0);

% --------- Brake Settings ---------

WordText(ActXWord, 'Brake definition', head, [0, 1]);
WordIndent(ActXWord, 7.5, -7.5, 'No Spacing')

if vehicle.RetardMethod < 2
    Fbr = vehicle.totMass(loadCond) * vehicle.retard_max;
    if vehicle.RetardMethod == 1
        Fbr = Fbr - vehicle.davies_a(loadCond);
    end
    if vehicle.utiliseEBrake
        spd = vehicle.retard_min / vehicle.retard_max * vehicle.retard_min_speed;
        WordText(ActXWord, ['Brake Effort - Constant effort:', char(9), num2str(Fbr, '%2.1f kN (up to Brake Base Speed)')], [], [0, 1], [], 1);
        WordText(ActXWord, ['Total Brake Base Speed:', char(9), num2str(spd * [3.6, 3.6 / 1.609344, 1], '%2.1f km/h (%2.1f mph, %2.2f m/s)')], [], [0, 1], [], 1);
        WordText(ActXWord, ['Total Brake Power:', char(9), num2str(spd*Fbr, '%2.1f kW')], [], [0, 1], [], 1);
    else
        WordText(ActXWord, ['Brake Effort - Constant effort:', char(9), num2str(Fbr, '%2.1f kN (all speeds)')], [], [0, 1], [], 1);
    end
else
    if vehicle.utiliseEBrake
        spd = vehicle.retard_min / vehicle.retard_max * vehicle.retard_min_speed;
        WordText(ActXWord, ['Constant retardation:', char(9), num2str(vehicle.retard_max, '%2.3f m/s^2')], [], [0, 1], [], 1);
        WordText(ActXWord, ['Retardation Corner Speed:', char(9), num2str(spd * [3.6, 3.6 / 1.609344, 1], '%2.1f km/h (%2.1f mph, %2.2f m/s)')], [], [0, 1], [], 1);
        WordText(ActXWord, ['Retardation constant (R = C x v):', char(9), num2str(vehicle.retard_max * spd, '%2.1f m^2/s^3')], [], [0, 1], [], 1);
    else
        WordText(ActXWord, ['Constant retardation:', char(9), num2str(vehicle.retard_max, '%2.3f m/s^2 (at all speeds)')], [], [0, 1], [], 1);
    end

end
if size(vehicle.TotBrake_f, 2) > 2
    WordText(ActXWord, '', 'No Spacing', [0, 1], [], 0);
    WordText(ActXWord, 'The table below shows the "Total Braking curve":', [], [0, 1])
    data = {'Speed [km/h]', 'Effort [kN]', 'Power [kW]'};
    for ii = 1:size(vehicle.TotBrake_f, 2)
        if ~isnan(vehicle.TotBrake_f(ii))
            data{ii+1, 1}=num2str(vehicle.TotBrake_v(1, ii)*3.6, '%8.2f');
            data{ii+1, 2}=num2str(vehicle.TotBrake_f(1, ii), '%8.2f');
            data{ii+1, 3}=num2str(vehicle.TotBrake_p(1, ii), '%8.1f');
        end
    end
    WordCreateTable(ActXWord, data, {'Centered'}, [0, 0]);
    WordCaption(ActXWord, 'Table', 'Total brake curve', 1);
end
WordText(ActXWord, '', 'No Spacing', [0, 1], [], 0);
WordIndent(ActXWord, 7.5, -7.5, 'No Spacing')
if vehicle.jerkrateB < inf
    WordText(ActXWord, ['Brake Jerk rate:', char(9), num2str(vehicle.jerkrateB, '%2.2f m/s^3')], [], [0, 1], [], 1);
end
if vehicle.residual_retard > 0
    WordText(ActXWord, ['Residual Retardation (just before stop):', char(9), num2str(vehicle.residual_retard, '%2.2f m/s^2')], [], [0, 1], [], 1);
end

WordText(ActXWord, ['Brake fading speed, start / finish:', char(9), num2str([vehicle.dynBrkFadeStart, vehicle.dynBrkFadeEnd] * 3.6, '%2.1f / %2.1f km/h')], [], [0, 1], [], 1);
WordText(ActXWord, '', [], [0, 1], [], 0);

% --------- Vehicle Masses ---------

WordText(ActXWord, 'Vehicle Mass', head, [0, 1]);
WordText(ActXWord, '', 'No Spacing');
WordIndent(ActXWord, 7.5, -7.5)

WordText(ActXWord, ['Vehicle static mass:', char(9), num2str((vehicle.statmassLoad(loadCond) - vehicle.LoadMassUnit{loadCond}), '%2.3f ton')], [], [0, 1], [], 1);
WordText(ActXWord, ['Vehicle load mass:', char(9), num2str(vehicle.LoadMassUnit{loadCond}, '%2.3f ton') ' (' ...
    vehicle.LoadMassName{loadCond} ': "' vehicle.LoadMassDef{loadCond} '")'], [], [0, 1], [], 1);
WordText(ActXWord, ['Rotation Inertia, equivalent mass:', char(9), num2str(vehicle.totMass(loadCond) - vehicle.statmassLoad(loadCond), '%2.3f ton')], [], [0 2], [], 1);

if level > 1
    WordText(ActXWord, '', 'No Spacing'); %Restore the tabs
    
    data       = vehicle.LoadMassName';
    data(:, 2) = vehicle.LoadMassDef';
    data(:, 3) = cellstr(num2str([vehicle.LoadMassUnit{:}]', '%2.1f'));
    data(:, 4) = cellstr(num2str(vehicle.statmassLoad(loadCond) + [vehicle.LoadMassUnit{:}]', '%2.1f'));
    data(:, 5) = cellstr(num2str(vehicle.totMass', '%2.1f'));
    data(2:end + 1, :) = data;
    data(1, :) = {'Cond.', 'Definition', 'Load weight', 'Total static mass', 'Total dynamic mass'};
    WordCreateTable(ActXWord, data, {'Left', 'Left', 'Centered', 'Centered', 'Centered'}, [0, 0]); % enter before table
    WordCaption(ActXWord, 'Table', 'Defined load weights and masses', 1);
    % ActXWord.Selection.InsertCaption(-2);
%     WordText(ActXWord,': Defined load weights and masses',[],[0, 1],[],0);
end

% --------- Train Resistance ---------

WordText(ActXWord, 'Train Resistance', head, [0, 1]);
WordText(ActXWord, ['The train resistance is defined for the train in ' char(vehicle.LoadMassName(vehicle.LoadCondition)) ' load condition. ']);
WordText(ActXWord, 'The total train resistance is calculated using Davies formula as shown below:', [], [0, 1], [], 1);
WordTypeDavies(ActXWord)
% WordText(ActXWord,'','Normal');
% WordText(ActXWord, '', 'No Spacing');
WordIndent(ActXWord, 7.5, -7.5, 'No Spacing')
WordText(ActXWord, ['Davies A:', char(9), num2str(vehicle.davies_a(vehicle.LoadCondition) * 1000, '%2.0f N')], [], [0, 1], [], 1);
if level > 1
    WordIndent(ActXWord, 7.5, -7)
    WordText(ActXWord, ['Mass dependency, A'':' char(9), num2str(vehicle.davies_a_dash, '%2.1f N/ton')], [], [0, 1], [], 1);
    WordIndent(ActXWord, 7.5, -7.5)
end
WordText(ActXWord, ['Davies B:', char(9) num2str(vehicle.davies_b(vehicle.LoadCondition) * 1000, '%2.2f N/(km/h)')], [], [0, 1], [], 1);
if level > 1
    WordIndent(ActXWord, 7.5, -7)
    WordText(ActXWord, ['Mass dependency, B_1'':' char(9), num2str(vehicle.davies_b1_dash, '%2.3f N/(t x (km/h))')], [], [0, 1], [], 1);
    WordText(ActXWord, ['Fixed part, B_2:' char(9), num2str(vehicle.davies_b2, '%2.1f N/(km/h)')], [], [0, 1], [], 1);
    WordIndent(ActXWord, 7.5, -7.5)
end
WordText(ActXWord, ['Davies C, open condition:', char(9), num2str(vehicle.davies_c1(1) * 1000, '%2.3f N/(km/h)^2')], [], [0, 1], [], 1);
if level > 1
    WordIndent(ActXWord, 7.5, -7)
    WordText(ActXWord, ['Front contribution, C''_{Front}:' char(9), num2str(vehicle.davies_c_front, '%2.3f N/(km/h)^2')], [], [0, 1], [], 1);
    WordText(ActXWord, ['Length contribution, D_{Length}:' char(9), num2str(vehicle.davies_c_len * 1000, '%2.2f x 10^{-3} N/(m x (km/h)^2)')], [], [0, 1], [], 1);
    WordIndent(ActXWord, 7.5, -7.5)
end
WordText(ActXWord, ['Davies C, in tunnels:', char(9), num2str(vehicle.davies_c2(1) * 1000, '%2.3f N/(km/h)^2')], [], [0, 1], [], 1);
if level > 1
    WordIndent(ActXWord, 7.5, -7)
    if vehicle.davies_c1(1) > 0
        WordText(ActXWord, ['Tunnel factor (single unit):', char(9), num2str(vehicle.davies_c2(1) / vehicle.davies_c1(1), '%2.2f')], [], [0, 1], [], 1);
    else
        WordText(ActXWord, ['Tunnel factor (single unit):', char(9), 'N/A'], [], [0, 1], [], 1);
    end
    WordIndent(ActXWord, 7.5, -7.5)
end
if vehicle.curve_cr0 > 0 || level > 1
    WordText(ActXWord, ['Curve resistance, C_{R0}:', char(9), num2str(vehicle.curve_cr0, '%2.1f N/kgm')], [], [0, 1], [], 1);
    WordText(ActXWord, ['Curve radius reduction, C_{R1}:', char(9), num2str(vehicle.curve_cr1, '%2.1f m')], [], [0, 1], [], 1);
end
if vehicle.curve_cr2 > 0 || level > 1
    WordText(ActXWord, ['Speed dependent curve resistance, C_{R2}:', char(9), num2str(vehicle.curve_cr2, '%2.4f Ns^2/kgm')], [], [0, 1], [], 1);
end
WordText(ActXWord, '', [], [0, 1], [], 0);


% --------- Performance graphs ---------

if ishandle(ActXWord)
    if level > 1
        WordText(ActXWord, 'Performance Graphs', head, [0, 1]);
        WordText(ActXWord, ['The following section contains the figures describing the train performance in terms of Efforts, Line- & Rail Power and resulting Acceleration in '...
            vehicle.LoadMassName{vehicle.LoadCondition} ' load condition.'], [], [0, 1]);
        % ActXWord.Selection.TypeText(['The following section contains the figures describing the train train performance in terms of Efforts, Line- & Rail Power and resulting Acceleration'...
        %     char(vehicle.LoadMassName(vehicle.load_cond)) ' load condition.']);
        % ActXWord.Selection.TypeParagraph;
        ShSim_Plot_Vehicle('a', data_main, [], sim_config); % Efforts
        figure(200);
        drawnow;
        print(gcf, '-dmeta')
        close(200)
        pause(1)
        ActXWord.Selection.Paste
        WordCaption(ActXWord, 'Figure', 'Tractive-, Dynamic Brake and Friction Brake Efforts vs. Speed', 1);
        %     WordText(ActXWord,': Tractive-, Dynamic Brake and Friction Brake Efforts vs. Speed',[],[0, 1]);
        
        ShSim_Plot_Vehicle('c', data_main, [], sim_config); % Power vs Speed
        figure(202);
        drawnow;
        print(gcf, '-dmeta')
        close(202)
        pause(1)
        ActXWord.Selection.Paste
        WordCaption(ActXWord, 'Figure', 'Line and Rail Power vs. Speed at Powering and Braking', 1);
        %     WordText(ActXWord,': Line and Rail Power vs. Speed at Powering and Braking',[],[0, 1]);
        
        ShSim_Plot_Vehicle('b', data_main, [], sim_config); % Acceleration
        figure(201);
        drawnow;
        print(gcf, '-dmeta')
        close(201)
        pause(1)
        ActXWord.Selection.Paste
        WordCaption(ActXWord, 'Figure', 'Acceleration vs. Speed at Powering and Braking', 1);
        %     WordText(ActXWord,': Acceleration vs. Speed at Powering and Braking',[],[0, 1]);
    end
end

% --------- Drive System ---------
if level > 1
    WordText(ActXWord, 'Drive System', head, [0, 1]);
    WordText(ActXWord, 'The drive system consists of the wheels, gearbox and the traction motor.', [], [0, 1], [], 1);
    WordText(ActXWord, '', 'No Spacing', [0, 0]);
    WordIndent(ActXWord, 7.5, -7.5)
    if vehicle.uFreq > 0
        WordText(ActXWord, ['DC-link voltage:', char(9), num2str(vehicle.uDcLink, '%2.0f VDC')], [], [0, 1], [], 1);
    end
    WordText(ActXWord, ['Wheel diameters, Nominal (New/Worn):', char(9), num2str([vehicle.wheelsize{1}, vehicle.wheelsize{3}, vehicle.wheelsize{2}] * 1000, '%2.0f mm (%2.0f/%2.0f mm)')], [], [0, 1], [], 1);
    sel = {'Nominal', 'Worn', 'New'};
    WordText(ActXWord, ['Default wheel size: ' char(9), char(sel(vehicle.ws_selection)), ' (', num2str(vehicle.wheelsize{vehicle.ws_selection} * 1000, '%2.0f mm'), ')'], [], [0, 1], [], 1);
    WordText(ActXWord, ['Gear ratio:', char(9), num2str(vehicle.gearRatio, '1:%2.3f')], [], [0, 1], [], 1);
    WordText(ActXWord, '', [], [0, 1], [], 0);
end

% --------- Drive System Losses ---------
if level > 2
    WordText(ActXWord, 'Drive System Losses', head2, [0, 1]);
    WordText(ActXWord, '', 'No Spacing', [0, 0]);
    WordIndent(ActXWord, 7.5, -7.5)
    if vehicle.gearLoss0 ~= 0
        WordText(ActXWord, ['Gear loss, constant torque:', char(9), num2str(vehicle.gearLoss0, '%2.3f Nm')], [], [0, 1], [], 1);
    end
    if vehicle.gearLoss1 ~= 0
        WordText(ActXWord, ['Gear loss, torque relative:', char(9), num2str(vehicle.gearLoss1*100, '%2.2f%%')], [], [0, 1], [], 1);
    end
    if vehicle.gearLoss2 ~= 0
        WordText(ActXWord, ['Gear loss, speed relative:', char(9), num2str(vehicle.gearLoss2*1000, '%2.2f mNm/rpm')], [], [0, 1], [], 1);
    end

    % --------- Traction motor ---------

    switch vehicle.motorMethod
        case 1 % Fixed losses
            WordText(ActXWord, ['Traction motor fix losses:', char(9), num2str(vehicle.TM_FixLoss * 100, '%2.2f%%')], [], [0, 1], [], 1);
        case 2 % Motor map
            WordText(ActXWord, '', 'No Spacing', [0, 1], 0);
            WordText(ActXWord, 'Traction Motor - Motor map', head2, [0, 1]);
            WordText(ActXWord, '', 'No Spacing');
            WordIndent(ActXWord, 7.5, -7.5)

            WordText(ActXWord, ['Max speed of motor:', char(9), num2str(vehicle.motorMap.VMax, '%2.1f rpm')], [], [0, 1], [], 1);
            WordText(ActXWord, ['DC-link voltage in motor operation:', char(9), num2str(vehicle.motorMap.UDC_motor, '%2.1f V')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Base speed in motor operation:', char(9), num2str(vehicle.motorMap.rpmBaseMotor, '%2.1f rpm')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Max torque in motor operation:', char(9), num2str(vehicle.motorMap.TorqueMT(end, end), '%2.1f Nm')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Max power in motor operation:', char(9), num2str(vehicle.motorMap.TorqueMT(end, end)*vehicle.motorMap.rpmBaseMotor*2*pi/60000, '%2.1f kW')], [], [0, 1], [], 1);
            WordText(ActXWord, ['DC-link voltage in generator mode:', char(9), num2str(vehicle.motorMap.UDC_generator, '%2.1f V')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Base speed in generator mode:', char(9), num2str(vehicle.motorMap.rpmBaseGenerator, '%2.1f rpm')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Max torque in generator mode:', char(9), num2str(vehicle.motorMap.TorqueGT(end, end), '%2.1f Nm')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Max power in generator mode:', char(9), num2str(-vehicle.motorMap.TorqueGT(end, end)*vehicle.motorMap.rpmBaseGenerator*2*pi/60000, '%2.1f kW')], [], [0, 1], [], 1);
            if ishandle(ActXWord)
                ShSim_Plot_Vehicle('d', data_main, [], sim_config); % Efforts
                figure(203);
                drawnow;
                print(gcf, '-dmeta')
                close(203)
                pause(1)
                ActXWord.Selection.Paste
                WordCaption(ActXWord, 'Figure', 'Torque vs. Speed', 1);
            end
    end % switch
    if vehicle.motorMethod > 1
        WordText(ActXWord, '', [], [0, 1]);
        WordText(ActXWord, 'Other losses', head, [0, 1], [], 1);
        WordIndent(ActXWord, 7.5, -7.5, 'No Spacing')
    end
    if vehicle.R_MotCable > 0
        WordText(ActXWord, ['Motor cable equivalent resistance:', char(9), num2str(vehicle.R_MotCable * 1000, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
    end
    if isempty(vehicle.INV_FixLoss)
        WordText(ActXWord, ['Inverter resistive loss:', char(9), num2str(vehicle.INV_RLoss * 1000, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
        WordText(ActXWord, ['Inverter switching loss:', char(9), num2str(vehicle.INV_FLoss, '%2.2f J/switch')], [], [0, 1], [], 1);
        WordText(ActXWord, ['Inverter switching resistivity:', char(9), num2str(vehicle.INV_FRLoss * 1e6, ['%2.2f ', char(956), char(8486), '/Hz'])], [], [0, 1], [], 1);
    else
        WordText(ActXWord, ['Traction Inverter fix losses:', char(9), num2str(vehicle.INV_FixLoss * 100, '%2.2f%%')], [], [0, 1], [], 1);
    end
    if vehicle.uFreq > 0
        if isempty(vehicle.CONV_FixLoss)
            WordText(ActXWord, ['Converter resistive loss:', char(9), num2str(vehicle.CONV_RLoss * 1000, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
            WordText(ActXWord, ['Converter switching loss:', char(9), num2str(vehicle.CONV_FLoss, '%2.2f J/switch')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Converter switching resistivity:', char(9), num2str(vehicle.CONV_FRLoss * 1e6, ['%2.2f ', char(956), char(8486), '/Hz'])], [], [0, 1], [], 1);
        else
            WordText(ActXWord, ['Converter fix losses:', char(9), num2str(vehicle.CONV_FixLoss * 100, '%2.2f%%')], [], [0, 1], [], 1);
        end
        if vehicle.R_ConCable > 0
            WordText(ActXWord, ['Cable resistance Converter to Transformer:', char(9), num2str(vehicle.R_ConCable  *1000, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
        end
        if isempty(vehicle.Trafo_FixLoss)
            WordText(ActXWord, 'Main transformer:', 'Underline', [0, 1], [], 1);
            WordIndent(ActXWord, 7.5, -7.5, 'No Spacing')
            WordText(ActXWord, ['Main ratio:', char(9), num2str(vehicle.MT_ratio, '1:%2.2f')], [], [0, 1], [], 1);
            WordText(ActXWord, ['No load losses:', char(9), num2str(vehicle.trafoIdleLoss / 1000, '%2.2f kW')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Harmonic losses:', char(9), num2str(vehicle.trafoHarmonicLoss / 1000, '%2.2f kW')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Iron loss:', char(9), num2str(vehicle.trafoIronLoss * 100, '%2.2f%%')], [], [0, 1], [], 1);
            WordText(ActXWord, ['Primary winding resistance:', char(9), num2str(vehicle.rTrafoPrim * 1000, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
            WordText(ActXWord, ['Secondary winding resistance:', char(9), num2str(vehicle.rTrafoSec * 1000, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
            WordText(ActXWord, ['Transformer temperature:', char(9), num2str(vehicle.TrafoTemp, ['%2.0f', char(176), 'C'])], [], [0, 1], [], 1);
            WordText(ActXWord, ['Equivalent resistance line side:', char(9), num2str(vehicle.rTrafoWind, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
        else
            WordText(ActXWord, ['Main transformer fix losses:', char(9), num2str(vehicle.Trafo_FixLoss * 100, '%2.2f%%')], [], [0, 1], [], 1);
        end
    else
        if vehicle.rLineInd>0
            WordText(ActXWord, ['Line filter resistance:', char(9), num2str(vehicle.rLineInd * 1000, ['%2.2f m', char(8486)])], [], [0, 1], [], 1);
        end
    end

    WordText(ActXWord, '', [], [0, 1], [], 0);

    % ------------- Insert Efficiency figures ------------
    if ishandle(ActXWord)
        WordText(ActXWord, 'Efficiency plot', head, [0, 1]);
        ShSim_Plot_Vehicle('e', data_main, [], sim_config); % Efforts
        figure(204);
        drawnow;
        print(gcf, '-dmeta')
        close(204)
        pause(1)
        ActXWord.Selection.Paste
        ActXWord.Selection.InsertCaption(-1);
        WordText(ActXWord, ': Total System Efficiency vs. Speed', [], [0, 1]);
    end
end
return % ShSim_VehicleDetailedReport

function WordTypeDavies(ActXWord)
if ishandle(ActXWord)
    ActXWord.Selection.Style = 'Normal';
    ActXWord.Selection.Font.Size = 14;
    ActXWord.Selection.Font.Italic = 1;
    ActXWord.Selection.ParagraphFormat.Alignment = 1; % Centre
    WordText(ActXWord, ['F_{Res} ', char(61), 'A + B ', char(215), ' v + C ', char(215), ' v^2'], [], [0, 1], [], 1)
    WordText(ActXWord, 'Where ', 'Normal');
    ActXWord.Selection.Font.Size = 14;
    ActXWord.Selection.Font.Italic = 1;
    WordText(ActXWord, 'v');
    ActXWord.Selection.Font.Size = 11;
    ActXWord.Selection.Font.Italic = 0;
    WordText(ActXWord, ' is the train speed in km/h and ');
    ActXWord.Selection.Font.Size = 14;
    ActXWord.Selection.Font.Italic = 1;
    WordText(ActXWord, 'F_{Res}')
    ActXWord.Selection.Font.Size = 11;
    ActXWord.Selection.Font.Italic = 0;
    WordText(ActXWord, ' is the result total train resistance.', [], [0, 1]);
end
return % WordTypeDavies