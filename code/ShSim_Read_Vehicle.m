function retVehicle = ShSim_Read_Vehicle(filename)
% ShSim_Read_Vehicle
% 0.1.1.2
%
% ShSim_Read_Vehicle reads the vehicle data defined by the excel file
% 'filename'.
% input
%  pathname - path name to vehicle file
%  filename - vehicle file name
%  alternative - 1 if it is to be seen as alternative vehicle
%  sim_config - configuration of simulation
%    sim_config.workingDirectory - selected working directory
% output
%  r - vehicle structure
%
% Ver       Date        Sign    Descr
% 0.1.1.1   2025-12-11  JR      removed path - filename contains both
% 0.1.1.2   2025-12-19  JR      Bugfix
%
%
retVehicle = [];
if exist(filename, 'file')
    try
        load(filename, 'vehicle_data')
        retVehicle = ShSim_calcVehicle(vehicle_data);
    catch
        dlg = warndlg('Could not load the vehicle data');
        uiwait(dlg);
    end
    retVehicle = ShSim_CalcMassRes(retVehicle);
end
return % ShSim_Read_Vehicle

