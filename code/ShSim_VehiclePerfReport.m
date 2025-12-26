function ShSim_VehiclePerfReport(data_main,sim_config)
% ShSim_VehiclePerfReport
% 0.1.1.2
%
% Ver       Date        Sign    Descr
% 0.1.1.2   2025-12-09  JR      Clean code, remove vehicle path
%
msg = 'Select output for your Performance Report';
answer = questdlg(msg, 'Performance Report', 'Command window', 'Word file', 'Command window');
if strcmp(answer, 'Command window')
    ActXWord = [];
else
    FileSpec = fullfile(pwd, 'Templates', 'Train Performance Template.docx');
    if exist(FileSpec, 'file')
        [ActXWord, WordHandle] = WordStart(FileSpec);
    else
        return
    end
    % Ask for document number
    defFile = ['TRG', char(datetime, 'yy'), '5xxxPROJ'];
    defFile = char(inputdlg({'Enter document No'}, 'Document number', [1, 35], {defFile}));
    if isempty(defFile)
        defFile = 'DocNo';
    end
    % Save the word document under a new name first
    old = cd(sim_config.workingDirectory);
    if exist('Documentation', 'dir')
        cd('Documentation');
    end
    [f, p] = uiputfile('*.docx', 'Save performance reports as...', [defFile, '_R0 Traction Performance Report.docx']);
    cd(old);
    
    if f
        if exist(fullfile(p, f), 'file')
            delete(fullfile(p, f));
        end
    else
        ActXWord.Quit
        delete(ActXWord);
        return
    end
    if exist(fullfile(p, f), 'file')
        warndlg('The selected word file could not be deleted. Check if it is opened by the application or locked')
        ActXWord.Quit;
        delete(ActXWord);
        return
    else
        invoke(WordHandle, 'SaveAs2', fullfile(p, f));
    end
end

vehicle = ShSim_Read_Vehicle(data_main.vehicle_file, 0);
if ishandle(ActXWord)
    WordGoToBkMrk(ActXWord, 'Header_subtitle');
    WordText(ActXWord, char(vehicle.class), [], [0, 0]);
    WordGoToBkMrk(ActXWord, 'Front_subhead1');
    WordText(ActXWord, char(vehicle.class), [], [0, 0]);
    WordGoToBkMrk(ActXWord, 'Trainset');
    WordText(ActXWord,['"', char(vehicle.class), '"'], [], [0, 0]);
    WordGoToBkMrk(ActXWord, 'Trainset1');
    WordText(ActXWord, ['"', char(vehicle.class), '"'], [], [0, 0]);
    WordGoToBkMrk(ActXWord, 'Footer_DocNo');
    WordText(ActXWord, char(defFile), [], [0, 0]);
end

% --------- References ---------

WordGoToBkMrk(ActXWord, 'References');
% ActXWord.ActiveDocument.Bookmarks.Item('References').Select
WordText(ActXWord, ['Vehicle File Name:', [char(9), char(9), char(9)], strrep(data_main.vehicle_file.fName, '_', '\_')], [], [0, 1], [], 1);
WordText(ActXWord, ['Location:', char(9), strrep(data_main.vehicle_file.path, '_', '\_')], [], [0, 1], [], 1);
dr = dir(fullfile(data_main.vehicle_file.path, data_main.vehicle_file.fName));
WordText(ActXWord, ['Modification date:', char(9), dr.date], [], [0, 1], [], 1);
WordText(ActXWord, ['Created with ShogunSim version:', char(9), sim_config.ShSimVer], [], [0, 1], [], 0);

% --------- General data ---------

WordGoToBkMrk(ActXWord, 'Perf_Report');
ShSim_VehicleDetailedReport(ActXWord, vehicle,data_main, sim_config, 3, '1.1 Heading', '1.1.1 Heading')

% ----------- Update TOC -------------

if ishandle(ActXWord)
    ActXWord.Selection.GoTo(7, 3, 1, 'TOC');
    ActXWord.Selection.Delete;
    ActXWord.ActiveDocument.TablesOfContents.Add(ActXWord.Selection.Range, 1, 1, 3);
    WordHandle.Save
    ActXWord.Quit
    delete(ActXWord);
    
    winopen(fullfile(p, f))
else
    input('Press enter to continue', 's');
end
return % ShSim_VehiclePerfReport

