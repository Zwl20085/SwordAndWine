% ----------------------------------------------
% Script Recorded by ANSYS Electronics Desktop Version 2019.2.0
% 1:00:24  2ÔÂ 23, 2020
% 
% Script Revised by W.T.Zhang 2020.2.23 14:24
% 2020.2.23 14:24
% ----------------------------------------------

iMaxwell = actxserver('AnsoftMaxwell.MaxwellScript');
Desktop = iMaxwell.GetAppDesktop();
Desktop.RestoreWindow;
Project = Desktop.SetActiveProject('1210-210 - straight');
Project = Desktop.SetActiveProject('1210-210 - straight');
Design = Project.SetActiveDesign('check1');
Module = Design.GetModule('Optimetrics');
invoke(Module,'ImportSetup', 'OptiParametric', {'NAME:ParametricSetup1','D:\Scientific Research\icem 2020\MaxwellAutomation\ParametricSetup9_Profile1.csv'});
Project.Save;
invoke(Module,'SolveSetup', 'ParametricSetup1');
Module = Design.GetModule('ReportSetup');
invoke(Module,'ExportToFile', 'Torque','D:/Scientific Research/icem 2020/MaxwellAutomation/Torque.csv');
invoke(Module,'ExportToFile', 'Winding Voltage','D:/Scientific Research/icem 2020/MaxwellAutomation/Winding Voltage.csv');
invoke(Module,'ExportToFile', 'FluxLinkage','D:/Scientific Research/icem 2020/MaxwellAutomation/FluxLinkage.csv');

while(~exist('D:/Scientific Research/icem 2020/MaxwellAutomation/FluxLinkage.csv','file'))
    pause(5);
end
Module = Design.GetModule('AnalysisSetup');
invoke(Module,'ResetSetupToTimeZero','Setup1');
Module = Design.GetModule('Optimetrics');
invoke(Module,'DeleteSetups' ,{'ParametricSetup1'})
delete(iMaxwell);