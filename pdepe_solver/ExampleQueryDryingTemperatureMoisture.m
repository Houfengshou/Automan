% 进入本文件所在目录，然后运行本示例。
% 与Process2类似，T和C的行对应时间、列对应半径。
queryTimeSeconds=[100,300,600,1800,3600,10800,14400,14400.5,18000];
queryRadiusCm=[0,0.5,1,1.5,2];
[T,C,solution]=SolveDryingTemperatureMoisturePDEPE(queryTimeSeconds,queryRadiusCm);
disp('温度T（摄氏度）：'); disp(T);
disp('含水率C（kg/kg）：'); disp(C);

% 求解后再次查询：此处时间1234.5秒由保存结果插值。
[temperatureAt1234Seconds,moistureAt1234Seconds]= ...
    EvaluateDryingTemperatureMoisture(solution,1234.5,[0,0.7,2]);

% 第三问完整运行示例（需要时取消注释）：
% settings=struct('StopWhenDry',true,'MeshRefinement',1);
% [Tdry,Cdry,dryingSolution]=SolveDryingTemperatureMoisturePDEPE( ...
%     0:60:144*3600,0:0.1:2,settings);
% % 144小时只是本次搜索上限，未达标时需要扩大，不是预设答案。
% disp(dryingSolution.dryingThresholdTimeHours);
% ExportDryingMoistureResultTables(dryingSolution,fullfile(pwd,'results'));
