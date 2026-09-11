function ExportDryingMoistureResultTables(solution,outputDirectory)
% 以草稿文件名导出，最终精度仍需网格收敛与独立对照确认。
assert(solution.thresholdReached,'尚未达到阈值；请延长求解时间后再导出表5。');
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
finishSeconds=solution.dryingThresholdTimeSeconds;
minuteTimes=(60:60:finishSeconds)';
minuteTimes=unique([minuteTimes;finishSeconds]);
[~,minuteMoisture]=EvaluateDryingTemperatureMoisture(solution,minuteTimes,0:0.1:2);
sixHourTimes=(21600:21600:finishSeconds)';
sixHourTimes=unique([sixHourTimes;finishSeconds]);
[~,sixHourMoisture]=EvaluateDryingTemperatureMoisture(solution,sixHourTimes,0:0.5:2);
writecell([{'时间(s)/半径(cm)'},num2cell(round(0:0.1:2,1)); ...
    num2cell([minuteTimes,round(minuteMoisture,4)])], ...
    fullfile(outputDirectory,'result3_pdepe_draft.xlsx'),'Sheet','Sheet1');
writecell([{'时间(h)/半径(cm)'},num2cell(0:0.5:2); ...
    num2cell(round([sixHourTimes/3600,sixHourMoisture],4))], ...
    fullfile(outputDirectory,'table5_pdepe_draft.xlsx'),'Sheet','Sheet1');
save(fullfile(outputDirectory,'pdepe_full_precision.mat'),'solution');
end
