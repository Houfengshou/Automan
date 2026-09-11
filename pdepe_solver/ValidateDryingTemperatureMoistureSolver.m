function report=ValidateDryingTemperatureMoistureSolver()
% 验证真实输入、任意查询、4小时衔接、阈值事件和数值加密。
rootDirectory=fileparts(mfilename('fullpath'));
timeSeconds=[0;100;1800;10800;14400;14400.5;18000];
[T,C,solution]=SolveDryingTemperatureMoisturePDEPE(timeSeconds,[0,0.5,1,1.5,2]);
assert(isequal(size(T),[7,5]) && isequal(size(C),[7,5]));
assert(max(abs(T(1,:)-28))<1e-12 && max(abs(C(1,:)-2.55))<1e-12);
assert(all(T(:)>=28-1e-5) && all(T(:)<51) && all(C(:)>0));
assert(all(diff(T(2,:),1,2)>=-1e-5) && all(diff(C(2,:),1,2)<=1e-5));
[Tq,Cq]=EvaluateDryingTemperatureMoisture(solution,1234.5,[0,0.7,2]);
[Td,Cd]=SolveDryingTemperatureMoisturePDEPE(1234.5,[0,0.7,2]);
assert(max(abs(Tq-Td))<0.01 && max(abs(Cq-Cd))<1e-4);
outOfRangeRejected=false;
try
    EvaluateDryingTemperatureMoisture(solution,18001,1);
catch
    outOfRangeRejected=true;
end
assert(outOfRangeRejected);
[initialT,initialC]=SolveDryingTemperatureMoisturePDEPE(0,[0,2]);
assert(isequal(initialT,[28,28]) && isequal(initialC,[2.55,2.55]));
[singleT,singleC]=SolveDryingTemperatureMoisturePDEPE(0.5,1);
assert(isscalar(singleT) && isscalar(singleC));
fprintf('Query, initial-state, split-time and range checks passed.\n');
report=struct('queryChecksPassed',true,'postQueryTemperatureError',max(abs(Tq-Td)), ...
    'postQueryMoistureError',max(abs(Cq-Cd)));
levels=[1,2,4];
report.refinement=zeros(numel(levels),4);
for index=1:numel(levels)
    settings=struct('StopWhenDry',true,'MeshRefinement',levels(index));
    [eventT,eventC,dryingSolution]=SolveDryingTemperatureMoisturePDEPE( ...
        0:60:144*3600,0:0.1:2,settings);
    assert(dryingSolution.thresholdReached);
    assert(abs(max(dryingSolution.thresholdMoistureDryBasis)-0.15)<1e-7);
    assert(size(eventC,1)==numel(dryingSolution.returnedTimeSeconds));
    assert(isequal(size(eventT),size(eventC)));
    assert(all(dryingSolution.returnedTimeSeconds<=dryingSolution.dryingThresholdTimeSeconds));
    report.refinement(index,:)=[levels(index),numel(dryingSolution.radiusMeters), ...
        dryingSolution.dryingThresholdTimeHours,dryingSolution.maximumRadialMoistureIncrease];
    fprintf('Mesh level %d: nodes=%d, threshold=%.9f h\n',levels(index), ...
        numel(dryingSolution.radiusMeters),dryingSolution.dryingThresholdTimeHours);
    if index==numel(levels)
        ExportDryingMoistureResultTables(dryingSolution,fullfile(rootDirectory,'validation_results'));
    end
end
% 单独收紧时间容差，分开检验空间误差和时间误差。
settings.RelativeTolerance=1e-7;
settings.AbsoluteTolerance=[1e-8;1e-10];
settings.MaximumInternalStepSeconds=150;
[~,~,tightSolution]=SolveDryingTemperatureMoisturePDEPE( ...
    [0,144*3600],0:0.1:2,settings);
assert(tightSolution.thresholdReached);
report.tightToleranceTimeHours=tightSolution.dryingThresholdTimeHours;
report.toleranceDifferenceSeconds=abs(tightSolution.dryingThresholdTimeHours- ...
    report.refinement(end,3))*3600;
% 对阈值前后各1秒重新求解；时间没有舍入到整分钟。
crossing=tightSolution.dryingThresholdTimeSeconds;
settings.StopWhenDry=false;
[~,crossingC]=SolveDryingTemperatureMoisturePDEPE([crossing-1,crossing+1],0:0.1:2,settings);
assert(max(crossingC(1,:))>0.15 && max(crossingC(2,:))<0.15);
report.crossingMoistureBeforeAfter=max(crossingC,[],2);
save(fullfile(rootDirectory,'validation_results','validation_report.mat'),'report');
fid=fopen(fullfile(rootDirectory,'validation_results','validation_report.json'),'w');
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true)); fclose(fid);
disp(report);
end
