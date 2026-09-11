function [temperatureCelsius, moistureDryBasis, solution] = ...
    SolveDryingTemperatureMoisturePDEPE(queryTimeSeconds, queryRadiusCm, settings)
% 一维圆柱药材温湿耦合求解器（附录3）。
% [T,C,solution] = SolveDryingTemperatureMoisturePDEPE(timesSeconds,radiusCm)
% T、C：行为查询时间，列为查询半径；T单位°C，C单位kg/kg。
% 时间可以为小数；从t=0计算到查询时间的最大值。
% settings可选字段见下方默认值。默认网格仅为起始网格，须检验收敛。
% StopWhenDry=true时，阈值之后的查询行不返回，实际时间见
% solution.returnedTimeSeconds。阈值瞬间的完整场单独保存于solution。
if nargin < 2 || isempty(queryRadiusCm), queryRadiusCm = 0:0.1:2; end
if nargin < 3, settings = struct(); end
validateattributes(queryTimeSeconds,{'numeric'},{'vector','real','finite','nonempty','nonnegative'});
validateattributes(queryRadiusCm,{'numeric'},{'vector','real','finite','nonempty','>=',0,'<=',2});
queryTimeSeconds = double(queryTimeSeconds(:));
queryRadiusCm = double(queryRadiusCm(:)');
assert(all(diff(queryTimeSeconds)>0),'查询时间必须严格递增。');
assert(all(diff(queryRadiusCm)>0),'查询半径必须严格递增。');
assert(isstruct(settings) && isscalar(settings),'settings必须是标量结构体。');
rootDirectory = fileparts(mfilename('fullpath'));
defaultSettings = struct('MeshRefinement',1,'SpatialMeshMeters',[], ...
    'RelativeTolerance',1e-6,'AbsoluteTolerance',[1e-7;1e-9], ...
    'MaximumInternalStepSeconds',300,'StorageIntervalSeconds',60, ...
    'StopWhenDry',false, ...
    'EnvironmentDataFile',fullfile(rootDirectory,'source','附件1去噪后.xlsx'));
providedFields = fieldnames(settings);
assert(all(ismember(providedFields,fieldnames(defaultSettings))),'settings包含未知字段。');
for fieldIndex=1:numel(providedFields)
    defaultSettings.(providedFields{fieldIndex})=settings.(providedFields{fieldIndex});
end
settings=defaultSettings;
validateattributes(settings.MeshRefinement,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(settings.RelativeTolerance,{'numeric'},{'scalar','positive','finite'});
validateattributes(settings.AbsoluteTolerance,{'numeric'},{'vector','positive','finite'});
assert(numel(settings.AbsoluteTolerance)==1 || numel(settings.AbsoluteTolerance)==2);
validateattributes(settings.MaximumInternalStepSeconds,{'numeric'},{'scalar','positive','finite'});
validateattributes(settings.StorageIntervalSeconds,{'numeric'},{'scalar','positive','finite'});
assert(isscalar(settings.StopWhenDry) && ismember(settings.StopWhenDry,[false,true]));

materialRadiusMeters=0.02;
if isempty(settings.SpatialMeshMeters)
    % 平方映射使网格逐渐向表面加密，同时纳入题目要求的输出位置。
    computationalCoordinate=linspace(0,1,120*settings.MeshRefinement+1);
    gradedMesh=materialRadiusMeters*(1-(1-computationalCoordinate).^2);
    spatialMeshMeters=unique(round([gradedMesh,0:0.001:materialRadiusMeters],14));
else
    spatialMeshMeters=double(settings.SpatialMeshMeters(:)');
end
validateattributes(spatialMeshMeters,{'numeric'},{'vector','real','finite','nonempty'});
assert(numel(spatialMeshMeters)>=3 && all(diff(spatialMeshMeters)>0));
assert(spatialMeshMeters(1)==0 && spatialMeshMeters(end)==materialRadiusMeters, ...
    '计算网格必须从0开始、在0.02米结束。');
environment=ReadDryingEnvironmentData(settings.EnvironmentDataFile);
maximumRequestedTime=max(queryTimeSeconds);
storageTimes=(0:settings.StorageIntervalSeconds:maximumRequestedTime)';
scheduledTimes=unique([0;storageTimes;queryTimeSeconds]);
allTimes=0;
allTemperature=28*ones(1,numel(spatialMeshMeters));
allMoisture=2.55*ones(1,numel(spatialMeshMeters));
thresholdTime=NaN;
thresholdTemperature=[];
thresholdMoisture=[];
absoluteTolerance=settings.AbsoluteTolerance(:);
if numel(absoluteTolerance)==2
    % PDEPE将各节点的[T;C]展平为ODE状态，容差也要按相同顺序展开。
    absoluteTolerance=repmat(absoluteTolerance,numel(spatialMeshMeters),1);
end
odeOptions=odeset('RelTol',settings.RelativeTolerance, ...
    'AbsTol',absoluteTolerance,'MaxStep',settings.MaximumInternalStepSeconds);
if settings.StopWhenDry
    odeOptions=odeset(odeOptions,'Events',@DetectAllNodesDryingThreshold);
end
% 在4小时处重启，以明确处理环境函数切换。继承完整场，保持空间网格。
segmentBoundaries=unique([0,min(14400,maximumRequestedTime),maximumRequestedTime]);
for segmentIndex=1:numel(segmentBoundaries)-1
    segmentStart=segmentBoundaries(segmentIndex);
    segmentEnd=segmentBoundaries(segmentIndex+1);
    segmentTimes=unique([segmentStart; ...
        scheduledTimes(scheduledTimes>segmentStart & scheduledTimes<segmentEnd);segmentEnd]);
    if numel(segmentTimes)<3
        segmentTimes=unique([segmentTimes;(segmentStart+segmentEnd)/2]);
    end
    if segmentStart==0
        initialFunction=@InitialUniformTemperatureMoisture;
    else
        previousTemperature=allTemperature(end,:);
        previousMoisture=allMoisture(end,:);
        initialFunction=@(radius) ContinuePreviousTemperatureMoisture( ...
            radius,spatialMeshMeters,previousTemperature,previousMoisture);
    end
    constantEnvironment=segmentStart>=14400;
    boundaryFunction=@(xl,ul,xr,ur,t) DryingConvectionBoundaryConditions( ...
        xl,ul,xr,ur,t,environment,constantEnvironment);
    if settings.StopWhenDry
        [segmentSolution,actualTimes,eventSolution,eventTimes]=pdepe(1, ...
            @CoupledTemperatureMoistureEquations,initialFunction,boundaryFunction, ...
            spatialMeshMeters,segmentTimes,odeOptions);
    else
        segmentSolution=pdepe(1,@CoupledTemperatureMoistureEquations, ...
            initialFunction,boundaryFunction,spatialMeshMeters,segmentTimes,odeOptions);
        actualTimes=segmentTimes; eventTimes=[]; eventSolution=[];
        assert(size(segmentSolution,1)==numel(segmentTimes), ...
            'pdepe提前退出，没有计算到指定时间。请检查求解器警告。');
    end
    actualTimes=actualTimes(:);
    newRows=actualTimes>allTimes(end);
    allTimes=[allTimes;actualTimes(newRows)]; %#ok<AGROW>
    allTemperature=[allTemperature;segmentSolution(newRows,:,1)]; %#ok<AGROW>
    allMoisture=[allMoisture;segmentSolution(newRows,:,2)]; %#ok<AGROW>
    if ~isempty(eventTimes)
        thresholdTime=eventTimes(1);
        thresholdTemperature=reshape(eventSolution(1,:,1),1,[]);
        thresholdMoisture=reshape(eventSolution(1,:,2),1,[]);
        if thresholdTime>allTimes(end)
            allTimes(end+1,1)=thresholdTime;
            allTemperature(end+1,:)=thresholdTemperature;
            allMoisture(end+1,:)=thresholdMoisture;
        end
        break
    end
    assert(abs(allTimes(end)-segmentEnd)<1e-7, ...
        '求解器未完成当前时间段，不返回不完整结果。');
end
assert(all(isfinite(allTemperature(:))) && all(isfinite(allMoisture(:))), ...
    '结果包含非有限值。');
assert(all(allMoisture(:)>0),'出现非正含水率，需检查网格和容差。');
solution=struct('timeSeconds',allTimes,'radiusMeters',spatialMeshMeters, ...
    'temperatureCelsius',allTemperature,'moistureDryBasis',allMoisture, ...
    'requestedTimeSeconds',queryTimeSeconds,'queryRadiusCm',queryRadiusCm, ...
    'dryingThresholdTimeSeconds',thresholdTime,'dryingThresholdTimeHours',thresholdTime/3600, ...
    'thresholdTemperatureCelsius',thresholdTemperature, ...
    'thresholdMoistureDryBasis',thresholdMoisture,'settings',settings, ...
    'thresholdReached',isfinite(thresholdTime),'spatialConvergenceVerified',false);
availableRows=queryTimeSeconds<=allTimes(end);
solution.returnedTimeSeconds=queryTimeSeconds(availableRows);
solution.requestedTimeAvailable=availableRows;
[temperatureCelsius,moistureDryBasis]=EvaluateDryingTemperatureMoisture( ...
    solution,solution.returnedTimeSeconds,queryRadiusCm);
% 事件表示到达0.15，不宣称事件瞬间严格小于0.15。
solution.maximumRadialMoistureIncrease=max(diff(allMoisture,1,2),[],'all');
end

function environment=ReadDryingEnvironmentData(dataFile)
data=readmatrix(dataFile);
assert(size(data,2)>=3,'环境数据至少需要时间、温度、水分浓度三列。');
data=data(:,1:3);
% 仅删除完全为空的表头/空行；部分缺失的数值行必须报错。
data=data(~all(isnan(data),2),:);
assert(all(isfinite(data(:))) && all(diff(data(:,1))>0),'环境数据缺失或时间不递增。');
assert(data(1,1)==0 && data(end,1)==14400,'环境数据应覆盖0至14400秒。');
environment.temperatureInterpolant=pchip(data(:,1),data(:,2));
environment.moistureInterpolant=pchip(data(:,1),data(:,3));
end

function [ambientTemperature,ambientMoisture]=EvaluateDryingEnvironment(timeSeconds,environment,constantEnvironment)
if constantEnvironment
    ambientTemperature=50; ambientMoisture=0.05;
else
    assert(timeSeconds>=0 && timeSeconds<=14400+1e-7);
    ambientTemperature=ppval(environment.temperatureInterpolant,min(timeSeconds,14400));
    ambientMoisture=ppval(environment.moistureInterpolant,min(timeSeconds,14400));
end
end

function [capacity,flux,source]=CoupledTemperatureMoistureEquations(~,~,state,radialDerivative)
temperatureCelsius=state(1); moisture=state(2);
if moisture<=0 || temperatureCelsius<=-273.15
    error('求解过程中出现非物理状态；请检查网格和求解容差。');
end
density=650+128*moisture;
specificHeat=1450+2736*moisture/(moisture+1);
conductivity=0.21+0.38*moisture/(moisture+1);
diffusivity=2.4e-3*exp(-0.45/moisture)*exp(-3850/(temperatureCelsius+273.15));
capacity=[density*specificHeat;1];
flux=[conductivity*radialDerivative(1);diffusivity*radialDerivative(2)];
source=[0;0];
end

function initialState=InitialUniformTemperatureMoisture(~)
initialState=[28;2.55];
end

function initialState=ContinuePreviousTemperatureMoisture(radius,mesh,temperature,moisture)
initialState=[pdeval(1,mesh,temperature,radius);pdeval(1,mesh,moisture,radius)];
end

function [leftP,leftQ,rightP,rightQ]=DryingConvectionBoundaryConditions( ...
    ~,~,~,surfaceState,timeSeconds,environment,constantEnvironment)
[ambientTemperature,ambientMoisture]=EvaluateDryingEnvironment(timeSeconds,environment,constantEnvironment);
leftP=[0;0]; leftQ=[1;1]; % m=1且左端为0时由pdepe自动施加对称性。
rightP=[25*(surfaceState(1)-ambientTemperature);8e-7*(surfaceState(2)-ambientMoisture)];
rightQ=[1;1];
end

function [value,isTerminal,direction]=DetectAllNodesDryingThreshold(~,~,mesh,meshState)
% pdepe传入按节点排列的列向量：[T1 C1 T2 C2 ...]'。
stateByNode=reshape(meshState,2,numel(mesh));
value=max(stateByNode(2,:))-0.15;
isTerminal=1;
direction=-1;
end
