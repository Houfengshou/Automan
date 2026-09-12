dr = 0.001/81;
times = 2;
rootDirectory = fileparts(mfilename('fullpath'));
[T,C,t,resultData,surfaceRadiusCm,table6Data,summary] = ...
    ExportShrinkingResult4(dr,times,rootDirectory);

function [T,C,t,resultData,surfaceRadiusCm,table6Data,summary] = ...
    ExportShrinkingResult4(dr,times,rootDirectory)

    validateattributes(dr, {'numeric'}, ...
        {'scalar','real','finite','positive'});

    validateattributes(times, {'numeric'}, ...
        {'scalar','real','finite','positive'});


    R0 = RadiusMeters(0);

    N = max(2,round(R0/dr));

    dt = 1/times;
    assert(isfinite(dt) && dt > 0, '时间步长必须为有限正数。');

    hT = 25;
    hm = 8e-7;


    templateFile = fullfile(fileparts(rootDirectory),'附件','附件3','result4.xlsx');
    outputFile = fullfile(rootDirectory,'result4.xlsx');
    tableFile = fullfile(rootDirectory,'table6.xlsx');
    assert(isfile(templateFile),'Template not found: %s',templateFile);
    assert(~isfile(outputFile) && ~isfile(tableFile), ...
        'Output exists. Rename the existing result4.xlsx/table6.xlsx before rerunning.');
    templateSheets = sheetnames(templateFile);
    assert(numel(templateSheets)==1 && templateSheets(1)=="Sheet1", ...
        'Expected the original result4.xlsx Sheet1 template.');
    dataFile = fullfile(rootDirectory,'附件1去噪后.xlsx');

    data = readmatrix(dataFile);

    assert(size(data,2) >= 3, ...
        '环境数据至少需要时间、温度和水分浓度三列。');

    data = data(:,1:3);
    data = data(~all(isnan(data),2),:);

    assert(size(data,1) >= 2 && all(isfinite(data(:))), ...
        '环境数据存在缺失或非有限值。');

    time_points = data(:,1);
    t_setting = data(:,2);
    c_setting = data(:,3);

    assert(all(diff(time_points) > 0), ...
        '环境数据时间必须严格递增。');

    assert(time_points(1) == 0 && time_points(end) == 14400, ...
        '环境数据应覆盖0至14400秒。');


    T = 28*ones(N+1,1);
    C = 2.55*ones(N+1,1);


    options.maxIterations = 50;
    options.relativeTolerance = 1e-8;
    options.temperatureAbsoluteTolerance = 1e-7;
    options.moistureAbsoluteTolerance = 1e-10;

    t = 0;
    outputRadiusCm = (0:19)/10;
    outputRadiusM = outputRadiusCm/100;
    nextOutputTime = 60;
    resultData = zeros(0,22); 
    surfaceRadiusCm = zeros(0,1);
    acceptedSteps = 0;
    retriedSteps = 0;
    minimumAcceptedStep = inf;

    while max(C) >= 0.15

        geometryOld = ShrinkingRadius(t,N);

        if t < dt/2
            theta = 1;
            tNext = dt/2;
        elseif t < dt
            theta = 1;
            tNext = dt;
        else
            theta = 0.5;
            tNext = t+dt;
        end

        if t < 14400
            tNext = min(tNext,14400);
        end

        tNext = min(tNext,nextOutputTime);


        [qTold,qCold] = QuerySetting( ...
            t,N+1,geometryOld.R,hT,hm, ...
            time_points,t_setting,c_setting,t>=14400);

        ok = false;

        for retry = 0:20

            tau = tNext-t;

            if tau <= 0
                break
            end

            geometryNext = ShrinkingRadius(tNext,N);

            [qTnext,qCnext] = QuerySetting( ...
                tNext,N+1,geometryNext.R,hT,hm, ...
                time_points,t_setting,c_setting);

            [Tnext,Cnext,ok] = Picard( ...
                T,C,geometryOld,geometryNext,hT,hm, ...
                tau,theta,qTold,qCold,qTnext,qCnext,options);

            if ok
                break
            end

            tNext = t+tau/2;
        end

        if ~ok
            error('Picard失败', '计算停在 %.4f h', t/3600);
        end

        previousTime = t;
        previousMaximumMoisture = max(C);
        acceptedSteps = acceptedSteps+1;
        retriedSteps = retriedSteps+(retry>0);
        minimumAcceptedStep = min(minimumAcceptedStep,tau);

        T = Tnext;
        C = Cnext;
        t = tNext;
        if t == nextOutputTime
            resultData(end+1,:) = MakeShrinkingOutputRow(t,C,geometryNext,outputRadiusM);
            surfaceRadiusCm(end+1,1) = geometryNext.R*100;
            nextOutputTime = nextOutputTime+60;
            if mod(t,21600)==0
                fprintf('t=%.0f h, R=%.6f cm, max C=%.8f\n', ...
                    t/3600,geometryNext.R*100,max(C));
            end
        end
    end

    if isempty(resultData) || resultData(end,1)<t
        resultData(end+1,:) = MakeShrinkingOutputRow(t,C,geometryNext,outputRadiusM);
        surfaceRadiusCm(end+1,1) = geometryNext.R*100;
    end
    assert(max(C)<0.15 && all(isfinite(resultData(:,[1,2,end])),'all'), 'Incomplete drying result.');

    resultHeader = [{'时间\到药材中心的距离'},num2cell(outputRadiusCm),{'药材表面'}];
    [okCopy,copyMessage] = copyfile(templateFile,outputFile);
    assert(okCopy,'Template copy failed: %s',copyMessage);
    WriteShrinkingOutputSheet(resultHeader,resultData,outputFile,'Sheet1');

    tableRows = find(mod(resultData(:,1),21600)==0);
    tableRows = unique([tableRows;size(resultData,1)]);
    table6Data = [resultData(tableRows,1)/3600, ...
        resultData(tableRows,[2,7,12,17,end])];
    tableHeader = [{'时间(h)\距离(cm)'},num2cell([0,0.5,1,1.5]),{'药材表面'}];
    WriteShrinkingOutputSheet(tableHeader,table6Data,tableFile,'Sheet1');

    summary = struct('initialStepMeters',R0/N,'radialIntervals',N, ...
        'nominalStepSeconds',dt,'minimumAcceptedStepSeconds',minimumAcceptedStep, ...
        'acceptedSteps',acceptedSteps,'stepsWithRetry',retriedSteps, ...
        'dryingTimeBracketSeconds',[previousTime,t], ...
        'previousMaximumMoisture',previousMaximumMoisture, ...
        'finalMaximumMoisture',max(C),'firstDryStepHours',t/3600);
    fprintf('干燥时间：%.6f 小时\n',t/3600);
    fprintf('Exported: %s\n',outputFile);
    fprintf('Exported: %s\n',tableFile);
end
function outputRow = MakeShrinkingOutputRow(t,C,geometry,outputRadiusM)
    inside = outputRadiusM<=geometry.R;
    values = nan(size(outputRadiusM));
    values(inside) = interp1(geometry.r,C,outputRadiusM(inside),'linear');
    assert(all(isfinite(values(inside))) && isfinite(C(end)), ...
        'Interpolation failed at an interior radius.');
    outputRow = [t,values,C(end)];
end

function WriteShrinkingOutputSheet(header,values,fileName,sheetName)
    roundedValues = values;
    roundedValues(:,2:end) = round(values(:,2:end),4);
    body = num2cell(roundedValues);
    body(isnan(roundedValues)) = {''};
    writecell(header,fileName,'Sheet',sheetName,'Range','A1');
    writecell(body,fileName,'Sheet',sheetName,'Range','A2');
end