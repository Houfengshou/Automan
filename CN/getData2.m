dr = 0.001/81;
times = 2;
rootDirectory = fileparts(mfilename('fullpath'));
[T,C,t,temperatureData,moistureData] = ExportCNResult2(dr,times,rootDirectory);

function [T,C,t,temperatureData,moistureData] = ExportCNResult2(dr,times,rootDirectory)
    validateattributes(dr,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(times,{'numeric'},{'scalar','real','finite','positive'});
    endTimeSeconds = 10800;
    templateFile = fullfile(fileparts(rootDirectory),'附件','附件3','result2.xlsx');
    outputFile = fullfile(rootDirectory,'result2.xlsx');
    assert(isfile(templateFile),'Template not found: %s',templateFile);
    assert(~isfile(outputFile),'Output exists; rename it or change outputFile: %s',outputFile);
    expectedSheets = ["温度";"水分浓度"];
    actualSheets = sheetnames(templateFile);
    assert(numel(actualSheets)==2 && all(ismember(expectedSheets,actualSheets)), ...
        'Expected two template sheets: 温度 and 水分浓度.');

    R = 0.02;
    N = max(1,round(R/dr));
    dr = R/N;
    dt = 1/times;
    assert(isfinite(dt) && dt>0,'Invalid time step.');

    hT = 25;
    hm = 8e-7;

    data = readmatrix(fullfile(rootDirectory,'附件1去噪后.xlsx'));
    assert(size(data,2)>=3,'Environment data needs three columns.');
    data = data(:,1:3);
    data = data(~all(isnan(data),2),:);
    assert(size(data,1)>=2 && all(isfinite(data(:))),'Invalid environment data.');
    assert(all(diff(data(:,1))>0) && data(1,1)==0 && data(end,1)>=endTimeSeconds, ...
        'Invalid environmental time coverage.');
    time_points = data(:,1);
    t_setting = data(:,2);
    c_setting = data(:,3);

    r = (0:N)'*dr;
    T = 28*ones(N+1,1);
    C = 2.55*ones(N+1,1);

    V = r*dr;
    V(1) = dr^2/8;
    V(end) = (R^2-(R-dr/2)^2)/2;


    options.maxIterations = 50;
    options.relativeTolerance = 1e-8;
    options.temperatureAbsoluteTolerance = 1e-7;
    options.moistureAbsoluteTolerance = 1e-10;

    t = 0;

    outputRadiusCm = 0:0.1:2;
    outputRadiusM = outputRadiusCm/100;
    
    nextOutputTime = 1;
    temperatureData = zeros(endTimeSeconds,22);
    moistureData = zeros(endTimeSeconds,22);
    outputRow = 0;
    while t < endTimeSeconds

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

        tNext = min([tNext,nextOutputTime,endTimeSeconds]);
        [qTold,qCold] = QuerySetting( ...
            t,N+1,R,hT,hm, ...
            time_points,t_setting,c_setting,t>=14400);

        ok = false;

        for retry = 0:20

            tau = tNext-t;

            if tau <= 0
                break
            end

            [qTnext,qCnext] = QuerySetting( ...
                tNext,N+1,R,hT,hm, ...
                time_points,t_setting,c_setting);

            [Tnext,Cnext,ok] = Picard( ...
                T,C,r,V,dr,hT,hm,tau,theta, ...
                qTold,qCold,qTnext,qCnext,options);

            if ok
                break
            end

            tNext = t+tau/2;
        end

        if ~ok
            error('ExportCNResult2:StepFailed', ...
                '计算停在 %.6f 小时，本步未收敛；未导出不完整结果。\n',t/3600);
        end

        T = Tnext;
        C = Cnext;
        t = tNext;
        if t == nextOutputTime
            outputRow = outputRow+1;
            Toutput = interp1(r,T,outputRadiusM,'linear');
            Coutput = interp1(r,C,outputRadiusM,'linear');
            temperatureData(outputRow,:) = [t,Toutput];
            moistureData(outputRow,:) = [t,Coutput];
            nextOutputTime = nextOutputTime+1;
            if mod(outputRow,1800)==0
                fprintf('已完成 %.1f h / 3 h\n',t/3600);
            end
        end
    end

    resultHeader = [{'时间\到药材中心的距离'},num2cell(outputRadiusCm)];
    writecell(resultHeader,outputFile,'Sheet','温度','Range','A1');
    writematrix([temperatureData(:,1),round(temperatureData(:,2:end),4)], ...
        outputFile,'Sheet','温度','Range','A2');
    writecell(resultHeader,outputFile,'Sheet','水分浓度','Range','A1');
    writematrix([moistureData(:,1),round(moistureData(:,2:end),4)], ...
        outputFile,'Sheet','水分浓度','Range','A2');
end