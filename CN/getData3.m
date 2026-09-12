ProcessCNdata(0.001/648, 2);
function [T,C,t] = ProcessCNdata(dr,times)

    R = 0.02;
    N = max(1,round(R/dr));
    dr = R/N;
    dt = 1/times;

    hT = 25;
    hm = 8e-7;

    data = readmatrix("附件1去噪后.xlsx");
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
    
    nextOutputTime = 60;
    resultData = zeros(0,22); 
    while max(C) >= 0.15

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
            fprintf('计算停在 %.6f 小时，本步未收敛。\n',t/3600);
            return
        end

        T = Tnext;
        C = Cnext;
        t = tNext;
        if t == nextOutputTime
            Coutput = interp1(r,C,outputRadiusM,'linear');
            resultData(end+1,:) = [t,Coutput];
            nextOutputTime = nextOutputTime+60;
        end
    end
    if isempty(resultData) || resultData(end,1) < t
        Coutput = interp1(r,C,outputRadiusM,'linear');
        resultData(end+1,:) = [t,Coutput];
    end
    
    resultHeader = [{'时间\到药材中心的距离'}, ...
                    num2cell(outputRadiusCm)];
    
    writecell([resultHeader; num2cell(resultData)], ...
        'result3.xlsx','Sheet','Sheet1','WriteMode','overwritesheet');
    
    tableRows = find(mod(resultData(:,1),21600) == 0);
    tableRows = unique([tableRows; size(resultData,1)]);
    
    table5Data = [resultData(tableRows,1)/3600, ...
                  resultData(tableRows,2:5:22)];
    
    table5Header = [{'时间(h)\距离(cm)'}, num2cell(0:0.5:2)];
    
    writecell([table5Header; num2cell(table5Data)], ...
        'table5.xlsx','Sheet','Sheet1','WriteMode','overwritesheet');
    
    disp(array2table(table5Data,'VariableNames', ...
        {'Time_h','C_0cm','C_0p5cm','C_1cm','C_1p5cm','C_2cm'}));
    
    fprintf('结束时最大含水率：%.4f kg/kg\n',max(C));

    fprintf('干燥时间：%.4f 小时\n',t/3600);
end
