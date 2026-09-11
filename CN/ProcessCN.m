function [T,C,t] = ProcessCN(dr,times)

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

    % 空间网格与初值
    r = (0:N)'*dr;
    T = 28*ones(N+1,1);
    C = 2.55*ones(N+1,1);

    V = r*dr;
    V(1) = dr^2/8;
    V(end) = (R^2-(R-dr/2)^2)/2;

    % Picard参数
    options.maxIterations = 50;
    options.relativeTolerance = 1e-8;
    options.temperatureAbsoluteTolerance = 1e-7;
    options.moistureAbsoluteTolerance = 1e-10;

    t = 0;

    while max(C) >= 0.15

        % 两个后向欧拉半步，之后使用CN
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

        % 时间步准确落在4小时切换点
        if t < 14400
            tNext = min(tNext,14400);
        end

        % 起点环境：4小时后使用恒定环境
        [qTold,qCold] = QuerySetting( ...
            t,N+1,R,hT,hm, ...
            time_points,t_setting,c_setting,t>=14400);

        ok = false;

        % 不收敛时缩小步长，最多重试20次
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

        % 收敛后才更新状态和时间
        T = Tnext;
        C = Cnext;
        t = tNext;
    end

    fprintf('首次达到干燥条件的步末时间：%.6f 小时\n',t/3600);
end
