function [T,C,t] = Process4(dr,times)
% 第四问：在原ProcessCN基础上加入按比例径向收缩。
%
% dr：初始空间步长，单位米
% times：时间细化参数，正常时间步长dt=1/times秒
%
% T、C：首次达标步末的温度、含水率列向量
% t：首次达标步末时间，单位秒

    validateattributes(dr, {'numeric'}, ...
        {'scalar','real','finite','positive'});

    validateattributes(times, {'numeric'}, ...
        {'scalar','real','finite','positive'});

    %% 初始几何与时间步长

    R0 = RadiusMeters(0);

    % 区间数只在初始时确定，之后保持不变
    N = max(2,round(R0/dr));

    dt = 1/times;
    assert(isfinite(dt) && dt > 0, '时间步长必须为有限正数。');

    hT = 25;
    hm = 8e-7;

    %% 环境数据

    rootDirectory = fileparts(mfilename('fullpath'));
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

    %% 初始状态

    T = 28*ones(N+1,1);
    C = 2.55*ones(N+1,1);

    %% Picard参数：沿用原CN程序

    options.maxIterations = 50;
    options.relativeTolerance = 1e-8;
    options.temperatureAbsoluteTolerance = 1e-7;
    options.moistureAbsoluteTolerance = 1e-10;

    t = 0;

    while max(C) >= 0.15

        % 旧时间层几何
        geometryOld = ShrinkingRadius(t,N);

        % 沿用原程序：两个后向欧拉半步，之后使用CN
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

        % 时间步准确落在4小时环境切换点
        if t < 14400
            tNext = min(tNext,14400);
        end

        % 旧时间层边界使用旧半径
        [qTold,qCold] = QuerySetting( ...
            t,N+1,geometryOld.R,hT,hm, ...
            time_points,t_setting,c_setting,t>=14400);

        ok = false;

        % 首次尝试，加最多20次缩步重试
        for retry = 0:20

            tau = tNext-t;

            if tau <= 0
                break
            end

            % 每次重试都重新计算新时间层几何
            geometryNext = ShrinkingRadius(tNext,N);

            % 新时间层边界使用新半径
            [qTnext,qCnext] = QuerySetting( ...
                tNext,N+1,geometryNext.R,hT,hm, ...
                time_points,t_setting,c_setting);

            [Tnext,Cnext,ok] = Picard( ...
                T,C,geometryOld,geometryNext,hT,hm, ...
                tau,theta,qTold,qCold,qTnext,qCnext,options);

            if ok
                break
            end

            % 本步未收敛：从原来的T、C重新尝试半步
            tNext = t+tau/2;
        end

        if ~ok
            error('Process4:StepFailed', ...
                '计算停在 %.6f 小时，本步在缩步重试后仍未收敛。', ...
                t/3600);
        end

        % 收敛后才更新状态和时间
        T = Tnext;
        C = Cnext;
        t = tNext;
    end

    fprintf('首次达到干燥条件的步末时间：%.6f 小时\n',t/3600);
end