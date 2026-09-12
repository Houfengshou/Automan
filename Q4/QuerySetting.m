function [qT, qC, Ta, Ca] = QuerySetting( ...
    t, nodeCount, R, hT, hm, ...
    time_points, t_setting, c_setting, ...
    usePostSwitchAtBoundary)

    if nargin < 9
        usePostSwitchAtBoundary = false;
    end

    assert(isscalar(t) && isfinite(t) && t >= 0, ...
        '查询时间必须是有限的非负标量。');

    switchTime = 14400;
    Tinf = 50; %最终环境温度恒定值
    Cinf = 0.05; %最终环境湿度恒定值

    if t > switchTime || ...
       (t == switchTime && usePostSwitchAtBoundary)

        Ta = Tinf;
        Ca = Cinf;

    else
        Ta = interp1(time_points, t_setting, t, 'pchip');
        Ca = interp1(time_points, c_setting, t, 'pchip');
    end

    assert(isfinite(Ta) && isfinite(Ca), ...
        '环境查询失败，请检查附件数据及其时间覆盖范围。');

    qT = zeros(nodeCount, 1);
    qC = zeros(nodeCount, 1);

    qT(end) = R*hT*Ta;
    qC(end) = R*hm*Ca;
end