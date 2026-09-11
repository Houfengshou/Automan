function R = RadiusMeters(timeSeconds)
% 根据附件2的双指数拟合，计算当前药材外半径。
%
% 输入：
%   timeSeconds：时间，单位秒，非负标量
%
% 输出：
%   R：药材外半径，单位米

    validateattributes(timeSeconds, {'numeric'}, ...
        {'scalar','real','finite','nonnegative'});

    % 秒转换为小时
    timeHours = timeSeconds/3600;

    % 拟合公式中的半径单位为厘米
    radiusCm = 1.1998 + 0.8002*( ...
        0.2082*exp(-timeHours/1.1588) ...
        + 0.7918*exp(-timeHours/4.6596));

    % 厘米转换为米
    R = radiusCm/100;
end