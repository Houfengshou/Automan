function geometry = ShrinkingRadius(timeSeconds,N)
% 计算按比例径向收缩时，当前时间层的几何参数。
%
% 输入：
%   timeSeconds：时间，单位秒
%   N：径向区间数，整个计算过程中保持不变
%
% 输出结构体 geometry：
%   R ：当前外半径，单位米
%   dr：当前真实空间步长，单位米
%   r ：当前真实节点位置，(N+1)×1列向量，单位米
%   V ：径向控制体积权重，(N+1)×1列向量，单位平方米
%
% 实际控制体积 = 2*pi*药材长度*V

    validateattributes(N, {'numeric'}, ...
        {'scalar','real','finite','integer','>=',2});

    % 当前外半径
    R = RadiusMeters(timeSeconds);

    % N固定，真实步长随半径变化
    dr = R/N;

    % 节点从圆心到表面，共N+1个
    r = (0:N)'*dr;

    % 内部节点的控制体积权重
    V = r*dr;

    % 圆心控制体积：0至dr/2
    V(1) = dr^2/8;

    % 表面控制体积：R-dr/2至R
    % 与 (R^2-(R-dr/2)^2)/2 等价
    V(end) = R*dr/2 - dr^2/8;

    geometry.R = R;
    geometry.dr = dr;
    geometry.r = r;
    geometry.V = V;
end