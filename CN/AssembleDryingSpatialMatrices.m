function [MT, MC, KT, KC] = AssembleDryingSpatialMatrices( ...
    T, C, r, V, dr, hT, hm)
% 根据当前温度、含水率组装径向有限体积矩阵。
%
% 输入：
%   T、C：节点温度（摄氏度）、干基含水率，均为列向量
%   r、V：节点半径、径向控制体积权重，均为列向量
%   dr：空间步长
%   hT、hm：表面换热系数、传质系数
%
% 输出：
%   MT、MC：温度、水分质量矩阵
%   KT、KC：温度、水分传输矩阵
%
% 对应方程：
%   MT * dT/dt + KT * T = qT
%   MC * dC/dt + KC * C = qC

    % 统一为列向量
    T = T(:);
    C = C(:);
    r = r(:);
    V = V(:);

    nodeCount = numel(r);
    N = nodeCount - 1;       % 区间数
    R = r(end);

    assert(nodeCount >= 2, '至少需要两个空间节点。');
    assert(numel(T) == nodeCount && ...
           numel(C) == nodeCount && ...
           numel(V) == nodeCount, ...
           'T、C、r、V的长度必须一致。');

    assert(all(isfinite(T)) && all(T > -273.15), ...
           '温度必须有限且高于绝对零度。');
    assert(all(isfinite(C)) && all(C > 0), ...
           '含水率必须为有限正数。');
    assert(isfinite(dr) && dr > 0 && ...
           all(isfinite(V)) && all(V > 0), ...
           '空间步长和体积权重必须为有限正数。');

    % 1. 节点物性
    rho = 650 + 128.*C;
    cp = 1450 + 2736.*C./(C + 1);
    H = rho.*cp;

    k = 0.21 + 0.38.*C./(C + 1);

    D = 0.0024.*exp(-0.45./C) ...
        .*exp(-3850./(T + 273.15));

    % 2. 界面位置与调和平均物性
    rFace = (r(1:end-1) + r(2:end))/2;

    kLeft = k(1:end-1);
    kRight = k(2:end);
    kFace = 2.*kLeft.*kRight./(kLeft + kRight);

    DLeft = D(1:end-1);
    DRight = D(2:end);

    % 等价的调和平均写法，避免直接计算两个很小的D之积
    DSmall = min(DLeft, DRight);
    DLarge = max(DLeft, DRight);

    DFace = zeros(N, 1);
    validFace = DLarge > 0;
    DFace(validFace) = 2.*DSmall(validFace) ...
        ./ (1 + DSmall(validFace)./DLarge(validFace));

    % 3. 界面传输系数
    gT = rFace.*kFace./dr;
    gC = rFace.*DFace./dr;

    % 4. 质量矩阵
    MT = spdiags(V.*H, 0, nodeCount, nodeCount);
    MC = spdiags(V,    0, nodeCount, nodeCount);

    % 5. 每个界面对两个相邻节点贡献一个2×2矩阵：
    %       [ g  -g
    %        -g   g ]
    leftNode = (1:N)';
    rightNode = leftNode + 1;

    rowIndex = [leftNode; leftNode; rightNode; rightNode];
    colIndex = [leftNode; rightNode; leftNode; rightNode];

    heatValues = [gT; -gT; -gT; gT];
    moistureValues = [gC; -gC; -gC; gC];

    % sparse会自动累加重复位置的元素
    KT = sparse(rowIndex, colIndex, heatValues, ...
                nodeCount, nodeCount);

    KC = sparse(rowIndex, colIndex, moistureValues, ...
                nodeCount, nodeCount);

    % 6. 表面交换项加入最后一个对角元素
    KT(end,end) = KT(end,end) + R*hT;
    KC(end,end) = KC(end,end) + R*hm;
end