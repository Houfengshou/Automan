function [MT, MC, KT, KC] = AssembleDryingSpatialMatrices( ...
    T, C, r, V, dr, hT, hm)
% 根据当前温度、含水率组装径向有限体积矩阵。
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

    rho = 760 + 90.*C;
    cp = 1850 + 2150.*C./(C + 1);
    H = rho.*cp;
    k = 0.12 + 0.20.*C./(C + 1);
    D = 4.2e-4.*exp(-0.30./C) ...
        .*exp(-3850./(T + 273.15));


    rFace = (r(1:end-1) + r(2:end))/2;
    kLeft = k(1:end-1);
    kRight = k(2:end);
    kFace = 2.*kLeft.*kRight./(kLeft + kRight);

    DLeft = D(1:end-1);
    DRight = D(2:end);
    DSmall = min(DLeft, DRight);
    DLarge = max(DLeft, DRight);
    DFace = zeros(N, 1);
    validFace = DLarge > 0;
    DFace(validFace) = 2.*DSmall(validFace) ...
        ./ (1 + DSmall(validFace)./DLarge(validFace));


    gT = rFace.*kFace./dr;
    gC = rFace.*DFace./dr;

    MT = spdiags(V.*H, 0, nodeCount, nodeCount);
    MC = spdiags(V,    0, nodeCount, nodeCount);


    %       [ g  -g
    %        -g   g ]
    leftNode = (1:N)';
    rightNode = leftNode + 1;

    rowIndex = [leftNode; leftNode; rightNode; rightNode];
    colIndex = [leftNode; rightNode; leftNode; rightNode];

    heatValues = [gT; -gT; -gT; gT];
    moistureValues = [gC; -gC; -gC; gC];

    KT = sparse(rowIndex, colIndex, heatValues, ...
                nodeCount, nodeCount);

    KC = sparse(rowIndex, colIndex, moistureValues, ...
                nodeCount, nodeCount);

    KT(end,end) = KT(end,end) + R*hT;
    KC(end,end) = KC(end,end) + R*hm;
end