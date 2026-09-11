function ProcessCN(dr, times)
    R = 0.02;
    N = max(1,round(R/dr)); % 区间数; 节点数 = 区间数 + 1
    dr = R/N;
    dt = 1/times;

    hT = 25; %对流换热
    hm = 8e-7; %对流传质

    data = readmatrix("附件1去噪后.xlsx");
    time_points = data(:,1).';
    t_setting = data(:,2).';
    c_setting = data(:,3).';

    r = (0:N)'*dr;
    T = 28.*ones(N + 1,1);
    C = 2.55.*ones(N + 1,1);
    V = zeros(N + 1,1);
    V(1,1) = dr^2/8;
    for index = 2:N
        V(index,1) = (index-1)*dr^2;
    end
    V(N+1,1) = (R^2 - (R-dr/2)^2)/2;
    [MT,MC,KT,KC] = AssembleDryingSpatialMatrices(T, C, r, V, dr, hT, hm);
 
end 
