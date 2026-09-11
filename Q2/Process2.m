function Process2(dr, times)
    radius = 0.02;
    span = 432000;
    M = floor(0.02/dr) + 1;
    N = times * span + 1; %附件一总时长
    dt = 1/times;

    hT = 25; %对流换热
    hm = 8e-7; %对流传质
    Tinf = 50;
    Cinf = 0.05;

    data = readmatrix("附件1去噪后.xlsx");
    time_points = data(:,1).';
    t_setting = data(:,2).';
    c_setting = data(:,3).';
    TimePoints = 0:1/times:14400;
    TSettings = interp1(time_points, t_setting, TimePoints, "pchip");
    CSettings = interp1(time_points, c_setting, TimePoints, "pchip");

    % ----------------------------------------------
    rskinned = radius - dr/2;
    VNstar =  (radius^2 - rskinned^2)/2;
    
    TResult = zeros(M, N);
    TResult(:,1) = 28;
    CResult = zeros(M, N);
    CResult(:,1) = 2.55;

    for n = 2:N
       qT = (4*dt*HarMean(k2(CResult(1,n-1)), k2(CResult(2,n-1))))/(H2(CResult(1,n-1))*dr^2);
       qC = (4*dt*HarMean(D2(TResult(1,n-1), CResult(1,n-1)), D2(TResult(2,n-1), CResult(2,n-1))))/(dr^2);
       TResult(1,n) = (1 - qT)*TResult(1,n-1) + qT*TResult(2,n-1);
       CResult(1,n) = (1 - qC)*CResult(1,n-1) + qC*CResult(2,n-1);
       for m = 2:M-1
           AT = (dt*((m-1)-1/2)*HarMean(k2(CResult(m-1,n-1)), k2(CResult(m,n-1))))/(H2(CResult(m,n-1))*(m-1)*dr^2);
           BT = (dt*((m-1)+1/2)*HarMean(k2(CResult(m,n-1)), k2(CResult(m+1,n-1))))/(H2(CResult(m,n-1))*(m-1)*dr^2);    
           AC = (dt*((m-1)-1/2)*HarMean(D2(TResult(m-1,n-1), CResult(m-1,n-1)), D2(TResult(m,n-1), CResult(m,n-1))))/((m-1)*dr^2);
           BC = (dt*((m-1)+1/2)*HarMean(D2(TResult(m,n-1), CResult(m,n-1)), D2(TResult(m+1,n-1), CResult(m+1,n-1))))/((m-1)*dr^2);
           TResult(m,n) = AT*TResult(m-1,n-1) + (1 - AT - BT)*TResult(m,n-1) + BT*TResult(m+1,n-1);
           CResult(m,n) = AC*CResult(m-1,n-1) + (1 - AC - BC)*CResult(m,n-1) + BC*CResult(m+1,n-1);
       end
       aT = (dt*radius*hT)/(H2(CResult(M,n-1))*VNstar);
       bT = (dt*rskinned*HarMean(k2(CResult(M-1,n-1)), k2(CResult(M,n-1))))/(H2(CResult(M,n-1))*VNstar*dr);   
       aC = dt*radius*hm/VNstar;
       bC = (dt*rskinned*HarMean(D2(TResult(M-1,n-1), CResult(M-1,n-1)), D2(TResult(M,n-1), CResult(M,n-1))))/(VNstar*dr);
       if n > length(TimePoints)
           TSetting = Tinf;
           CSetting = Cinf;
       else
           TSetting = TSettings(n-1);
           CSetting = CSettings(n-1);   
       end
       TResult(M,n) = (1 - aT - bT)*TResult(M,n-1) + bT*TResult(M-1,n-1) + aT*TSetting;
       CResult(M,n) = (1 - aC - bC)*CResult(M,n-1) + bC*CResult(M-1,n-1) + aC*CSetting;
    end

    T = TResult(floor([0,0.5,1,1.5,2]*0.01/dr+1),floor([0.5,1,1.5,2,2.5,3]*60*60*times + 1)).'
    C = CResult(floor([0,0.5,1,1.5,2]*0.01/dr+1),floor([0.5,1,1.5,2,2.5,3]*60*60*times + 1)).'

end