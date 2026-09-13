% Run this script: Q2 and Q3 share the first three hours of evolution.
dr = 0.001/1;
times = 1;
rootDirectory = fileparts(mfilename('fullpath'));
addpath(rootDirectory);
summary = ExportExplicitResult2AndResult3(dr,times,fullfile(rootDirectory,'results'));


function summary = ExportExplicitResult2AndResult3(dr, times, outDir)
% Q2/Q3: rolling time layers, explicit equations from Q3/Process2.m.
% Example: ExportExplicitResult2AndResult3(0.001/9,81)
% Temperature is stored in degC; moisture is kg/kg (dry basis).
if nargin < 1, dr = 0.001/9; end
if nargin < 2, times = 81; end
root = fileparts(mfilename('fullpath'));
if nargin < 3, outDir = fullfile(root,'results'); end
if ~exist(outDir,'dir'), mkdir(outDir); end
assert(~isfile(fullfile(outDir,'result2.xlsx')) && ~isfile(fullfile(outDir,'result3.xlsx')), 'Output files exist; choose another output directory.');
R = 0.02; dt = 1/times; Nr = round(R/dr); M = Nr+1;
assert(times>0 && times==round(times),'times must be a positive integer.');
assert(Nr>=2 && abs(Nr*dr-R)<1e-12,'dr must divide R.');
outR = (0:20)'*0.001;
idx = round(outR/dr)+1;
assert(max(abs((idx-1)*dr-outR))<1e-12,'Output radii must be mesh nodes.');
data = readmatrix(fullfile(root,'附件1去噪后.xlsx'));
data = data(all(isfinite(data(:,1:3)),2),1:3);
assert(data(1,1)==0 && data(end,1)==14400 && all(diff(data(:,1))>0));
envTime = (0:14400*times)'*dt;
Ta = interp1(data(:,1),data(:,2),envTime,'pchip');
Ca = interp1(data(:,1),data(:,3),envTime,'pchip');
temperatureData=zeros(10800,22); moistureData=zeros(10800,22);
T = 28*ones(M,1); C = 2.55*ones(M,1);
hT = 25; hm = 8e-7;
i = (1:Nr-1)'; rs = R-dr/2; Vs = (R^2-rs^2)/2;
left = dt*(i-0.5)./(i*dr^2);
right = dt*(i+0.5)./(i*dr^2);
saveEvery = 60*times; sixHours = 21600*times;
results = zeros(10000,22); table5 = zeros(200,6); ns=0; nt=0;
worstWeight=0; maxRadialIncrease=0; maxBalanceError=0;
r=(0:Nr)'*dr; edges=[0; (r(1:end-1)+r(2:end))/2; R];
vol=diff(edges.^2)/2;
timer=tic;
n=0;
while true
    n=n+1;
    Told=T; Cold=C;
    rho=650+128*Cold;
    cp=1450+2736*Cold./(Cold+1);
    H=rho.*cp;
    k=0.21+0.38*Cold./(Cold+1);
    D=2.4e-3*exp(-0.45./Cold).*exp(-3850./(Told+273.15));
    kf=2*k(1:end-1).*k(2:end)./(k(1:end-1)+k(2:end));
    Df=2*D(1:end-1).*D(2:end)./(D(1:end-1)+D(2:end));
    if n <= numel(Ta)
        ambientT=Ta(n); ambientC=Ca(n);
    else
        ambientT=50; ambientC=0.05;
    end
    qT=4*dt*kf(1)/(H(1)*dr^2); qC=4*dt*Df(1)/dr^2;
    AT=left.*kf(1:end-1)./H(2:end-1);
    BT=right.*kf(2:end)./H(2:end-1);
    AC=left.*Df(1:end-1); BC=right.*Df(2:end);
    aT=dt*R*hT/(H(end)*Vs); bT=dt*rs*kf(end)/(H(end)*Vs*dr);
    aC=dt*R*hm/Vs; bC=dt*rs*Df(end)/(Vs*dr);
    weight=max([qT;qC;AT+BT;AC+BC;aT+bT;aC+bC]);
    assert(weight<=1+1e-12,'Explicit scheme unstable: increase times.');
    worstWeight=max(worstWeight,weight);
    T(1)=Told(1)+qT*(Told(2)-Told(1));
    C(1)=Cold(1)+qC*(Cold(2)-Cold(1));
    T(2:end-1)=Told(2:end-1)+AT.*(Told(1:end-2)-Told(2:end-1))+BT.*(Told(3:end)-Told(2:end-1));
    C(2:end-1)=Cold(2:end-1)+AC.*(Cold(1:end-2)-Cold(2:end-1))+BC.*(Cold(3:end)-Cold(2:end-1));
    T(end)=Told(end)+bT*(Told(end-1)-Told(end))+aT*(ambientT-Told(end));
    C(end)=Cold(end)+bC*(Cold(end-1)-Cold(end))+aC*(ambientC-Cold(end));
    assert(all(isfinite(T)) && all(isfinite(C)) && all(C>0),'Nonphysical solution.');
    maxRadialIncrease=max(maxRadialIncrease,max(diff(C)));
    balance=abs(sum(vol.*(C-Cold))+dt*R*hm*(Cold(end)-ambientC));
    maxBalanceError=max(maxBalanceError,balance);
    if n<=10800*times && mod(n,times)==0
        outputSecond=n/times;
        temperatureData(outputSecond,:)=[outputSecond,T(idx)'];
        moistureData(outputSecond,:)=[outputSecond,C(idx)'];
        if outputSecond==10800
            WriteExplicitResult2(root,outDir,temperatureData,moistureData);
        end
    end
    if mod(n,saveEvery)==0
        ns=ns+1; results(ns,:)=[n*dt,C(idx)'];
    end
    if mod(n,sixHours)==0
        nt=nt+1; table5(nt,:)=[n*dt/3600,C(idx(1:5:21))'];
        fprintf('t=%.0f h, max C=%.8f, wall=%.1f s\n',n*dt/3600,max(C),toc(timer));
    end
    if max(C)<0.15, break; end
end
assert(n>=10800*times,'Drying ended before Q2 output finished.');
stopSeconds=n*dt;
if mod(n,saveEvery)~=0
    ns=ns+1; results(ns,:)=[stopSeconds,C(idx)'];
end
if mod(n,sixHours)~=0
    nt=nt+1; table5(nt,:)=[stopSeconds/3600,C(idx(1:5:21))'];
end
results=results(1:ns,:); table5=table5(1:nt,:);
summary=struct('dr_m',dr,'dt_s',dt,'stop_s',stopSeconds, ...
    'stop_h',stopSeconds/3600,'previous_s',(n-1)*dt, ...
    'previous_max_C',max(Cold),'final_max_C',max(C), ...
    'max_weight_sum',worstWeight,'max_radial_increase',maxRadialIncrease, ...
    'max_discrete_balance_error',maxBalanceError,'wall_s',toc(timer));
save(fullfile(outDir,'q3_full_precision.mat'),'summary','results','table5','T','C','r','temperatureData','moistureData');
fid=fopen(fullfile(outDir,'summary.json'),'w'); fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true)); fclose(fid);
template=readcell(fullfile(root,'templates','result3.xlsx'),'Sheet','Sheet1','Range','A1:V1');
header=template;
writecell(header,fullfile(outDir,'result3.xlsx'),'Sheet','Sheet1','Range','A1');
writematrix([results(:,1),round(results(:,2:end),4)],fullfile(outDir,'result3.xlsx'),'Sheet','Sheet1','Range','A2');
writecell([{'时间/h'},num2cell(0:0.5:2)],fullfile(outDir,'table5.xlsx'),'Range','A1');
writematrix(round(table5,4),fullfile(outDir,'table5.xlsx'),'Range','A2');
disp(summary); disp(table5);
end

function WriteExplicitResult2(root,outDir,temperatureData,moistureData)
file=fullfile(outDir,'result2.xlsx');
for sheet=["温度","水分浓度"]
    header=readcell(fullfile(root,'templates','result2.xlsx'),'Sheet',sheet,'Range','A1:V1');
    if sheet=="温度", values=temperatureData; else, values=moistureData; end
    writecell(header,file,'Sheet',sheet,'Range','A1');
    writematrix([values(:,1),round(values(:,2:end),4)],file,'Sheet',sheet,'Range','A2');
end
fprintf('result2.xlsx exported at 3 hours.\n');
end
