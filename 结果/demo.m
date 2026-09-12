%% 两组收敛数据
y1 = [
148.030694
114.238333
90.831667
76.678611
68.109861
63.425694
60.875278
59.443750
58.648472
58.183333
57.913472
57.750833
57.650556
57.588056
57.548194
57.523056
57.506528
57.495556
57.488333
57.483333
57.480000
57.477778
57.476111
57.475000
57.474167
57.473750
57.473333
57.473056
57.472778
57.472778
];

y2 = [
52.688750
52.106389
51.763611
51.539167
51.385417
51.279583
51.211389
51.164444
51.132917
51.110972
51.095972
51.085556
51.078333
51.073194
51.069583
51.067083
51.065417
51.064167
51.063194
51.062639
51.062222
51.061806
51.061667
51.061528
51.061389
51.061250
51.061250
51.061250
51.061250
51.061111
];

x = exp(linspace(1,6,30));

%% ========== 图1：独立大图 方案1 ==========
figure(1);
set(gcf,'Color','w','Position',[50,50,850,600]);
plot(x,y1,'o-','LineWidth',1.3,'MarkerSize',5,'Color','#0072BD');
xlabel('exp(linspace(1,6,30))');
ylabel('首次达到干燥条件的步末时间 / h');
title('方案1：干燥终止时间收敛曲线');
grid on; grid minor;
% set(gca,'XScale','log');
savefig('fig1_model1.fig');
print(gcf,'fig1_model1','-deps','-r300');

%% ========== 图2：独立大图 方案2 ==========
figure(2);
set(gcf,'Color','w','Position',[120,120,850,600]);
plot(x,y2,'s-','LineWidth',1.3,'MarkerSize',5,'Color','#D95319');
xlabel('exp(linspace(1,6,30))');
ylabel('首次达到干燥条件的步末时间 / h');
title('方案2：干燥终止时间收敛曲线');
grid on; grid minor;
% set(gca,'XScale','log');
savefig('fig2_model2.fig');
print(gcf,'fig2_model2','-deps','-r300');

%% ========== 图3：同一坐标系两条曲线叠加 ==========
figure(3);
set(gcf,'Color','w','Position',[190,190,850,600]);
hold on;
plot(x,y1,'o-','LineWidth',1.3,'MarkerSize',5,'Color','#0072BD','DisplayName','方案1');
plot(x,y2,'s-','LineWidth',1.3,'MarkerSize',5,'Color','#D95319','DisplayName','方案2');
hold off;
xlabel('exp(linspace(1,6,30))');
ylabel('首次达到干燥条件的步末时间 / h');
title('两种方案干燥终止时间收敛对比（同一坐标系）');
legend('Location','best');
grid on; grid minor;
% set(gca,'XScale','log');
savefig('fig3_compare_both.fig');
print(gcf,'fig3_compare_both','-deps','-r300');
