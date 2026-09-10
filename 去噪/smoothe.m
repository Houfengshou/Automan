clc; clear; close all;


filePath = 'D:\HuaweiMoveData\Users\Darian\Desktop\MathModeling\附件\附件1.xlsx';

if ~isfile(filePath)
    error('未找到文件：%s', filePath);
end

data = readtable(filePath);
time     = data{:,1};   % 时间
temp     = data{:,2};   % 温度
moisture = data{:,3};   % 水分浓度

fprintf('成功读取文件，共 %d 行数据\n', height(data));


win = 15;

temp_smooth     = smoothdata(temp,     'gaussian', win);
moisture_smooth = smoothdata(moisture, 'gaussian', win);


outputFile = fullfile(pwd, '附件1_去噪后.xlsx');

outTable = table(time, temp_smooth, moisture_smooth, ...
    'VariableNames', {'时间', '温度', '水分浓度'});

writetable(outTable, outputFile);
fprintf('去噪数据已保存到：%s\n', outputFile);

fig = figure('Color', 'w', 'Position', [100 100 1000 750]);

subplot(2,1,1);
plot(time, temp, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2); hold on;
plot(time, temp_smooth, '-', 'Color', [0.85 0.20 0.20], 'LineWidth', 2);
xlabel('时间', 'FontSize', 12);
ylabel('温度', 'FontSize', 12);
title('温度：去噪前后对比', 'FontSize', 13, 'FontWeight', 'bold');
legend('原始数据', '去噪后数据', 'Location', 'best', 'FontSize', 10);
grid on; box on;

subplot(2,1,2);
plot(time, moisture, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2); hold on;
plot(time, moisture_smooth, '-', 'Color', [0.10 0.45 0.85], 'LineWidth', 2);
xlabel('时间', 'FontSize', 12);
ylabel('水分浓度', 'FontSize', 12);
title('水分浓度：去噪前后对比', 'FontSize', 13, 'FontWeight', 'bold');
legend('原始数据', '去噪后数据', 'Location', 'best', 'FontSize', 10);
grid on; box on;

figFile = fullfile(pwd, '去噪前后对比图.png');
exportgraphics(fig, figFile, 'Resolution', 200);
fprintf('对比图已保存到：%s\n', figFile);

rmse_temp = sqrt(mean((temp - temp_smooth).^2));
rmse_moist = sqrt(mean((moisture - moisture_smooth).^2));
fprintf('\n=== 去噪残差（RMS）评估 ===\n');
fprintf('温度      RMS = %.6f\n', rmse_temp);
fprintf('水分浓度  RMS = %.6f\n', rmse_moist);