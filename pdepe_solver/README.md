# 温度与含水率PDEPE求解器

从烘干开始求解任意指定时刻的T、C。采用已确认的附录3模型、固定半径2厘米；4小时前使用仓库去噪数据，4小时后空气50°C、水分浓度0.05 kg/kg。在4小时处重新启动求解，初值接上一段的完整场。

## 文件及入口

- `SolveDryingTemperatureMoisturePDEPE.m`：主求解器。内部包含明确命名的物理方程、环境、初值、边界和事件函数。
- `EvaluateDryingTemperatureMoisture.m`：查询已经求得的解。
- `ExampleQueryDryingTemperatureMoisture.m`：可直接运行的查询示例。
- `ExportDryingMoistureResultTables.m`：导出第三问结果草稿。
- `source/附件1去噪后.xlsx`：已附带数据，运行时不依赖其他Q1/Q2文件。

在MATLAB中进入本目录，或将本目录加入路径。

## 像Process2一样使用

```matlab
timeSeconds = [100, 300, 600, 1800, 3600, 10800];
radiusCm = [0, 0.5, 1, 1.5, 2];
[T, C, solution] = SolveDryingTemperatureMoisturePDEPE(timeSeconds, radiusCm);
```

T和C均为6行5列：行对应timeSeconds，列对应radiusCm。T单位°C，C单位kg/kg，返回完整浮点精度。不要在递推或阈值判断前舍入。

时间单位始终是秒：查询6小时写`6*3600`。半径单位是厘米：表面写`2`。查询数组须严格递增。允许小数时间及半径，例如：

```matlab
[T, C] = SolveDryingTemperatureMoisturePDEPE(1234.5, [0, 0.7, 2]);
```

指定时间没有0也会从0求解，最大指定时间为求解终点。主函数指定的查询时刻直接加入pdepe输出时刻；半径由pdeval作空间求值。输入任意半径不会改变计算网格。

## 求解后再查询

```matlab
[Tnew, Cnew] = EvaluateDryingTemperatureMoisture(solution, [1234.5, 5678.9], [0, 0.7, 2]);
```

保存时刻直接取值；其他时刻采用PCHIP时间插值，再以pdeval求空间值。这属于对保存结果的插值，不是ODE求解器的连续输出。对精度敏感的查询，应在首次求解时就指定，或减小StorageIntervalSeconds。禁止在已求解时段或0–2厘米之外外推。

## 第三问阈值搜索

```matlab
settings = struct('StopWhenDry', true, 'MeshRefinement', 1);
[T, C, solution] = SolveDryingTemperatureMoisturePDEPE( ...
    0:60:144*3600, 0:0.1:2, settings);
solution.dryingThresholdTimeHours
ExportDryingMoistureResultTables(solution, fullfile(pwd, 'results'));
```

144小时为示例搜索上限，不是预设烘干时间。若`solution.thresholdReached`为false，时间字段为NaN，需扩大范围后再求解，不能把上限当成答案。

当StopWhenDry为true时，返回的T、C只保留终止之前可查询的行；行对应`solution.returnedTimeSeconds`。原请求是否有效由`solution.requestedTimeAvailable`标记。阈值瞬间的场见`thresholdTemperatureCelsius`和`thresholdMoistureDryBasis`，以及完整时间数组的末行。

事件定位的是`max(C)=0.15`的下降穿越点，不是严格小于0.15的瞬间。严格达标可在关闭停止选项后，对阈值稍后的时刻求解并核对，时间增量和容差应与所需精度相符。不要把四位小数显示的0.1500当作事件判断值。

导出文件为`result3_pdepe_draft.xlsx`与`table5_pdepe_draft.xlsx`。表格补充阈值时刻，按四位小数舍入；MAT文件保留完整精度。文件名带draft表示尚未作为最终提交认定。

## 数值设置

|设置字段|默认值|含义|
|---|---|---|
|MeshRefinement|1|平方映射网格的加密倍数，正整数；比较1、2、4用于检查空间误差|
|SpatialMeshMeters|空|可直接提供严格递增、起点0终点0.02的空间网格（米）；提供后覆盖自动网格|
|RelativeTolerance|1e-6|时间积分相对容差|
|AbsoluteTolerance|[1e-7;1e-9]|温度、含水率的绝对容差|
|MaximumInternalStepSeconds|300|内部时间步上限；非固定步长|
|StorageIntervalSeconds|60|保存完整空间场的时间间隔；不是内部步长|
|StopWhenDry|false|是否在所有节点到达含水率阈值时停止|
|EnvironmentDataFile|随包提供的数据|可更改数据文件路径|

自动网格为r=0.02[1-(1-s)^2]，s等距划分120×MeshRefinement个区间，叠加题目0.1厘米间隔输出节点。默认值用于启动数值对照，不是已经认证精度的最终网格。最终应以网格、容差加密，以及独立有限体积结果比较确定精度。

## 与Process2的区别

不再手写时间循环和圆心/内部/表面递推，不调用HarMean。PDEPE自行离散并以刚性ODE求解器推进。温度、含水率控制方程和物性关系仍与已确认模型一致，未加入收缩或潜热。现有127.9392小时与65.0198小时来自另一个未空间收敛的显式程序，不能作为本程序的正确答案基准。

官方接口依据：[pdepe](https://www.mathworks.com/help/matlab/ref/pdepe.html)、[pdeval](https://www.mathworks.com/help/matlab/ref/pdeval.html)。

## 本机验证结果

已实际运行MATLAB验证：初值、单时刻查询、小数时间、小数半径、4小时衔接、越界拒绝、阈值停止及Excel输出均通过。

|MeshRefinement|节点数|阈值时间(h)|
|---|---:|---:|
|1|139|57.4737546862|
|2|259|57.4726415273|
|4|499|57.4723244010|

在第4级设置上进一步收紧时间容差及步长上限，得到57.4723987891小时，变化约0.2678秒。对该阈值前后各1秒重新求解，最大含水率分别为0.1500002913与0.1499997087，确认下降穿越。

这些检查支持程序可用于下一步计算。当前空间加密末两级仍差约1.14秒，不应宣称烘干时间已精确到小时的四位小数；同一模型的独立有限体积交叉核对尚未完成。`spatialConvergenceVerified`因此保持false，当前约57.47小时仍是待最终确认的数值估计。

`ValidateDryingTemperatureMoistureSolver`可复现检查，详细指标见validation_results/validation_report.json。打包代码中的该报告记录本次检查；完整MAT和Excel另保存在工作区validation_results目录。
