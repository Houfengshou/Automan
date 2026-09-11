function [temperatureCelsius,moistureDryBasis]= ...
    EvaluateDryingTemperatureMoisture(solution,queryTimeSeconds,queryRadiusCm)
% 查询已保存解。返回矩阵行为时间、列为半径。
% 已保存时刻直接取值；新增时刻使用PCHIP时间插值。
% 空间位置使用pdeval。禁止时间或空间外推。
validateattributes(queryTimeSeconds,{'numeric'},{'real','finite','nonnegative'});
validateattributes(queryRadiusCm,{'numeric'},{'vector','real','finite','nonempty','>=',0,'<=',2});
queryTimeSeconds=double(queryTimeSeconds(:));
queryRadiusMeters=double(queryRadiusCm(:)')/100;
temperatureCelsius=zeros(numel(queryTimeSeconds),numel(queryRadiusMeters));
moistureDryBasis=temperatureCelsius;
if isempty(queryTimeSeconds), return; end
assert(all(queryTimeSeconds>=solution.timeSeconds(1) & queryTimeSeconds<=solution.timeSeconds(end)), ...
    '查询时间超出已求解范围。请重新求解到所需时间，或关闭StopWhenDry。');
for row=1:numel(queryTimeSeconds)
    savedRow=find(solution.timeSeconds==queryTimeSeconds(row),1);
    if ~isempty(savedRow)
        temperatureProfile=solution.temperatureCelsius(savedRow,:);
        moistureProfile=solution.moistureDryBasis(savedRow,:);
    else
        temperatureProfile=interp1(solution.timeSeconds,solution.temperatureCelsius,queryTimeSeconds(row),'pchip');
        moistureProfile=interp1(solution.timeSeconds,solution.moistureDryBasis,queryTimeSeconds(row),'pchip');
    end
    temperatureCelsius(row,:)=pdeval(1,solution.radiusMeters,temperatureProfile,queryRadiusMeters);
    moistureDryBasis(row,:)=pdeval(1,solution.radiusMeters,moistureProfile,queryRadiusMeters);
end
end
