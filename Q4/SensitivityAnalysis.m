function result = SensitivityAnalysis(func, ParamStruct, ploting, saving, ax)
    if nargin == 2
        ploting = 0;
        saving = 0;
        ax = [];
    elseif nargin == 3
        saving = 0;
        ax = [];
    elseif nargin ==4
        ax = [];
    end

    fields = getStructFields2Vec1N(ParamStruct);

    constnum = 0;
    definiens = []; % {字段, ... } [字段对应的值; ... ]
    definiendoms = string([]);

    varnum = 0;
    variables = []; % {字段, ... } [起点, 终点, 个数; ... ]
    varnames = string([]);

    for field = fields
        value = ParamStruct.(field);
        if all(size(value) == [1,1])
            constnum = constnum + 1;
            definiens(constnum,1) = value;
            definiendoms(constnum,1) = field;
        elseif all(size(value) == [1,3])
            varnum = varnum + 1;
            variables(varnum,1:3) = [value(1), value(2), value(3)];
            varnames(varnum,1) = field; 
        else 
            error("错误输入格式");
        end
    end
    if varnum == 0
        error("变量个数为0");
    elseif varnum > 2
        warning("变量超过两个, 绘图功能失效");
    end
%--------------------------------------------------------------------------
    ArgumStruct = struct();
    for constidx = 1:constnum
        ArgumStruct.(definiendoms(constidx,1)) = definiens(constidx,1);
    end    
    AntiRadixIterator(Vector1N(variables(:,3)));
    result = zeros(prod(Vector1N(variables(:,3))), varnum + 1); % 实验次数*(变量数+结果)
    domin = cell(1, varnum);

    for varidx = 1:varnum
        domin{varidx} = exp(linspace(variables(varidx, 1),  variables(varidx, 2), variables(varidx, 3)));
    end
    index = AntiRadixIterator();
    i = 1;
    while ~isempty(index)
        for varidx = 1:varnum
            domin1N = domin{varidx};
            ArgumStruct.(varnames(varidx,1)) = domin1N(1,index(1,varidx));
            result(i,varidx) = domin1N(1,index(1,varidx));
        end
        result(i,varnum+1) = func(ArgumStruct);
        index = AntiRadixIterator();
        i = i + 1;
    end

    if ploting
        if varnum == 1
            if isempty(ax)
                figure;
                ax = gca;
            end
            hold(ax,'on');    % 指定在ax坐标轴开启hold
            plot(ax, result(:,1), result(:,2));
            xlabel(ax, varnames(1,1));
            ylabel(ax, "result value");
        end
        if varnum == 2
            if isempty(ax)
                figure;
                ax = gca;
            end
            hold(ax,'on');
            [X,Y] = meshgrid(domin{1}, domin{2});
            Z = reshape(result(:,3), variables(2,3), variables(1,3));
            contourf(ax,X,Y,Z,10);
            shading(ax,'interp');
            cb = colorbar(ax);
            xlabel(ax, varnames(1,1));
            ylabel(ax, varnames(2,1));
            cb.Label.String = '目标函数值';
        end
    end


    if saving
        colNames = [cellstr(varnames); {'result'}];
        T = array2table(result, 'VariableNames', colNames);
        writetable(T, 'SensitivityResult.csv');
        fprintf('SensitivityResult.csv\n');
    end

end
