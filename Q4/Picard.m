function [Tnew,Cnew,converged,iter] = Picard( ...
    Told,Cold,geometryOld,geometryNext,hT,hm, ...
    tau,theta,qTold,qCold,qTnext,qCnext,options)

    % 输入均使用列向量
    Tnew = Told;
    Cnew = Cold;
    converged = false;
    iter = 0;

    [MT,MC,KT,KC] = AssembleDryingSpatialMatrices( ...
        Told,Cold,geometryOld.r,geometryOld.V,geometryOld.dr,hT,hm);

    FTold = MT \ (qTold-KT*Told);
    FCold = MC \ (qCold-KC*Cold);

    baseT = Told + (1-theta)*tau*FTold;
    baseC = Cold + (1-theta)*tau*FCold;

    Tg = Told;
    Cg = Cold;
% 以旧状态作为物性初猜，但使用新时间层的几何
    [MT,MC,KT,KC] = AssembleDryingSpatialMatrices( ...
        Tg,Cg,geometryNext.r,geometryNext.V,geometryNext.dr,hT,hm);
    for iter = 1:options.maxIterations

        % 固定当前物性，求解候选状态
        Tc = (MT+theta*tau*KT) \ (MT*baseT+theta*tau*qTnext);
        Cc = (MC+theta*tau*KC) \ (MC*baseC+theta*tau*qCnext);

        if any(~isfinite([Tc;Cc])) || ...
           any(Tc <= -273.15) || any(Cc <= 0)
            return
        end

        scaleT = options.temperatureAbsoluteTolerance ...
            + options.relativeTolerance.*max(abs(Tc),abs(Told));
        scaleC = options.moistureAbsoluteTolerance ...
            + options.relativeTolerance.*max(abs(Cc),abs(Cold));

        % 用候选状态更新物性，检查非线性残差
        [MT,MC,KT,KC] = AssembleDryingSpatialMatrices( ...
            Tc,Cc,geometryNext.r,geometryNext.V,geometryNext.dr,hT,hm);

        FT = MT \ (qTnext-KT*Tc);
        FC = MC \ (qCnext-KC*Cc);

        resT = Tc-Told-tau*((1-theta)*FTold+theta*FT);
        resC = Cc-Cold-tau*((1-theta)*FCold+theta*FC);

        if any(~isfinite([resT;resC]))
            return
        end

        err = max([abs(Tc-Tg)./scaleT;
                   abs(Cc-Cg)./scaleC;
                   abs(resT)./scaleT;
                   abs(resC)./scaleC]);

        if err <= 1
            Tnew = Tc;
            Cnew = Cc;
            converged = true;
            return
        end

        Tg = Tc;
        Cg = Cc;
    end
end