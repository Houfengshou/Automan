function [Tnew,Cnew,converged,iter] = Picard( ...
    Told,Cold,r,V,dr,hT,hm,tau,theta,qTold,qCold,qTnext,qCnext,options)

    Tnew = Told;
    Cnew = Cold;
    converged = false;
    iter = 0;

    [MT,MC,KT,KC] = AssembleDryingSpatialMatrices( ...
        Told,Cold,r,V,dr,hT,hm);

    FTold = MT \ (qTold-KT*Told);
    FCold = MC \ (qCold-KC*Cold);

    baseT = Told + (1-theta)*tau*FTold;
    baseC = Cold + (1-theta)*tau*FCold;

    Tg = Told;
    Cg = Cold;

    for iter = 1:options.maxIterations

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

        [MT,MC,KT,KC] = AssembleDryingSpatialMatrices( ...
            Tc,Cc,r,V,dr,hT,hm);

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