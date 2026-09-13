function geometry = ShrinkingRadius(timeSeconds,N)

    validateattributes(N, {'numeric'}, ...
        {'scalar','real','finite','integer','>=',2});

    R = RadiusMeters(timeSeconds);
    dr = R/N;
    r = (0:N)'*dr;

    V = r*dr;
    V(1) = dr^2/8;
    V(end) = R*dr/2 - dr^2/8;

    geometry.R = R;
    geometry.dr = dr;
    geometry.r = r;
    geometry.V = V;
end