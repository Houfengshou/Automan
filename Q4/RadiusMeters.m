function R = RadiusMeters(timeSeconds)

    validateattributes(timeSeconds, {'numeric'}, ...
        {'scalar','real','finite','nonnegative'});

    timeHours = timeSeconds/3600;

    radiusCm = 1.1998 + 0.8002*( ...
        0.2082*exp(-timeHours/1.1588) ...
        + 0.7918*exp(-timeHours/4.6596));

    R = radiusCm/100;
end