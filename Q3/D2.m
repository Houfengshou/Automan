function retval = D2(T, C)
    retval = 2.4e-3*exp(-0.45/C)*exp(-3850/(T + 273.15));
end