% 51.166944
function retval = Process(s)
    [T,C,t] = Process4(0.001/9,2,s.hT,s.hm,s.D0);
    retval = t/3600;
end