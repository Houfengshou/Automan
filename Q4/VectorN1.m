function v = VectorN1(x)
% 强制转列向量
    [s1,s2] = size(x);
    if s2 == 1
        v = x;
    elseif s1 == 1
        v = x.';
    else
        error('VectorN1:输入不是一维向量');
    end
end
