function v = Vector1N(x)
% 强制转行向量
    [s1,s2] = size(x);
    if s1 == 1
        v = x;
    elseif s2 == 1
        v = x.';
    else
        error('Vector1N:输入不是一维向量');
    end
end
