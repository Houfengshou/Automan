function out = AntiRadixIterator(dims)
%ANTIRADIXITERATOR 迭代器模式，末位优先变化
% 调用方式
%   AntiRadixIterator([2,3,4])   % 传入维度，初始化迭代器，无有效输出
%   out = AntiRadixIterator      % 无输入参数，获取下一个组合；返回[]代表迭代结束
    persistent current dims_val n is_done
    
    if nargin > 0
        if ~isvector(dims) || size(dims,1)~=1
            error('输入必须为行向量');
        end
        dims_val = dims(:)';
        n = length(dims_val);
        current = ones(1,n);
        is_done = false;
        out = [];
        return;
    end
    
    if isempty(is_done)
        error("迭代器未初始化，先调用 AntiRadixIterator(dims)");
    end
    
    if is_done
        out = [];
        return;
    end
    
    out = current;
    
    carry = 1;
    for i = n:-1:1
        if carry == 0
            break;
        end
        current(i) = current(i) + 1;
        if current(i) > dims_val(i)
            current(i) = 1;
            carry = 1;
        else
            carry = 0;
        end
    end
    
    if all(current == 1)
        is_done = true;
    end

end
