function [cor] = CORRELATION_func(X1,X2)
    if(sum((X1-X2).^2)>1)
        cor = 0;
    else
        cor = 1-(sum((X1-X2).^2));
    end
end