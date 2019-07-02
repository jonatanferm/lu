function [w1, w2] = payoff(s1, s2, n, U)
    i = min([find(s1 == 2), find(s1 == 2)]);
    if numel(i) == 0
        i = n;
    end
    w1 = sum(U(sub2ind(size(U), s1(1:i), s2(1:i))));
    w2 = sum(U(sub2ind(size(U), s2(1:i), s1(1:i))));
end
