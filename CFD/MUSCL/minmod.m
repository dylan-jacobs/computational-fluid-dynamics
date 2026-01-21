function [output] = minmod(a, b)
    output = 0.5*(sign(a)+sign(b)).*min(abs(a), abs(b));
end