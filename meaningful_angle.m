function out = meaningful_angle(in)

out = zeros(size(in));

for i = 1:length(in)
    
    if rem(in(i),2*pi) < 0
        out(i) = rem(in(i), 2*pi) + 2*pi;
    else
        out(i) = rem(in(i), 2*pi);
    end
end
