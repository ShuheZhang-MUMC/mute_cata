function F = get_F(n)
F = ones(1,n+2);
for n = 3:n+2
    F(n) = F(n-1) + F(n-2);
end
end