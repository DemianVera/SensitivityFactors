function sal = dK(m,n,q)
if nargin >= 3
if m == n || m == q
    sal = 1;
else
    sal = 0;
end



end

if nargin <3 
if m == n
    sal = 1;
else
    sal = 0;
end

end
end