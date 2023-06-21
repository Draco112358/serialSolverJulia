function    [out] = round_ud(a,b)

if nargin<2
    out=round(a);
else
    out=round(a*10^b)/10^b;
end

end