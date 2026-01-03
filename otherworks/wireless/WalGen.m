% This program generates the Walsh codes for a given length N
function W=WalGen(N)
h2=[1 1;1 -1];
W=h2;
m=log10(N)/log10(2);

for i=1:m-1,
    tempvar=[W W;W -W];
    clear W;
    W=tempvar;
    clear tempvar;
end