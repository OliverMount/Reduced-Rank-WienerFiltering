function num=count(x,val)

k=length(x);
num=0;
for i=1:k,
    if x(i)==val,
    num=num+1;
    end
end