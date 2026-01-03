function x=ForSubs(a,b)

% Algorithm for Forward substitution for solving Ax=b;
[m,n]=size(a);

for i=1:m,
    temp=0;
    for j=1:i-1,
        temp=temp+a(i,j)*x(j);
    end
    x(i)=(b(i)-temp)/a(i,i);
end

x=x';

% Test matrices for verification

% A =[  3     0     0     0
%     -1     1     0     0
%      3    -2    -1     0
%      1    -2     6     2];
%  
%  b= [5 6 4 2];
%  
% x=[ 1.6667  7.6667 -14.3333   50.8333]; % answer