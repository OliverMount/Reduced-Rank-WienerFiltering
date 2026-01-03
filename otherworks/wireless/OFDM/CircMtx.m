function H=CircMtx(h,N)

% h-impulse response vector
% N-No of subcarriers 
L=length(h)-1; %L-Channel Memory
tem=[h zeros(1,N-(L+1))];
for i=1:N,
    H(i,:)=wshift(1,tem,-(i-1));
end