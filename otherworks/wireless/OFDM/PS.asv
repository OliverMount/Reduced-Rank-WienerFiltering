function  Y=PS(y,N)

% x- data sequence
% N - no of OFDM symbols

[m,n]=size(y);
for i=1:n,
    Y(1+(i-1)*N:i*N)=fliplr(y(:,i).');
end
