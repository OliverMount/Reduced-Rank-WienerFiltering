function  X=SP(x,N)

% x- data sequence
% N - no of OFDM symbols

bl=length(x)/N; % Block length
for i=1:bl,
    X(:,i)=flipud(x((1+(i-1)*N):(i*N)).');
end
