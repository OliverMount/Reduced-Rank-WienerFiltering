function r=convolmtx(x,M)

% CONVOLMTX(X,M) Returns the Length(x)+M-1 X M convolution matrix 
% x- Input data sequence -  in row wise input.
% M- length of the filter

for i=1:M,
    r(:,i)=[zeros(1,i-1) x  zeros(1,M-i)].';  % column wise implementation 
end 