function y=corr(x,h)

x=flipud(x.');
X=convolmtx(x.',length(h));
y=X*h.'; % Correlated result out as a column vector;