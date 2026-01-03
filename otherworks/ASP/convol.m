function y=convol(x,h)

X=convolmtx(x,length(h));
y=X*h'; % convoled result out as a column vector;