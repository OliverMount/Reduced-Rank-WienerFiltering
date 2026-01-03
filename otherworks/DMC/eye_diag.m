

beta = 0.35;
N = 100;
t = -10.001:0.1:10;
h = sinc(t).*(cos(pi*beta*t)./(1 - (2*beta*t).^2));

figure;

for ind=1:1000
  d = 2*round(rand(1,N))-1;
  d1 = reshape([d; zeros(9,N);],1,10*N);
  y = conv(d1,h);
  plot(0:0.1:0.1*(length(y)-1),y);
  hold on;
end
axis([40 42 -3 3]);
grid on;
