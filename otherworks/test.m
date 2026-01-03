w=-pi:pi/20:pi;
C=cos(w);S=sin(w);

a=-.5;
R=(1-(a.*C))./(1-(2*a.*C)+(a^2));
I=(-a*S)./(1-(2*a.*C)+(a^2));
M=sqrt(1./(1-(2*a.*C)+(a^2)));
A=atan2(-a*S,(1-(a.*C)));
plot(w,I);