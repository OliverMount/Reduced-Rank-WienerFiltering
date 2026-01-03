% ***********SIMULATION OF COEFFICIENTS OF RAYLEIGH FADING CHANNEL************
%
% VEHICULAR SPEED, DOPPLER SHIFT
% FADING COEFFICIENTS
%******************************************************************************

Num_path=2000;							% Number of paths
t=0.0001:10/Num_path:10;				% Time range
f=150*10.^6; 							% Carrier frequency (150 Mhz, 900 Mhz)
wc=2*pi*f;								
vehicle_speed=150;                      % Speed of car[km/hrs]
v=vehicle_speed*5/18;					% Receiver speed[m/hrs]
c=300*10^3;								% Speed of light
wm=wc*(v/c);							% Maximum shift
fm=wm/(2*pi);							% Doppler shift

% SIMULATING ENSEMBLES OF SINUSOIDS
for i=1:Num_path
    A(i)=(2*pi/Num_path)*i;				%Azimuthal angles
    wn(i)=wm*cos(A(i));
    phi(i)=(pi*i)/(Num_path+1);
    xc(i)=2*cos(wn(i)*t(i)).*cos(phi(i))+cos(wm*t(i));
    xs(i)=2*cos(wn(i)*t(i)).*sin(phi(i));
    T(i)=(1/(2*Num_path+1)^0.5).*(xc(i)+j*xs(i));% Complex envelope
end

M=mean(abs(T));							% Mean
MdB=20*log10(M);
TdB=floor(20*log10(abs(T)));		    % Field [dB]

% PLOTTING THE HISTOGRAM
z1=hist(abs(T));
z=hist(TdB,9);
n=0;
for k=1:9
   n=n+z(k);
   end
   	for j=1:9
    	P(j)=z(j)/n;
 	end
  f(1)=P(1);
	for x=2:9
		f(x)=f(x-1)+P(x);
  		F(10-x)=f(x);
    end
plot(z1); 								 % Distribution chart
title('Rayleigh distribution');
semilogy(t,abs(T)/max(abs(T)),'r')       % Fading graphic
title('Received field');
ylabel('Received field intensity');
xlabel('time');
grid on;