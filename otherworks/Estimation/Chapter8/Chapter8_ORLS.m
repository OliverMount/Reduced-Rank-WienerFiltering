% This program for Order Recursive Least squares (ORLS)

% source Generation
clear all;clc;close all;
n=0:99;
N=length(n);
theta_org=[1 .03]'; % Original parameters

H_m=[ones(N,1) n']; % Original source H for data generation
H=[ones(N,1) n' n'.^2 n'.^3]; % H for recursive estimation

% Model for Source (data to be modeled)
x=H_m*theta_org + (.3)*randn(1,N)'; % variance of the Noise .0625

p=3; % No of parameters to estimate

% Initialization of Least Squares

D=inv(H(:,1)'*H(:,1));
theta=D*H(:,1)'*x;

theta_comp{1,1}=theta;
Jmin(1)=x'*(x-H(:,1)*theta(:,1));
P=H(:,1)*D*H(:,1)';
I=eye(N,N);

for k=1:p-1,
    
    a=H(:,k+1)'*[I-P]*H(:,k+1);
    alpha=H(:,k+1)'*[I-P]*x/(a);
    b=H(:,(1:k))'*H(:,k+1);
    c=-D*b/a;
    d=-b'*D/a;
    
    theta=[theta-D*b*alpha ; alpha]; % Estimator update
    theta_comp{1,k+1}=theta;   % storing of Intermediate estimators starting from lower order
    D=[D+(D*b*b'*D)/a  c ; c' 1/a];  % Recursive update for D
    Jmin(k+1)= Jmin(k)-a*(alpha^2); % Jmin updation
    
    P=P+([I-P]*H(:,k+1)*H(:,k+1)'*[I-P]/a); % updation of the projection matrix
end

figure,scatter(n,x);hold on;
hnd=title('Experimental data (Original Two parameter)');
set(hnd,'fontsize',15); 
hnd=xlabel('Time (t)');
set(hnd,'fontsize',15); 

for i=1:p,
mtx=[1 0 0;0 0 1;0 0 0];
hnd=plot(n,H(:,(1:i))*theta_comp{1,i}); % If u vary the parameter change here..
set(hnd,'color',mtx(i,:));
end

legend('One Parameter Fit','Two Parameter Fit','Three Parameter Fit','Original data');

figure,plot(1:length(Jmin),Jmin);