function [p, histo, ptheory]=wantpdf(sig2esti,Ncells,m,variance)

% Function WANTPDF Calculates the PDF of the signal sig2esti
% N cells Specifies the number of cells for Histogram estimation
% and if theoretical PDF (Gaussian) is needed, input the mean and variance 

% Finding the PDF
M=length(sig2esti);
a=max(sig2esti);
b=min(sig2esti);
range=a-b;
Hepsi=(range/(2*Ncells));
values=(b+Hepsi):(2*Hepsi):a;

for j=1:Ncells,
    temp=sig2esti-values(j);
    counter =abs(temp) <= Hepsi;
    histo(j)=sum(counter)/M; % Histogram Estimator
    p(j)=histo(j)/(2*Hepsi); % PDF estimator
end

% Theoretical  Asymptotic PDF of the estimator for A=1

ptheory=inv(sqrt(2*pi*variance))*exp(-(values-m).^2/(2*variance)); 