function A=LSDataform(x,M,opt)

% LSDataform(x,M) Returns the Data matrix for LS filtering
% x- Input data sequence -  in row wise input.
% M- length of the LS filter
% Option specifies the type of forming  data matrix
% 1-Auto correlation method 2- Auto Covariance method 3- Pre windowing method 4 - Post
% windowig method

for i=1:M,
    B(:,i)=[zeros(1,i-1) x  zeros(1,M-i)]';  % column wise implementation 
end 

a=length(x);

switch opt
    case 1
    disp(' Auto Correlation Method');
    A = B;
    case 2
    disp(' Auto Covaiance Method');
    A = B(M:a,:);  
    case 3
    disp(' PreWindowing Method ');
    A = B(1:a,:);       
    case 4
    disp(' Post Windowing Method ');
    A = B(M:end,:);  
end

        
    