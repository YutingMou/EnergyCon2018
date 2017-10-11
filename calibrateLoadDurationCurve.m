%%% This script is used to calibrate load duration curve using a nonlinear
%%% function, in order to get valuation function.
%%% inputs:
%%%	inputLDC - the LDC you want to calibrate
%%% p0 - average electricity price
%%% elasticity - average elasticty during the whole horzion
%%% outputs: parameters a0, a1, b, labmda
%%% According to the theory in EnergyCon2018Appendix.pdf
%%% author: Yuting Mou
%%% date: 2017-10-10
%%% email: yuting.mou(AT)outlook.com

%% step 0: specify your input data
inputLDC = 720:-1:1;    % inputLDC is a row vector
p0 = 213; 
elasticity = -0.13;   % adapt to your data


%% step 1: find the best lambda using least square
H = max(size(inputLDC)); % the horizon of the LDC
t = 1:H;
lambdaVec = 0.1:0.1:2;
lambdaN = size(lambdaVec,2);
error = zeros(lambdaN,1);
% traverse over all the instances of lambda, do polyfit and calculate error
for i = 1:lambdaN
    lambda = lambdaVec(i);
    temp1 = (t./H).^(1/lambda);
    temp2 = inputLDC;
    temp = polyfit(temp1,temp2,1);   
    
    tempk = temp(1);
    tempb = temp(2);
    
    tempx= linspace(1,H, H); % Adapt n for resolution of graph
    tempy = tempk*(tempx./H).^(1/lambda)+tempb;
    
    error(i) = sum((inputLDC - tempy).^2);  
end
% plot lambda-error
figure;
plot(lambdaVec,error);
xlabel('\lambda');
ylabel('Error');

% find the labmda that error is minimized
[~,lambdaIdx] = min(error); % the index of the best lambda
lambda = lambdaVec(lambdaIdx);

temp1 = (t./H).^(1/lambda);
temp2 = inputLDC;
temp = polyfit(temp1,temp2,1);

tempk = temp(1);
tempb = temp(2);

%% step 2: calculate othere parameters
% tempk = -a1/b; tempb = -p0/b + a0/b;
% => a1 = -tempk*b; a0 = tempb*b+p0;

syms x b % 
solutionb = vpasolve( int(  1/( b*tempb + b*tempk*x^(1/lambda) ) , x,0,1) == elasticity/(-p0),b);
% if the calculation above is too slow, you can reduce the resolution of
% lambda. For example, if lambad = 1.99, you can  use 2 instead.
a1 =  -tempk*solutionb;
a0 = p0 + solutionb*tempb;

%% step 3: compare the calibrated LDC with the actual one
tempx= linspace(1,H, H); % Adapt n for resolution of graph
tempy = tempk*(tempx./H).^(1/lambda) + tempb;
tempy2 = -a1/solutionb*(tempx./H).^(1/lambda) - p0/solutionb + a0/solutionb;
figure;
plot(inputLDC,'linewidth',1);
hold on;
plot(tempy2,'linewidth',1);
legend('Actual','Calibrated')
xlabel('Duration (Hours)');
ylabel('Power (MW)')
