%%% This script is used to disign a price menu accrodinig to my submission
%%% to EnergyCon2018
%%% inputs: 
%%% value function (parameters from load duration curve)
%%% cost function
%%% h(\omega)
plotFigures = 1; % swith to enable ploting figures
%% step 0: specify your input data
h_omega = ;
omegaIdxN = max(size(h_omega)); % the number of omega
Lmax = ; % max L index, when you do ecomonic dispatch
L0Calibrated = ;    % calibrated LDC when duration becomes 0
Hr_vec = zeros(1,omegaIdxN);
for i = 1:omegaIdxN
    Hr_vec(i) = sum(h_omega(1:i))/omegaIdxN;
end
cOmega=; % cost at load level L in senario omega

%% step 1: c(L,r) and R(L)
cOmega(MC==3000)=3000; % get the index when MC==3000 and set c to 3000 meaning this slice is not served
c = cOmega;
% check R(L) of each L level
RL = zeros(Lmax,1);
for L = 1:Lmax
    
    tempR = find(c(L,:)<3000,1,'last');
    if ~isempty(tempR)
    RL(L) = tempR;
    end
    for j =RL(L):-1:1
        c(L,j) = sum(c(L,1:j),2)/omegaIdxN;
    end
    
end
%% step 2: r(L) and t(L)
% solve this problem. traverse over all the cases of reliablity_vec and find the best objective value
% parameters of value function
%VLt = -a1*(1/T)^(1/lambda)/(1/lambda+1)*t^(1/lambda+1)+(a0-b*L)*t
% simplied as VLt = m*t^(1/lambda+1)+(a0-b*L)*t;
lambda =  ;  %1.449;
a0 =  ;
a1 =  ;
b =    ;
p0 = ;
T = ;
m = -a1*(1/T)^(1/lambda)/(1/lambda+1);

reliablity_vec = 1/omegaIdxN:1/omegaIdxN:1;

% optimal solution
duration = zeros(Lmax,1);
reliablity = zeros(Lmax,1);
welfare = zeros(Lmax,1);
welfareIdx = zeros(Lmax,1); % the index of best welfare in each level, i.e.,  the best reliability

duration_temp = zeros(Lmax,omegaIdxN);
reliablity_temp = zeros(Lmax,omegaIdxN);
welfare_temp = zeros(Lmax,omegaIdxN);

for L = 1:Lmax
    for i = 1:RL(L)
        reliablity_temp(L,i) = reliablity_vec(i);
        t_root_temp = (c(L,i) + Hr_vec(i)*(b*L-a0))/(Hr_vec(i)*-a1*(1/T)^(1/lambda));
        if t_root_temp>T^(1/lambda)
            duration_temp(L,i) = T;
        else
            if t_root_temp<0
                duration_temp(L,i) = 0;
            else
                duration_temp(L,i) = t_root_temp^lambda;
            end
        end      
        welfare_temp(L,i) = Hr_vec(i)*(m*duration_temp(L,i)^(1/lambda+1)+(a0-b*L)*duration_temp(L,i)) - c(L,i)*duration_temp(L,i);        
    end   
end

% process the results, find the best welfare
for L = 1:Lmax
    [welfare(L),welfareIdx(L)] = max(welfare_temp(L,:));
    
    if L>1 && welfareIdx(L)>welfareIdx(L-1)
        welfareIdx(L) = welfareIdx(L-1);
    end
    duration(L) = duration_temp(L,welfareIdx(L));
    reliablity(L) = reliablity_temp(L,welfareIdx(L));
end

% find L0
L0 = find(welfare>1e-6, 1,'last'); % find the last slice level with a welfare that is bigger than 0;
L_supplied = 1:L0;                 % the index of slices that are served

% plot optimal duration
% note, you may need an 'ironing' procedure, or just sort duration
% descendingly
duration = sort(duration,'descend');
if plotFigures
    figure;
    plot(L_supplied,duration(L_supplied),'linewidth',1);
    xlabel('L (MW)');
    ylabel('Duration (hours)');
end

%% step 3: design the menu and obtain p(L), f(L) and g(L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate PL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hr = zeros(L0,1);
PL = zeros(L0,1);
HrVL = zeros(L0,1); % the data to be integrated, in numeric integration
for L = 1:L0
    
    r_Idx_L = welfareIdx(L);
    PL(L) = m*duration(L)^(1+1/lambda) + (a0-b*L)*duration(L); % VLt = m*t^(1/lambda+1)+(a0-b*L)*t, the first term in the equation
    Hr(L) = Hr_vec(r_Idx_L);    
end

% using Integration of Numeric Data, cumtrapz function
% prepare the data
for L = 1:L0
    HrVL(L) = Hr(L)*duration(L)*(-b);
end
% calculate PL
temp1 = cumtrapz(HrVL);
for L = 1:L0
    PL(L) = PL(L) + 1/Hr(L)* (temp1(L0) - temp1(L));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate f(t) and  calculate g(r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v_t = -a1*(t/T)^(1/lambda) -b*L + a0
Lt = 1:L0;
Lt(1:find(duration==T,1,'last')) = 1; % find the last L whose duation is 8760
vL0t0 = m*duration(L0)^(1+1/lambda) + (a0-b*L0)*duration(L0);
fL = vL0t0*ones(L0,1);
ft = zeros(T,1);
% using Integration of Numeric Data, cumtrapz function
% prepare the data to be integrated
Vtlt_seeds = zeros(L0,1);   % not all the instances
for L = 1:L0
    Vtlt_seeds(L) = -a1*(duration(L)/T)^(1/lambda)  - b*Lt(L) + a0;
end

% discretized implementation
for L = 2:L0
    for j = L+1:L0
        fL(L) = fL(L) - Vtlt_seeds(j)*(duration(j)-duration(j-1));
    end
end
fL(1) = fL(2);

gL = PL - fL;
if plotFigures  
   figure 
    subplot(1,2,1);
    plot(PL,'linewidth',1);
    hold on;
    plot(fL,'linewidth',1);
    xlabel('L(MW)');
    ylabel('Price (€/MW \cdot Month)');
    legend('P(L)','f(L)')   
    
    subplot(1,2,2);
    yyaxis left
    plot(gL,'linewidth',1);
    xlabel('L(MW)');
    ylabel('Price (€/MW \cdot Month)');
    hold on
    yyaxis right
    plot(reliablity(L_supplied),'linewidth',1);
    ylabel('Reliablity');
    ylim([0,1.01]);
    legend('g(L)','Reliability')
end

%% cost, price and valuation
cost = zeros(L0,1);
valueFunction = zeros(L0,1);

for  L = 1:L0
    cost(L) = c(L,welfareIdx(L))*duration(L);
    valueFunction(L) = m*duration(L)^((1+1/lambda)) + (a0-b*L)*duration(L);
end
price = Hr.*PL;
valueFunction = Hr.*valueFunction;

if plotFigures
    figure
    plot(L_supplied,cost,L_supplied,valueFunction,L_supplied,price,'linewidth',1);   
    xlabel('L(MW)');
    ylabel('Price (€/MW)');
    legend('Cost','Valuation','Price');
end
