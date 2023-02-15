clear all
close all
load('silver.dat')
br = 9.32; %background radiation in bq
co = 37;
silver(:,2) = silver(:,2) - br * 3; %Activity from silver samples
silverlog1 = log(silver(:,2)); %logaritmized activity
figure(1);
scatter(silver(:,1), silverlog1) %plot of logaritmic activity over time

fit1 = polyfit(silver(co:end, 1), log(silver(co:end,2)), 1); %gives ln(A) = ln(A0) - labda*t

%ln(A0) = 7.2655 - 6.7577, -lambda*t = -0.0052 - -0.0045
silver1 = exp(silver(:,1) * fit1(1)) * exp(fit1(2)); %activity for long decay silver A = A0*e^(-lambda*t)
figure(2);
plot(silver(:,1), (silver1));

silver2 = (silver(:,2)) - (silver1); %approximation for quick decay silver
figure(3); 
scatter(silver(:,1), (silver2))


fit2 = polyfit(silver(1:co,1), log(silver2(1:co)),1); %fit for quick decay silver
%ln(A0) = 1.5415, -lambda*t = -0.0114
A01 = exp(fit1(2));
A02 = exp(fit2(2));
lambda1 = -fit1(1);
lambda2 = -fit2(1);
T1 = (log(2) / lambda1);
T2 = (log(2) / lambda2);


%% del 2
clear all
close all
t = 1E-3; %thickness of protective plates in m
conf95 = 12.71; %number of standard deviations from t-table for two tailed confidence interval of 95%
conf60 = 1.376; %number of standard deviations from t-table for two tailed confidence interval of 60%
conf70 = 1.963; %number of standard deviations from t-table for two tailed confidence interval of 70%
conf80 = 3.078; %number of standard deviations from t-table for two tailed confidence interval of 80%
Ma = 6.022E23; %Avogadros number
Np = 17836 - 9.32*120; %number of decays without plates
Na = 18253 - 9.32*120; %number of decays with aluminum plates
Nc = 2389 - 9.32*120; %number of decays with cadmium plates
nc = 8.65E3*Ma/0.112411; %number of cadmium particles per m^3
na = 2.7E3*Ma/0.02698; %number of aluminum particles per m^3
csa = -log(Na/Np)/(na*t); %cross section for aluminum m^2/particle
csc = -log(Nc/Np)/(nc*t); %cross section for cadmium m^2/particle

stda = sqrt(Na); %standard deviation of decays with aluminium plates
stdc = sqrt(Nc); %standard deviation of decays with cadmium plates
stdp = sqrt(Np); %standard deviation of decays with no plates

Np2max = Np + conf95*stdp; %Maximum decays no plates 95%conf
Np2min = Np - conf95*stdp; %Minimum decays no plates 95%conf
Na2min = Na - conf95*stda; %Minimum decays aluminium plates 95%conf
Na2max = Na + conf95*stda; %Maximum decays aluminium plates 95%conf
Nc2min = Nc - conf95*stdc; %Minimum decays cadmium plates 95%conf
Nc2max = Nc + conf95*stdc; %Maximum decays cadmium plates 95%conf

csa2min = -log(Na2max/Np2min)/(na*t); %minimizes cross section for aluminium
csa2max = -log(Na2min/Np2max)/(na*t); %maximizes cross section for aluminium
csc2min = -log(Nc2max/Np2min)/(nc*t); %minimizes cross section for cadmium
csc2max = -log(Nc2min/Np2max)/(nc*t); %maximizes cross section for aluminium







