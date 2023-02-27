close all
clear all

load silver.dat

conf95 = 6.314; %number of standard deviations from t-table for two tailed confidence interval of 95%

silver(:,2) = silver(:,2)-3*9.32;

p1 = polyfit(silver(80:146,1),log(silver(80:146,2)),1)
pfit1 = exp(silver(:,1).*p1(1))*exp(p1(2));


s2 = silver(:,2)-pfit1;
p2 = polyfit(silver(1:30,1),log(s2(1:30)),1)
pfit2 = exp(silver(:,1)*p2(1))*exp(p2(2));

sumFit = pfit1 + pfit2;

plot(silver(:,1),silver(:,2), '.',silver(:,1),sumFit,'--')
hold on
plot(silver(:,1),pfit1)
plot(silver(:,1),pfit2)
title('Decay as a function of time.')
xlabel('Time (s)') 
ylabel('Activity (decay/s)') 
legend('Activity of both isotopes', 'Mean decay of both isotopes','Fast isotope', 'Slow isotope')
hold off

figure(2);
plot(silver(:,1),log(silver(:,2)), '.',silver(:,1),log(sumFit),'--')
hold on
plot(silver(:,1),log(pfit1))
plot(silver(:,1),log(pfit2))
title('Decay as a function of time logarithmised.')
xlabel('Time (s)') 
ylabel('Activity (decay/s)')
legend('Activity of both isotopes', 'Mean decay of both isotopes','Fast isotope', 'Slow isotope')

halfTime1 = -log(2)/p1(1)
halfTime2 = -log(2)/p2(1)

AmountAtStart1 = pfit1(1)/silver(1,2) % långsam
AmountAtStart2 = pfit2(1)/silver(1,2) % snabb
AmountTotal = AmountAtStart1 + AmountAtStart2


%%


clear all

load silver.dat

NoP = 17836-120*9.32;
AlP = 18253-120*9.32;
CaP = 2389-120*9.32;

DAl = 2.7*10^6;
MAl = 26.98;
DCa = 8.65*10^6;
MCa = 112.411;


avo = 6.022*10^23;
nAl = DAl*avo/MAl;  
nCa = DCa*avo/MCa;

CSAl = log(NoP/AlP)/(nAl*0.001*10^-28)
CSCa = log(NoP/CaP)/(nCa*0.001*10^-28)

SDAl = sqrt(AlP)
SDCa = sqrt(CaP)
SDNP = sqrt(NoP)

maxAl = AlP+6.314*SDAl 
minAl = AlP-6.314*SDAl
maxNS = NoP+6.314*SDNP
minNS = NoP-6.314*SDNP
maxCa = CaP+6.314*SDCa
minCa = CaP-6.314*SDCa

maxCSAl = log(maxNS/minAl)/(nAl*0.001*10^-28)
maxCSCa = log(maxNS/minCa)/(nCa*0.001*10^-28)

minCSAl = log(minNS/maxAl)/(nAl*0.001*10^-28)
minCSCa = log(minNS/maxCa)/(nCa*0.001*10^-28)

%Np2max = Np + conf95*stdp; %Maximum decays no plates 95%conf
