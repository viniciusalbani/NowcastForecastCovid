clear all; close all; %clc;

%%%% Chicago

SummaryChicago1 = zeros(25,6);
SummaryChicago2 = zeros(25,6);
CI90Chicago1= zeros(15,6);
CI90Chicago2= zeros(15,6);
time=[10,15,30,45,60,90];
for zz = 1:6
load(['data_Chicago_20220101_',num2str(time(zz)),'dias.mat']);

NSamples = size(Err1,1);

aux = sort(Err1);
aux2 = round(0.25*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI50 = [min(aux);max(aux)];

aux = sort(Err1);
aux2 = round(0.15*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI90 = [min(aux);max(aux)];

%%%% Summary Statistics
aux = [min(Err1);CI50(1,:);median(Err1);CI50(2,:);max(Err1)];
SummaryChicago1(:,zz) = reshape(aux,1,25)';

aux = [CI90(1,:);median(Err1);CI90(2,:)];
CI90Chicago1(:,zz) = reshape(aux,1,15)';
%%%%%%%%%%
aux = sort(Err2);
aux2 = round(0.25*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI50 = [min(aux);max(aux)];

aux = sort(Err2);
aux2 = round(0.15*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI90 = [min(aux);max(aux)];

%%%% Summary Statistics
aux = [min(Err2);CI50(1,:);median(Err2);CI50(2,:);max(Err2)];
SummaryChicago2(:,zz) = reshape(aux,1,25)';

aux = [CI90(1,:);median(Err1);CI90(2,:)];
CI90Chicago2(:,zz) = reshape(aux,1,15)';
if zz==1
aux = sort(coefMR(:,1));
aux2 = round(0.15*size(coefMR,1));
aux = aux(aux2+1:end-aux2,:);
CI70A = [min(aux),max(aux)];
aux = sort(coefMR(:,2));
aux2 = round(0.15*size(coefMR,2));
aux = aux(aux2+1:end-aux2,:);
CI70B = [min(aux),max(aux)];
aux = sort(coefMR(:,3));
aux2 = round(0.15*size(coefMR,3));
aux = aux(aux2+1:end-aux2,:);
CI70C = [min(aux),max(aux)];
aux = sort(coefMean);
aux2 = round(0.15*length(coefMean));
aux = aux(aux2+1:end-aux2,:);
CI70D = [min(aux),max(aux)];
AUX2 = [median(coefMean),CI70D,median(coefMR(:,1)),CI70A,median(coefMR(:,2)),CI70B,median(coefMR(:,3)),CI70C];
disp(num2str(AUX2))
end
end

H = [100,100,700,400];

tt = [10,15,30,45,60,90];
CI90Chicago1 = 100*CI90Chicago1;

figure
hold on
grid on
box on
title('Chicago - US')
plot(tt,CI90Chicago1(2,:),'-sb','LineWidth',2)
plot(tt,CI90Chicago1(5,:),'-og','LineWidth',2)
plot(tt,CI90Chicago1(8,:),'-+r','LineWidth',2)
plot(tt,CI90Chicago1(11,:),'-^m','LineWidth',2)
plot(tt,CI90Chicago1(14,:),'-*c','LineWidth',2)
legend('Average','Linear Regression','Fourier','Fourier w. Noise','Mean-Reverting','Location','NorthWest')
xlabel('Number of days ahead')
ylabel('Median Prediction Error (%)')
xticks([0,tt,tt(end)+10])
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'Chicago_MEDIANErr1.fig');
print('-dpng','Chicago_MEDIANErr1');

AUX = [CI90Chicago1(2:3:end,1),CI90Chicago1(1:3:end,1),CI90Chicago1(3:3:end,1),CI90Chicago1(2:3:end,2),CI90Chicago1(1:3:end,2),CI90Chicago1(3:3:end,2),CI90Chicago1(2:3:end,3),CI90Chicago1(1:3:end,3),CI90Chicago1(3:3:end,3),...
CI90Chicago1(2:3:end,4),CI90Chicago1(1:3:end,4),CI90Chicago1(3:3:end,4),CI90Chicago1(2:3:end,5),CI90Chicago1(1:3:end,5),CI90Chicago1(3:3:end,5),CI90Chicago1(2:3:end,6),CI90Chicago1(1:3:end,6),CI90Chicago1(3:3:end,6)]/100;