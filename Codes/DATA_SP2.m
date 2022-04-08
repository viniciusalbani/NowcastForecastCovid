%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reading Daily Cases and Deaths in Texas

aux = importdata('Dados-covid-19-estado_20211119.txt');
dailyCasesDeaths = aux.data;
datesCD = string(aux.textdata(2:end,1));
datesCD = datetime(datesCD,'InputFormat','dd/MM/yyyy');

t_span = datesCD;

dataA = dailyCasesDeaths(1:end,:);
for jj=7:size(dataA,1)
for ii = 1:size(dataA,2)
dataA(jj,ii) = mean(dailyCasesDeaths(jj-6:jj,ii));
end
end

dataB = dataA;
for jj=7:size(dataA,1)
for ii = 1:size(dataA,2)
dataB(jj,ii) = mean(dataA(jj-6:jj,ii));
end
end

Population = 46649132;      


dailyData = [dataB(:,2),dataB(:,3)];
data2 = dailyData(1:end-3,:);
t_span = t_span(1:size(data2,1));
t_actual = 0:length(t_span);
