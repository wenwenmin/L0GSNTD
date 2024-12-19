function plt0=plotplt(trigger,acc,rdx,NMIs)

mess={};
plt0=0;
num=length(trigger);
for i=1:num
    mess=pltplt(num,trigger(i),mess);
end

accmean=mean(acc)';
accstd=std(acc)';
rdxmean=mean(rdx)';
rdxstd=std(rdx)';
nmimean=mean(NMIs)';
nmistd=std(NMIs)';

x=1:length(trigger);
color= parula(12);


figure(1)
% GO=bar(x,accmean,0.3);

m = length(color);
for k = 1:length(accmean)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;
    h = bar(x(k), accmean(k),0.3);
    set(h,'FaceColor',color(i,:));
end

hE = errorbar(x, accmean,accstd,'LineWidth',2);
ylabel('ACC','FontSize',35);
set(gca,'xticklabel',mess,'FontSize',35);
xticks(1:length(trigger));
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
prettyAxes().gbase2()
hE.Color = [0 0 0];                            
hE.LineStyle = 'none'; 

figure(2)
for k = 1:length(rdxmean)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;
    h = bar(x(k), rdxmean(k),0.3);
    set(h,'FaceColor',color(i,:));
end
hE = errorbar(x, rdxmean,rdxstd,'LineWidth',2);

ylabel('ARI','FontSize',35);
xticks(1:length(trigger));
set(gca,'xticklabel',mess,'FontSize',35);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
prettyAxes().gbase2()
hE.Color = [0 0 0];                            
hE.LineStyle = 'none'; 

figure(3)
for k = 1:length(nmimean)
    i = mod(k-1,m); %%i is remainder after division of k-1 by m
    i = i+1;
    h = bar(x(k), nmimean(k),0.3);
    set(h,'FaceColor',color(i,:));
end
hE = errorbar(x, nmimean,nmistd,'LineWidth',2);
ylabel('NMI','FontSize',35);
xticks(1:length(trigger));
set(gca,'xticklabel',mess,'FontSize',35);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';

prettyAxes().gbase2()

hE.Color = [0 0 0];                            
hE.LineStyle = 'none';  


set(hE, 'LineStyle', 'none', 'Color', 'k','LineWidth', 1.2)
end

function mess=pltplt(num,trigger,mess)
figure(1);
for i=1:num
    if(trigger==1)
    mes='$\ell_0$-SNTD';
    elseif(trigger==2)
    mes='$\ell_1$-GSNTD';
    elseif(trigger==3)
    mes='$\ell_0$-GSNTD';
    end
    hold on;
end
ls=length(mess);
mess{1+ls}=mes;

end

