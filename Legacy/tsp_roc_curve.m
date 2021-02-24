%-------------------- Plot Recever operating characteristic curve---
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
load('tsp_Pd_vs_radar_SNR_v2.mat');
figure()
linestyle = ["-.","--","-.",":","--",'--'];
linecolor = 'bkrgcm';
edgecolor = 'rmybkc';
linemarker = ["none",".","*","x",'s','none'];
for ii = 1:length(SNR_rtr)
    plot(PFA_co(ii,:),PD_co(ii,:),'Color',linecolor(ii),'Marker',char(linemarker(ii)),'MarkerSize',8,'MarkerEdgeColor',edgecolor(ii),'LineStyle',char(linestyle(ii)),...
        'LineWidth',1.5,'DisplayName',['$\mathrm{SNR}_{\textrm{rtr}}=$',num2str(SNR_rtr(ii)),'\textrm{ dB}']);
    hold on
end
hold off
grid on
xlabel('$\mathrm{P}_{\textrm{fa}}$','FontSize',10,'Interpreter','latex');
ylabel('$\mathrm{P}_\mathrm{d}$','Interpreter','latex','FontSize',10);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',10)
%legend(['\omega = [',num2str(omega_1),']';'\omega = [',num2str(omega_2),']'])
%legend({strcat('\omega = [',num2str(omega_1'),']'),strcat('\omega = [',num2str(omega_2'),']')});
%legend(LegendCell,'Interpreter','latex');
legend('Interpreter','latex');
legend('Location','southeast');
legend('boxoff');
print('-depsc2','-r150','PD_vs_vu_v1.eps')