set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
%%
load('Pd_vs_vu_SNR_Btr_stpsz2.mat','PD','vu','SNR_Btr')
ii = 0;
figure()
% LegendCell = cell(length(SNR_Btr),1);
linestyle = ["-","--","-.",":","-"];
linecolor = 'bkrgc';
edgecolor = 'rmybk';
linemarker = ["none",".","*","x",'s'];
for snr_Btr = SNR_Btr
    ii = ii + 1;
    plot(vu,PD(ii,:),'Color',linecolor(ii),'Marker',char(linemarker(ii)),'MarkerSize',8,'MarkerEdgeColor',edgecolor(ii),'LineStyle',char(linestyle(ii)),...
        'LineWidth',1.5,'DisplayName',['$\mathrm{SNR}_{\textrm{Btr}}=$',num2str(snr_Btr),'\textrm{ dB}']);
    hold on
    %LegendCell{ii}=['$\textrm{SNR}_{\textrm{Btr}}=$',num2str(snr_Btr),'\text{ dB}'];
end
hold off
grid on
xlabel('Detection Threshold \nu','FontSize',10)
ylabel('$\mathrm{P}_\mathrm{d}$','Interpreter','latex','FontSize',10);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',10)
%legend(['\omega = [',num2str(omega_1),']';'\omega = [',num2str(omega_2),']'])
%legend({strcat('\omega = [',num2str(omega_1'),']'),strcat('\omega = [',num2str(omega_2'),']')});
%legend(LegendCell,'Interpreter','latex');
legend('Interpreter','latex');
legend('Location','northeast');
legend('boxoff');
print('-depsc2','-r150','PD_vs_nu.eps')