set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
%%
load('tsp_Pd_vs_radar_SNR_v3.mat')
for ii = 2
    figure()
    linestyle = ["-","--","-.",":","-","--"];
    linecolor = 'bkrgcm';
    edgecolor = 'rmybkc';
    %linemarker = ["+",".","*","x","s","o"];
    
    [PFA_co_interp,PD_co_interp] = smooth_data(PFA_co(ii,:),PD_co(ii,:));
    %plot(PFA_co_interp,PD_co_interp,'Color',linecolor(1),'Marker',char(linemarker(1)),'MarkerSize',5,'MarkerEdgeColor',edgecolor(1),'LineStyle',char(linestyle(1)),...
    plot(PFA_co_interp,PD_co_interp,'Color',linecolor(1),'LineStyle',char(linestyle(1)),...
        'LineWidth',1,'DisplayName','Proposed with cooperation');
    hold on
    [PFA_co_wc_interp,PD_co_wc_interp] = smooth_data(PFA_co_wc(ii,:),PD_co_wc(ii,:));
    %plot(PFA_co_wc_interp,PD_co_wc_interp,'Color',linecolor(2),'Marker',char(linemarker(2)),'MarkerSize',5,'MarkerEdgeColor',edgecolor(2),'LineStyle',char(linestyle(2)),...
    plot(PFA_co_wc_interp,PD_co_wc_interp,'Color',linecolor(2),'LineStyle',char(linestyle(2)),...
        'LineWidth',1,'DisplayName','Proposed without cooperation');
    hold on
    [PFA_uc_interp,PD_uc_interp] = smooth_data(PFA_uc(ii,:),PD_uc(ii,:));
    %plot(PFA_uc_interp,PD_uc_interp,'Color',linecolor(3),'Marker',char(linemarker(3)),'MarkerSize',5,'MarkerEdgeColor',edgecolor(3),'LineStyle',char(linestyle(3)),...
    plot(PFA_uc_interp,PD_uc_interp,'Color',linecolor(3),'LineStyle',char(linestyle(3)),...
        'LineWidth',1,'DisplayName','Uncoded with cooperation');
    hold on
    [PFA_uc_wc_interp,PD_uc_wc_interp] = smooth_data(PFA_uc_wc(ii,:),PD_uc_wc(ii,:));
    %plot(PFA_uc_wc_interp,PD_uc_wc_interp,'Color',linecolor(4),'Marker',char(linemarker(4)),'MarkerSize',5,'MarkerEdgeColor',edgecolor(4),'LineStyle',char(linestyle(4)),...
    plot(PFA_uc_wc_interp,PD_uc_wc_interp,'Color',linecolor(4),'LineStyle',char(linestyle(4)),... 
        'LineWidth',1,'DisplayName','Uncoded without cooperation');
    hold on
    [PFA_rc_interp,PD_rc_interp] = smooth_data(PFA_rc(ii,:),PD_rc(ii,:));
    %plot(PFA_rc_interp,PD_rc_interp,'Color',linecolor(5),'Marker',char(linemarker(5)),'MarkerSize',5,'MarkerEdgeColor',edgecolor(5),'LineStyle',char(linestyle(5)),...
    plot(PFA_rc_interp,PD_rc_interp,'Color',linecolor(5),'LineStyle',char(linestyle(5)),...
        'LineWidth',1,'DisplayName','Randomly coded with cooperation');
    hold on
    [PFA_rc_wc_interp,PD_rc_wc_interp] = smooth_data(PFA_rc_wc(ii,:),PD_rc_wc(ii,:));
    %plot(PFA_rc_wc_interp,PD_rc_wc_interp,'Color',linecolor(6),'Marker',char(linemarker(6)),'MarkerSize',5,'MarkerEdgeColor',edgecolor(6),'LineStyle',char(linestyle(6)),...
    plot(PFA_rc_wc_interp,PD_rc_wc_interp,'Color',linecolor(6),'LineStyle',char(linestyle(6)),...
        'LineWidth',1,'DisplayName','Randomly coded without cooperation');
    hold off
    grid on
    xlabel('$\mathrm{P}_\mathrm{fa}$','Interpreter','latex','FontSize',12)
    ylabel('$\mathrm{P}_\mathrm{d}$','Interpreter','latex','FontSize',12);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','FontSize',12)
    axis square
    legend('FontSize',8,'Location','southeast');
    legend('boxoff');
end
% LegendCell = cell(length(SNR_Btr),1);

%legend(['\omega = [',num2str(omega_1),']';'\omega = [',num2str(omega_2),']'])
%legend({strcat('\omega = [',num2str(omega_1'),']'),strcat('\omega = [',num2str(omega_2'),']')});
%legend(LegendCell,'Interpreter','latex');
%legend('Interpreter','latex');
print('-depsc2','-r150','tsp_ROC.eps')