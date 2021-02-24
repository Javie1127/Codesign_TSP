set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
%% Convergence Iteration vs total rate 
%load('convergence_complete_v2.mat')
load('tsp_convergence_new.mat')
iter = 1:radar.ell_max; 
%for ii = 1:length(SNR_rtr)
for ii = 
    figure()
    plot(iter,I_total_ew_ni(:,ii),'r--','LineWidth',1.5);
    hold on
%     plot(iter,I_UL_ew_ni(:,ii),'k:','LineWidth',1.5);
%     hold on
    I_FD_ew_ni_ii = I_UL_ew_ni(:,ii) + I_DL_ew_ni(:,ii);
    plot(iter,I_FD_ew_ni_ii,'b--','LineWidth',1.5);
    hold on
    plot(iter,I_radar_ew_ni(:,ii),'g:','LineWidth',1.5);
    hold on
    plot(iter,I_total_ew_ri(:,ii),'b-.','LineWidth',1.5);
    hold on
    I_FD_ew_ri_ii = I_UL_ew_ri_ii(:,ii) + I_DL_ew_ri_ii(:,ii); 
    %plot(iter,I_UL_ew_ri(:,ii),'k-.*','LineWidth',1.5);
    plot(iter,I_FD_ew_ri_ii,'k-','LineWidth',1.5);
    hold on
%     plot(iter,I_DL_ew_ri(:,ii),'g--','LineWidth',1.5);
%     hold on
    plot(iter,I_radar_ew_ri(:,ii),'m-.','LineWidth',1.5);
    hold off
    grid on
    xlabel('Iteration Index','FontSize',12)
    %ylabel('$\Xi_{\textrm{MSE}}$','Interpreter','latex','FontSize',10)
    %ylabel('$\mathrm{I}_{\textrm{CWSM}} \textrm{(bits per channel use)}$','Interpreter','latex','FontSize',10);
    ylabel('Weigthed MI (bits per channel use)','FontSize',12)
    %legend(['\omega = [',num2str(omega_1),']';'\omega = [',num2str(omega_2),']'])
    %legend({strcat('\omega = [',num2str(omega_1'),']'),strcat('\omega = [',num2str(omega_2'),']')});
%     legend('$I_{\textrm{CWSM}}$ initialization 1','$I_{\textrm{u}}$ initialization 1',...
%         '$I_{\textrm{d}}$ initialization 1','$I_{\textrm{r}}$ initialization 1','$I_{\textrm{CWSM}}$ initialization 2','$I_{\textrm{u}}$ initialization 2',...
%         '$I_{\textrm{d}}$ initialization 2','$I_{\textrm{r}}$ initialization 2','Interpreter','latex','FontSize',12);
        legend('$I_{\textrm{CWSM}}$ initialization 1','$I_{\textrm{fd}}$ initialization 1',...
        '$I_{\textrm{r}}$ initialization 1','$I_{\textrm{CWSM}}$ initialization 2','$I_{\textrm{fd}}$ initialization 2',...
        '$I_{\textrm{r}}$ initialization 2','Interpreter','latex','FontSize',11);
    legend('Location','best','NumColumns',2);
    legend('boxoff');
    ylim([0,18]);
    xlim ([1,30]);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',12)
end

% axes('position',[.35 .575 .25 .25], 'NextPlot', 'add')
% indexOfInterest1 = (I_total_op_ew_ni < 110.0823) & (I_total_op_ew_ni > 110.0736);
% indexOfInterest2 = (I_total_op_uew_ni < 110.0757) & (I_total_op_uew_ni > 110.0682);
% indexOfInterest3 = (iter>=7) & (iter<=12);
% plot(iter(indexOfInterest3), I_total_op_ew_ni(indexOfInterest1));
% plot(iter(indexOfInterest3), I_total_op_uew_ni(indexOfInterest2));
print('-depsc2','-r150','tsp_convergence_10dB.eps')
