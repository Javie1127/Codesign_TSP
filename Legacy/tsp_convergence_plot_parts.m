set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
%% Convergence Iteration vs total rate 
%load('convergence_complete_v2.mat')
load('tsp_convergence_total.mat')
iter = 1:radar.ell_max; 
for ii = 1:length(SNR_rtr)
    figure()
    plot(iter,I_total_ew_ni(:,ii),'r-','LineWidth',1.5);
    hold on
    plot(iter,I_total_uew_ni(:,ii),'c--','LineWidth',1.5);
    hold on
    plot(iter,I_total_ew_ri(:,ii),'b-d','LineWidth',1.5);
    hold on
    plot(iter,I_total_uew_ri(:,ii),'k-.d','LineWidth',1.5);
    grid on
    xlabel('Iteration Index','FontSize',10)
    %ylabel('$\Xi_{\textrm{MSE}}$','Interpreter','latex','FontSize',10)
    ylabel('$\mathrm{I}_{\textrm{CWSM}} \textrm{(bits per channel use)}$','Interpreter','latex','FontSize',10);
    %legend(['\omega = [',num2str(omega_1),']';'\omega = [',num2str(omega_2),']'])
    %legend({strcat('\omega = [',num2str(omega_1'),']'),strcat('\omega = [',num2str(omega_2'),']')});
    legend('Uniform weights initialization method 1','Non-uniform weights initialization method 1',...
        'Uniform weights initialization method 2', 'Non-uniform weights initialization method 2', 'FontSize',12);
    legend('Location','best');
    legend('boxoff');
    %ylim([8,15]);
    xlim ([1,30]);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',10)
end

% axes('position',[.35 .575 .25 .25], 'NextPlot', 'add')
% indexOfInterest1 = (I_total_op_ew_ni < 110.0823) & (I_total_op_ew_ni > 110.0736);
% indexOfInterest2 = (I_total_op_uew_ni < 110.0757) & (I_total_op_uew_ni > 110.0682);
% indexOfInterest3 = (iter>=7) & (iter<=12);
% plot(iter(indexOfInterest3), I_total_op_ew_ni(indexOfInterest1));
% plot(iter(indexOfInterest3), I_total_op_uew_ni(indexOfInterest2));
print('-depsc2','-r150','convergence.eps')
