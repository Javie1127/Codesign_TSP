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
figure()
C = {'k','b','r','g','y'};
% Legend=cell(N,1);
%  for iter=1:N
%    Legend{iter}=strcat('Your_Data number', num2str(iter));
%  end
%  legend(Legend)
%for ii = 1:length(SNR_rtr)
for ii = 1:length(SNR_rtr)
    plot(iter,I_total_ew_ni(:,ii),'color',C{ii},'LineStyle','-','LineWidth',1.5);
    hold on
    plot(iter,I_total_ew_ri(:,ii),'color',C{ii},'LineStyle','-.','LineWidth',1.5);
    hold on
    %legend('Initialization 1','Initialization 2')
%     plot(iter,I_UL_ew_ni(:,ii),'k:','LineWidth',1.5);
%     hold on
%     I_FD_ew_ni_ii = I_UL_ew_ni(:,ii) + I_DL_ew_ni(:,ii);
%     plot(iter,I_FD_ew_ni_ii,'b--','LineWidth',1.5);
%     hold on
%     plot(iter,I_radar_ew_ni(:,ii),'g:','LineWidth',1.5);
%     hold on
%     I_FD_ew_ri_ii = I_UL_ew_ri_ii(:,ii) + I_DL_ew_ri_ii(:,ii); 
%     %plot(iter,I_UL_ew_ri(:,ii),'k-.*','LineWidth',1.5);
%     plot(iter,I_FD_ew_ri_ii,'k-','LineWidth',1.5);
%     hold on
% %     plot(iter,I_DL_ew_ri(:,ii),'g--','LineWidth',1.5);
% %     hold on
%     plot(iter,I_radar_ew_ri(:,ii),'m-.','LineWidth',1.5);
%     hold off
end
hold off
grid on
xlabel('Iteration Index','FontSize',12)
ylabel('Weigthed MI (bits per channel use)','FontSize',12)
legend('initialization 1 $\textrm{SNR}_{\textrm{r}} = -5 \textrm{ dB}$',' initialization 2 $\textrm{SNR}_{\textrm{r}} = -5 \textrm{ dB}$',...
' initialization 1 $\textrm{SNR}_{\textrm{r}} = 0 \textrm{ dB}$',' initialization 2 $\textrm{SNR}_{\textrm{r}} = 0 \textrm{ dB}$',' initialization 1 $\textrm{SNR}_{\textrm{r}} = 5 \textrm{ dB}$',...
' initialization 2 $\textrm{SNR}_{\textrm{r}} = 5 \textrm{ dB}$',' initialization 1 $\textrm{SNR}_{\textrm{r}} = 10 \textrm{ dB}$',' initialization 2 $\textrm{SNR}_{\textrm{r}} = 10 \textrm{ dB}$','Interpreter','latex','FontSize',10);
legend('Location','best','NumColumns',2);
legend('boxoff');
ylim([11,18]);
xlim ([1,30]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',12)
% axes('position',[.35 .575 .25 .25], 'NextPlot', 'add')
% indexOfInterest1 = (I_total_op_ew_ni < 110.0823) & (I_total_op_ew_ni > 110.0736);
% indexOfInterest2 = (I_total_op_uew_ni < 110.0757) & (I_total_op_uew_ni > 110.0682);
% indexOfInterest3 = (iter>=7) & (iter<=12);
% plot(iter(indexOfInterest3), I_total_op_ew_ni(indexOfInterest1));
% plot(iter(indexOfInterest3), I_total_op_uew_ni(indexOfInterest2));
print('-depsc2','-r150','tsp_convergence_snr.eps')
