set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
%% 
%load('tsp_fd_vs_hd.mat')
%load(')
% load('tsp_fd_vs_hd_v5_wc.mat');
% DLcomm_wc_save_plot = DLcomm_save_wc;
% load('tsp_fd_vs_hd_v5.mat')
%load('tsp_fd_vs_hd_v6_wc.mat')
load('tsp_fd_vs_hd_wc_v7.mat')
% fdcomm_save_plot = fdcomm_save;
% DLcomm_save_plot = DLcomm_save;
% ULcomm_save_plot = ULcomm_save;
%SNR_rtr = SNR_rtr(1:end-1);
I_total_fdcomm = zeros(length(SNR_rtr),1);
I_total_fdcomm_wc = zeros(length(SNR_rtr),1);
I_total_ULcomm = zeros(length(SNR_rtr),1);
I_total_DLcomm = zeros(length(SNR_rtr),1);
I_total_DLcomm_wc = zeros(length(SNR_rtr),1);
for ii = 1:length(SNR_rtr)
    I_total_fdcomm(ii) = fdcomm_save{ii}.I_total_op(end);
    I_total_fdcomm_wc(ii) = fdcomm_wc_save{ii}.I_total_op(end);
    I_total_ULcomm(ii) = ULcomm_save{ii}.I_total_op(end);
    I_total_DLcomm(ii) = DLcomm_save{ii}.I_total_op(end);
    I_total_DLcomm_wc(ii) = DLcomm_save_wc{ii}.I_total_op(end);
end
figure()
plot(SNR_rtr,I_total_fdcomm,'m--o','LineWidth',1.5);
hold on
plot(SNR_rtr,I_total_fdcomm_wc,'b-.','LineWidth',1.5);
hold on
plot(SNR_rtr,I_total_ULcomm,'r-+','LineWidth',1.5);
hold on
plot(SNR_rtr,I_total_DLcomm,'k-.x','LineWidth',1.5);
hold on
plot(SNR_rtr,I_total_DLcomm_wc,'g:*','LineWidth',1.5);
hold off
grid on
xlabel('$\mathrm{SNR}_{\textrm{r}}\textrm{ (dB)}$','FontSize',12,'Interpreter','latex');
%ylabel('$\Xi_{\textrm{MSE}}$','Interpreter','latex','FontSize',10)
%ylabel('$\mathrm{MI} \textrm{(bits per channel use)}$','Interpreter','latex','FontSize',12);
ylabel('Weighted sum MI (bits per channel use)', 'FontSize',12)
%ylabel('Weighted MI (bits per channel use)','FontSize',12)
%xlim([-12,10])
xticks(SNR_rtr);
a = get(gca,'XTickLabel');
for ii = 2:2:length(a)
    a{ii} = ' ';
end
set(gca,'XTickLabel',a,'FontName','Times','fontsize',12)
%legend(['\omega = [',num2str(omega_1),']';'\omega = [',num2str(omega_2),']'])
%legend({strcat('\omega = [',num2str(omega_1'),']'),strcat('\omega = [',num2str(omega_2'),']')});
legend('FD+radar with cooperation','FD+radar without cooperation','UL+radar', 'DL+radar with cooperation','DL+radar without cooperation', 'FontSize',12);
legend('Location','northwest');
legend('boxoff');
print('-depsc2','-r150','fd_vs_hd.eps')