set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
load('tsp_cnr_sweep_new.mat')
I_total = zeros(length(CNR_sweep),1);
I_total_fdcomm = zeros(length(CNR_sweep),1);
I_total_radar = zeros(length(CNR_sweep),1);
I_total_ULcomm = zeros(length(CNR_sweep),1);
I_total_DLcomm = zeros(length(CNR_sweep),1);
% 
% for ii = 1:length(CNR_sweep)
%     I_total(ii) = fdcomm_save{ii}.I_total_op(end)+sum(radar_save{ii}.MI_radar);
%     I_total_fdcomm(ii) = fdcomm_save{ii}.I_total_op(end);
%     I_total_radar(ii) = sum(radar_save{ii}.MI_radar);
%     I_total_ULcomm(ii) = sum(fdcomm_save{ii}.MI_UL(:));
%     I_total_DLcomm(ii) = sum(fdcomm_save{ii}.MI_DL(:));
% end
for ii = 1:length(CNR_sweep)
    I_total(ii) = fdcomm_save{ii}.I_total(end);
    I_total_radar(ii) = fdcomm_save{ii}.I_radar_op(end);
    I_total_ULcomm(ii) = fdcomm_save{ii}.I_UL_op(end);
    I_total_DLcomm(ii) = fdcomm_save{ii}.I_DL_op(end);
    I_total_fdcomm(ii) = I_total_ULcomm(ii) + I_total_DLcomm(ii);
end
figure();
plot(CNR_sweep, I_total,'m--o','LineWidth',1.5)
hold on
plot(CNR_sweep, I_total_radar,'r-+','LineWidth',1.5)
hold on
plot(CNR_sweep, I_total_fdcomm,'g:*','LineWidth',1.5)
hold on
plot(CNR_sweep, I_total_ULcomm,'k-.x','LineWidth',1.5)
hold on
plot(CNR_sweep, I_total_DLcomm,'b:^','LineWidth',1.5)
hold off
grid on
axis square
xlabel('$\mathrm{CNR}_{\textrm{r}}\textrm{ (dB)}$','FontSize',12,'Interpreter','latex');
%ylabel('$\mathrm{MI}\textrm{ (bits per channel use)}$','Interpreter','latex','FontSize',12)
ylabel('Weighted MI (bits per channel use)', 'FontSize',12)
legend('total MI','radar MI', 'FD MI','UL MI', 'DL MI','FontSize',10)
legend('Location','best');
legend('boxoff');
print('-depsc2','-r150','CNR_sweep.eps')