set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
load('tsp_UL_SNR_avg.mat');
I_total_radar_new = I_total_radar;
I_total_UL_new = I_total_UL;
I_total_DL_new = I_total_DL;
I_total_FD_new = I_total_FD;
%load('tsp_UL_SNR.mat','I_total','I_total_radar','I_total_DL','I_total_FD','I_total_UL');
%I_total_radar_new = [I_total_radar_new,I_total_radar];
%I_total_UL_new = [I_total_UL_new,I_total_UL];
%I_total_DL_new = [I_total_DL_new,I_total_DL];
%I_total_FD_new = [I_total_FD_new,I_total_FD];
I_total_radar_avg = mean(I_total_radar_new,2);
I_total_UL_avg = mean(I_total_UL_new,2);
I_total_DL_avg = mean(I_total_DL_new,2);
I_total_FD_avg = mean(I_total_FD_new,2);
I_total_avg = I_total_radar_avg + I_total_UL_avg + I_total_DL_avg;
%load('tsp_UL_SNR_5_10dB.mat','')
figure()
% I_total = zeros(length(UL_SNR_sweep),1);
% I_total_fdcomm = zeros(length(UL_SNR_sweep),1);
% I_total_radar = zeros(length(UL_SNR_sweep),1);
% I_total_ULcomm = zeros(length(UL_SNR_sweep),1);
% I_total_DLcomm = zeros(length(UL_SNR_sweep),1);
% for ii = 1:length(UL_SNR_sweep)
%     I_UL_weighted = fdcomm_save{ii}.MI_UL.*fdcomm_save{ii}.alpha_UL;
%     I_DL_weighted = fdcomm_save{ii}.MI_DL.*fdcomm_save{ii}.alpha_DL;
%     I_radar_weighted = radar_save{ii}.alpha_r.*radar_save{ii}.MI_radar;
%     I_total_ULcomm(ii) = sum(I_UL_weighted(:));
%     I_total_DLcomm(ii) = sum(I_DL_weighted(:));
%     I_total_fdcomm(ii) = I_total_ULcomm(ii)+I_total_DLcomm(ii);
%     I_total_radar(ii) = sum(I_radar_weighted(:));
%     I_total(ii) = sum(I_DL_weighted(:))+I_total_ULcomm(ii)+sum(I_radar_weighted(:));
% end
plot(UL_SNR_sweep, I_total_avg,'m--o','LineWidth',1.5)
hold on
plot(UL_SNR_sweep, I_total_radar_avg,'r-+','LineWidth',1.5)
hold on
plot(UL_SNR_sweep, I_total_FD_avg,'g:*','LineWidth',1.5)
hold on
plot(UL_SNR_sweep, I_total_UL_avg,'k-.x','LineWidth',1.5)
hold on
plot(UL_SNR_sweep, I_total_DL_avg,'b:^','LineWidth',1.5)
hold off
grid on
axis square
xlabel('$\mathrm{SNR}_{\textrm{UL}}\textrm{ (dB)}$','FontSize',12,'Interpreter','latex');
%ylabel('$\textrm{Weighted MI}\textrm{ (bits per channel use)}$','Interpreter','latex','FontSize',12)
ylabel('Weighted MI (bits per channel use)', 'FontSize',12)
legend('total MI', 'radar MI', 'FD MI', 'UL MI','DL MI','FontSize',10)
legend('Location','best');
legend('boxoff');
print('-depsc2','-r150','UL_SNR_sweep.eps')