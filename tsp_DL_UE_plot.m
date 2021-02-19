set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
load('tsp_DL_UE_new.mat')
I_total = zeros(length(JJ),1);
I_total_fdcomm = zeros(length(JJ),1);
I_total_radar = zeros(length(JJ),1);
I_total_ULcomm = zeros(length(JJ),1);
I_total_DLcomm = zeros(length(JJ),1);
% 
% for JJ = 1:length(JJ)
%     I_total(JJ) = fdcomm_save{JJ}.I_total_op(end)+sum(radar_save{JJ}.MI_radar);
%     I_total_fdcomm(JJ) = fdcomm_save{JJ}.I_total_op(end);
%     I_total_radar(JJ) = sum(radar_save{JJ}.MI_radar);
%     I_total_ULcomm(JJ) = sum(fdcomm_save{JJ}.MI_UL(:));
%     I_total_DLcomm(JJ) = sum(fdcomm_save{JJ}.MI_DL(:));
% end
for jj = 1:length(JJ)
    I_total(jj) = fdcomm_save{jj}.I_total(end);
    I_total_radar(jj) = fdcomm_save{jj}.I_radar_op(end);
    I_total_ULcomm(jj) = fdcomm_save{jj}.I_UL_op(end);
    I_total_DLcomm(jj) = fdcomm_save{jj}.I_DL_op(end);
    I_total_fdcomm(jj) = I_total_ULcomm(jj) + I_total_DLcomm(jj);
end
figure();
plot(JJ, I_total,'m--o','LineWidth',1.5)
hold on
plot(JJ, I_total_radar,'r-+','LineWidth',1.5)
hold on
plot(JJ, I_total_fdcomm,'b:^','LineWidth',1.5)
hold on
plot(JJ, I_total_ULcomm,'k-.x','LineWidth',1.5)
hold on
plot(JJ, I_total_DLcomm,'g:*','LineWidth',1.5)
hold off
grid on
axis square
%xlabel('$\mathrm{CNR}_{\textrm{r}}\textrm{ (dB)}$','FontSize',10,'Interpreter','latex');
xlabel('Number of DL UEs', 'FontSize',12)
%ylabel('$\mathrm{MI}\textrm{ (bits per channel use)}$','Interpreter','latex','FontSize',12)
ylabel('Weighted MI (bits per channel use)', 'FontSize',12)
legend('total MI','radar MI', 'FD MI','UL MI', 'DL MI','FontSize',10)
legend('Location','east');
legend('boxoff');
print('-depsc2','-r150','tsp_DL_UE.eps')