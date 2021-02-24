set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
load('tsp_DL_UE_avg.mat')
I_total = zeros(length(JJ),M);
I_total_fdcomm = zeros(length(JJ),M);
I_total_radar = zeros(length(JJ),M);
I_total_ULcomm = zeros(length(JJ),M);
I_total_DLcomm = zeros(length(JJ),M);
% 
% for JJ = 1:length(JJ)
%     I_total(JJ) = fdcomm_save{JJ}.I_total_op(end)+sum(radar_save{JJ}.MI_radar);
%     I_total_fdcomm(JJ) = fdcomm_save{JJ}.I_total_op(end);
%     I_total_radar(JJ) = sum(radar_save{JJ}.MI_radar);
%     I_total_ULcomm(JJ) = sum(fdcomm_save{JJ}.MI_UL(:));
%     I_total_DLcomm(JJ) = sum(fdcomm_save{JJ}.MI_DL(:));
% end
for mm = 1:M
    for jj = 1:length(JJ)
        I_total(jj,mm) = fdcomm_save{jj,mm}.I_total(end);
        I_total_radar(jj,mm) = fdcomm_save{jj,mm}.I_radar_op(end);
        I_total_ULcomm(jj,mm) = fdcomm_save{jj,mm}.I_UL_op(end);
        I_total_DLcomm(jj,mm) = fdcomm_save{jj,mm}.I_DL_op(end);
        I_total_fdcomm(jj,mm) = I_total_ULcomm(jj,mm) + I_total_DLcomm(jj,mm);
    end
end
I_total_avg = mean(I_total,2);
I_total_radar_avg = mean(I_total_radar,2);
I_total_fdcomm_avg = mean(I_total_fdcomm,2);
I_total_ULcomm_avg = mean(I_total_ULcomm,2);
I_total_DLcomm_avg = mean(I_total_DLcomm,2);
figure();
plot(JJ, I_total_avg,'m--o','LineWidth',1.5)
hold on
plot(JJ, I_total_radar_avg,'r-+','LineWidth',1.5)
hold on
plot(JJ, I_total_fdcomm_avg,'b:^','LineWidth',1.5)
hold on
plot(JJ, I_total_ULcomm_avg,'k-.x','LineWidth',1.5)
hold on
plot(JJ, I_total_DLcomm_avg,'g:*','LineWidth',1.5)
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