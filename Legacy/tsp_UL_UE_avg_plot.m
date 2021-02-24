set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',10,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',8,...
'DefaultLineLineWidth',1,...
'DefaultLineMarkerSize',7.75);
load('tsp_UL_UE_avg.mat')
I_total = zeros(length(II),M);
I_total_fdcomm = zeros(length(II),M);
I_total_radar = zeros(length(II),M);
I_total_ULcomm = zeros(length(II),M);
I_total_DLcomm = zeros(length(II),M);
% 
% for ii = 1:length(II)
%     I_total(ii) = fdcomm_save{ii}.I_total_op(end)+sum(radar_save{ii}.MI_radar);
%     I_total_fdcomm(ii) = fdcomm_save{ii}.I_total_op(end);
%     I_total_radar(ii) = sum(radar_save{ii}.MI_radar);
%     I_total_ULcomm(ii) = sum(fdcomm_save{ii}.MI_UL(:));
%     I_total_DLcomm(ii) = sum(fdcomm_save{ii}.MI_DL(:));
% end
for m = 1:M
    for ii = 1:length(II)
        I_total(ii,m) = fdcomm_save{ii,m}.I_total(end);
        I_total_radar(ii,m) = fdcomm_save{ii,m}.I_radar_op(end);
        I_total_ULcomm(ii,m) = fdcomm_save{ii,m}.I_UL_op(end);
        I_total_DLcomm(ii,m) = fdcomm_save{ii,m}.I_DL_op(end);
        I_total_fdcomm(ii,m) = I_total_ULcomm(ii,m) + I_total_DLcomm(ii,m);
    end
end
I_total_avg = mean(I_total,2);
I_total_radar_avg = mean(I_total_radar,2);
I_total_fdcomm_avg = mean(I_total_fdcomm,2);
I_total_ULcomm_avg = mean(I_total_ULcomm,2);
I_total_DLcomm_avg = mean(I_total_DLcomm,2);
figure();
plot(II, I_total_avg,'m--o','LineWidth',1.5)
hold on
plot(II, I_total_radar_avg,'r-+','LineWidth',1.5)
hold on
plot(II, I_total_fdcomm_avg,'b:^','LineWidth',1.5)
hold on
plot(II, I_total_ULcomm_avg,'k-.x','LineWidth',1.5)
hold on
plot(II, I_total_DLcomm_avg,'g:*','LineWidth',1.5)
hold off
grid on
axis square
%xlabel('$\mathrm{CNR}_{\textrm{r}}\textrm{ (dB)}$','FontSize',10,'Interpreter','latex');
xlabel('Number of UL UEs', 'FontSize',12)
%ylabel('$\mathrm{MI}\textrm{ (bits per channel use)}$','Interpreter','latex','FontSize',12)
ylabel('Weighted MI (bits per channel use)', 'FontSize',12)
legend('total MI','radar MI', 'FD MI','UL MI', 'DL MI','FontSize',10)
legend('Location','east');
legend('boxoff');
print('-depsc2','-r150','II.eps')