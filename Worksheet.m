% h = waitbar(0,'setting up');
% for i=1:40
%     waitbar(i/40,h,sprintf('Fibril %d/40',i));
% %     E.FM{i}.create_fibril_mask;
% %     E.FM{i}.DropoutNet = MonteCarloDrop;
% %     E.FM{i}.level_height_map;
% %     E.FM{i}.show_e_mod_map;
% %     E.FM{i}.calc_fib_diam;
% %     E.FM{i}.estimate_cp_cnn('Dropout',100);
%     E.FM{i}.calculate_e_mod_hertz('CNN','parabolic');
% %     E.SPM{i}.create_fibril_mask(0.5);
% %     E.SPM{i}.calc_fib_diam_dband;
% %     E.SPM{i}.calc_glass_fib_potdiff;
% %     E.SPM{i}.manual_exclusion
% %     E.SPM{i}.calc_glass_fib_potdiff;
% %     subplot(3,1,2)
% %     imshow(E.SPM{i}.FPotMask+E.SPM{i}.GPotMask)
% %     subplot(3,1,1)
% %     imshow(E.SPM{i}.HeightBased,[min(E.SPM{i}.HeightBased,[],'all') max(E.SPM{i}.HeightBased,[],'all')])
% %     title(E.FM{i}.Name)
% %     subplot(3,1,3)
% %     imshow(E.SPM{i}.Potential,[min(E.SPM{i}.Potential,[],'all') max(E.SPM{i}.Potential,[],'all')])
% %     pause(1.5)
% %     FibPot(i) = E.SPM{i}.FibPot;
% end
% close(h)
% E.FM{1}.show_e_mod_map
% % k = randi(40);
% j = randi(1024);
% plot(E.FM{k}.THApp{j}, E.FM{k}.BasedApp{j})
% drawpoint('Position',[mean(E.FM{k}.CP_MonteCarlo(:,1,j)) mean(E.FM{k}.CP_MonteCarlo(:,2,j))],'Color','green');
% drawpoint('Position',[mean(E.FM{k}.CP(j,1)) mean(E.FM{k}.CP(j,2))],'Color','red');
% title(sprintf('L1 unceretainty MonteCarlo14:\n X-Axis: %e\n Y-Axis: %e',std((E.FM{k}.CP_MonteCarlo(:,1,j))),std(E.FM{k}.CP_MonteCarlo(:,2,j))));
% 
% 
% for i=1:1024
%     F = fit(E.FM{k}.HHApp{i}(end-4:end),E.FM{k}.App{i}(end-4:end)/0.48,'poly1');
%     DzSlope(i) = F.p1;
% end
% 
% [OrderedDz,Index] = sort(DzSlope,'descend');

% plot(E.FM{1}.HHApp{Index(k)},E.FM{1}.App{Index(k)}./E.FM{1}.Header{17,2})
% title(sprintf('Slope = %e',OrderedDz(k)))
% k = k + 1;
% 
% histogram(OrderedDz)

% for i=1:40
%     RoV(i,:) = E.FM{i}.EModOliverPharr_RoV(E.FM{i}.RectApexIndex);
%     for j=1:length(E.FM{i}.RectApexIndex)
%         if RoV(i,j) > (nanmedian(RoV(i,:))+2.5*iqr(RoV(i,:))) || ...
%                 RoV(i,j) < (nanmedian(RoV(i,:))-2.5*iqr(RoV(i,:))) || ...
%                 E.FM{i}.ExclMask(E.FM{i}.List2Map(E.FM{i}.RectApexIndex(j),1),E.FM{i}.List2Map(E.FM{i}.RectApexIndex(j),2)) == 0
%             RoV(i,j) = NaN;
%         end
%     end
% end
% for i=1:40
%     Old(i,:) = E.FM{i}.EModOliverPharr_Old(E.FM{i}.RectApexIndex);
%     for j=1:length(E.FM{i}.RectApexIndex)
%         if Old(i,j) > (nanmedian(Old(i,:))+2.5*iqr(Old(i,:))) || ...
%                 Old(i,j) < (nanmedian(Old(i,:))-2.5*iqr(Old(i,:))) || ...
%                 E.FM{i}.ExclMask(E.FM{i}.List2Map(E.FM{i}.RectApexIndex(j),1),E.FM{i}.List2Map(E.FM{i}.RectApexIndex(j),2)) == 0
%             Old(i,j) = NaN;
%         end
%     end
% end
% for i=1:40
%     CNN(i,:) = E.FM{i}.EModOliverPharr_CNN(E.FM{i}.RectApexIndex);
%     for j=1:length(E.FM{i}.RectApexIndex)
%         if CNN(i,j) > (nanmedian(CNN(i,:))+2.5*iqr(CNN(i,:))) || ...
%                 CNN(i,j) < (nanmedian(CNN(i,:))-2.5*iqr(CNN(i,:))) || ...
%                 E.FM{i}.ExclMask(E.FM{i}.List2Map(E.FM{i}.RectApexIndex(j),1),E.FM{i}.List2Map(E.FM{i}.RectApexIndex(j),2)) == 0
%             CNN(i,j) = NaN;
%         end
%     end
% end
% k=randi(40)
% RandApex = E.FM{k}.RectApexIndex(randperm(length(E.FM{k}.RectApexIndex)));
% i=RandApex(1)
% figure('Color','w');
% hold on
% plot((E.FM{k}.HHApp{i}-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9,E.FM{k}.BasedApp{i}*10e9,(E.FM{k}.HHRet{i}-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9,E.FM{k}.BasedRet{i}*10e9);
% plot((E.FM{k}.CP_OliverPharr_CNN(i,1)-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9, E.FM{k}.CP_OliverPharr_CNN(i,2)*10e9,'O',...
%     'LineWidth',1.5,...
%     'MarkerSize',7,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','g');
% plot((E.FM{k}.CP_OliverPharr_Old(i,1)-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9, E.FM{k}.CP_OliverPharr_Old(i,2)*10e9,'O',...
%     'LineWidth',1.5,...
%     'MarkerSize',7,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','y');
% plot((E.FM{k}.CP_OliverPharr_RoV(i,1)-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9, E.FM{k}.CP_OliverPharr_RoV(i,2)*10e9,'O',...
%     'LineWidth',1.5,...
%     'MarkerSize',7,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','r');
% xlabel('Z-Displacement [nm]');
% ylabel('Deflection [nm]');
% grid on
% grid minor
% dim = [0.4 0.5 0.3 0.3];
% str = {sprintf('E_{CNN} = %.3f MPa',E.FM{k}.EModOliverPharr_CNN(i)*1e-6),...
%     sprintf('E_{SD6} = %.3f MPa',E.FM{k}.EModOliverPharr_Old(i)*1e-6),...
%     sprintf('E_{RoV} = %.3f MPa',E.FM{k}.EModOliverPharr_RoV(i)*1e-6)};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% legend('Approach','Retract','CNN','SD6','RoV','Location','northwest')
% 
% 
% k=randi(40)
% RandApex = E.FM{k}.RectApexIndex(randperm(length(E.FM{k}.RectApexIndex)));
% i=RandApex(1)
% [~,Idx] = max(E.FM{k}.CP_OliverPharr_MonteCarlo_STD);
% figure('Color','w');
% hold on
% b = plot((E.FM{k}.HHApp{i}-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9,E.FM{k}.BasedApp{i}*10e9,(E.FM{k}.HHRet{i}-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9,E.FM{k}.BasedRet{i}*10e9,'LineWidth',1.5);
% a = plot((E.FM{k}.CP_OliverPharr_CNN(i,1)-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9, E.FM{k}.CP_OliverPharr_CNN(i,2)*10e9,'O',...
%     'LineWidth',1.5,...
%     'MarkerSize',7,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','g');
% for j=1:100
%     plot((E.FM{k}.CP_OliverPharr_MonteCarlo(j,1,i)-E.FM{k}.CP_OliverPharr_CNN(i,1))*10e9, E.FM{k}.CP_OliverPharr_MonteCarlo(j,2,i)*10e9,'O',...
%         'LineWidth',0.5,...
%         'MarkerSize',7,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor','y');
% end
% legend('Approach','Retract','Mean Prediction','Dropout Predictions','Location','northwest');
% uistack(b,'top');
% uistack(a,'top');
% xlabel('Z-Displacement [nm]');
% ylabel('Deflection [nm]');
% grid on
% grid minor
% dim = [0.4 0.3 0.3 0.3];
% str = {'     Network Uncertainty',sprintf('Standard Deviation = %.2f nm',E.FM{k}.CP_OliverPharr_MonteCarlo_STD(i)*1e9)};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% i = i + 1;
% 
% Act = activations(E.CP_CNN,X(:,:,1,randi(length(X))),'conv_15');
% figure()
% I = imtile(mat2gray(Act),'GridSize',[8 8]);
% imshow(I)

% % Dry-Wet plots
% subplot(2,2,1)
% plot(1,Dry(1:10),'rO',2,Dry(11:20),'rO',[1 2],[mean(Dry(1:10)) mean(Dry(11:20))],'r-X')
% xticks([1 2])
% ylim([(min(Dry(1:20))-10) (max(Dry(1:20))+10)])
% xlim([0 3])
% xticklabels({'Control Before','Control After'})
% ylabel('Dry Diameter [nm]')
% subplot(2,2,2)
% plot(1,Dry(21:30),'rO',2,Dry(31:40),'rO',[1 2],[mean(Dry(21:30)) mean(Dry(31:40))],'r-X')
% xticks([1 2])
% ylim([(min(Dry(21:40))-10) (max(Dry(21:40))+10)])
% xlim([0 3])
% xticklabels({'MGO Before','MGO After'})
% ylabel('Dry Diameter [nm]')
% subplot(2,2,3)
% plot(1,Wet(1:10),'bO',2,Wet(11:20),'bO',[1 2],[mean(Wet(1:10)) mean(Wet(11:20))],'b-X')
% xticks([1 2])
% ylim([(min(Wet(1:20))-10) (max(Wet(1:20))+10)])
% xlim([0 3])
% xticklabels({'Control Before','Control After'})
% ylabel('Wet Diameter [nm]')
% subplot(2,2,4)
% plot(1,Wet(21:30),'bO',2,Wet(31:40),'bO',[1 2],[mean(Wet(21:30)) mean(Wet(31:40))],'b-X')
% xticks([1 2])
% ylim([(min(Wet(21:40))-10) (max(Wet(21:40))+10)])
% xlim([0 3])
% xticklabels({'MGO Before','MGO After'})
% ylabel('Wet Diameter [nm]')

% Fast Batchprep Workplace: bilinear interpolation of HH-Based-Points;
% 
% ImgSize = 128;
% ImgSizeFinal = 128;
% X = zeros(ImgSizeFinal,ImgSizeFinal,1,1024);
% 
% for i=1:1024
%     Image = zeros(ImgSize,ImgSize);
%     
%     Points(:,1) = E.FM{1}.HHApp{i};
%     Points(:,2) = E.FM{1}.BasedApp{i};
%     Points(:,1) = (Points(:,1)-min(Points(:,1)))/range(Points(:,1))*(ImgSize-1);
%     Points(:,2) = (Points(:,2)-min(Points(:,2)))/range(Points(:,2))*(ImgSize-1);
%     
%     L = length(Points);
%     for j=1:L
%         
%         if j==1
%             Position = [floor(Points(j,1))+1,floor(Points(j,2))+1];
%         end
%         if j<L
%             NextPosition = [floor(Points(j+1,1))+1,floor(Points(j+1,2))+1];
%         end
%         Image((ImgSize+1)-Position(2),Position(1)) = 1;
%         
%         % fill out points between actual data points
%         if j<L
%             L1Norm = norm(Points(j+1,:)-Points(j,:));
%             IncrX = (Points(j+1,1) - Points(j,1))/L1Norm;
%             IncrY = (Points(j+1,2) - Points(j,2))/L1Norm;
%             FillerPos = Points(j,:);
%             while norm((FillerPos + [IncrX IncrY]) - Points(j,:)) < norm(Points(j+1,:) - Points(j,:))
%                 FillerPos(1) = FillerPos(1) + IncrX;
%                 FillerPos(2) = FillerPos(2) + IncrY;
%                 Image((ImgSize+1)-(floor(FillerPos(2))+1),floor(FillerPos(1))+1) = 1;
%             end
%         end
%         Position = NextPosition;
%     end
%     
%     Image = imresize(Image,[ImgSizeFinal ImgSizeFinal]);
%     X(:,:,1,i) = Image;
%     clear Points
% end

X = E.FM{1}.CP_batchprep_new_fast(E.FM,128);