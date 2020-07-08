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
E.FM{1}.show_e_mod_map
% k = randi(40);
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


% HEAD HEIGHT SCALING WTFFFFFFF?????