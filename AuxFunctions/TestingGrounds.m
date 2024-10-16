h = waitbar(0,'Yikes');
NMaps = 21;
for k=1:Nmaps
    waitbar(k/NMaps,h,sprintf('%i/%i processed',k,Nmaps))
%     FM{k}.estimate_cp_rov(20);
    FM{k}.DropoutNet = MonteCarloDrop;
%     FM{k}.estimate_cp_cnn('Fast');
    FM{k}.estimate_cp_cnn('Dropout',100);
%     FM{k}.estimate_cp_old
%     FM{k}.DropoutNet = DropoutNet;
%     FM{k}.estimate_cp_gof;
%     FM{k}.estimate_cp_combined;
%     [E_CNN{k},Hertz_CNN{k}] = FM{k}.calculate_e_mod('CNN','parabolic');
%     [E_old{k},Hertz_old{k}] = FM{k}.calculate_e_mod('OLD','parabolic');
%     [E_RoV{k},Hertz_RoV{k}] = FM{k}.calculate_e_mod('RoV','parabolic');
%     [E_GoF{k},Hertz_GoF{k}] = FM{k}.calculate_e_mod('GoF','parabolic');
%     [E_Combo{k},Hertz_Combo{k}] = FM{k}.calculate_e_mod('Combo','parabolic');
%     [E_Manual{k},Hertz_Manual{k}] = FM{k}.calculate_e_mod('manual','parabolic');
%     Mean(k,1) = mean(E_CNN{k}(find(FM{k}.SelectedCurves)));
%     Mean(k,2) = mean(E_old{k}(find(FM{k}.SelectedCurves)));
%     Mean(k,3) = mean(E_RoV{k}(find(FM{k}.SelectedCurves)));
%     Mean(k,4) = mean(E_GoF{k}(find(FM{k}.SelectedCurves)));
%     Mean(k,5) = mean(E_Combo{k}(find(FM{k}.SelectedCurves)));
%     Mean(k,6) = mean(E_Manual{k}(find(FM{k}.SelectedCurves)));
%     STD(k,1) = std(E_CNN{k}(find(FM{k}.SelectedCurves)));
%     STD(k,2) = std(E_old{k}(find(FM{k}.SelectedCurves)));
%     STD(k,3) = std(E_RoV{k}(find(FM{k}.SelectedCurves)));
%     STD(k,4) = std(E_GoF{k}(find(FM{k}.SelectedCurves)));
%     STD(k,5) = std(E_Combo{k}(find(FM{k}.SelectedCurves)));
%     STD(k,6) = std(E_Manual{k}(find(FM{k}.SelectedCurves)));
end
close(h);
f = figure('Name','Look at it!','Color','w');
%GoF
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
subplot(2,1,1)
plot(FM{k}.BasedApp{i})
drawpoint('Position',[260 0],'Color','w')
title('Force')
subplot(2,1,2)
plot(FM{k}.RoV{i})
title('R-Squared - GoF')

f = figure('Name','Look at it!','Color','w');
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
Which = randi(length(Curves));
i = Curves(Which);
plot(FM{k}.THApp{i},FM{k}.BasedApp{i})
for j=1:100
    drawpoint('Position',[FM{k}.CP_MonteCarlo(j,1,i) FM{k}.CP_MonteCarlo(j,2,i)],'Color','yellow');
    if j == 100
        drawpoint('Position',[FM{k}.Man_CP(i,1) FM{k}.Man_CP(i,2)],'Color','green');
        drawpoint('Position',[mean(FM{k}.CP_MonteCarlo(:,1,i)) mean(FM{k}.CP_MonteCarlo(:,2,i))],'Color','red');
    end
end

% RoV
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
subplot(2,1,1)
plot(FM{k}.BasedApp{i})
drawpoint('Position',[250 0],'Color','y')
title('Force')
subplot(2,1,2)
plot(FM{k}.RoV{i})
title('RoV')

% Combined
f = figure('Name','Look at it!','Color','w');
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
subplot(2,1,1)
plot(FM{k}.BasedApp{i})
drawpoint('Position',[356 0],'Color','red')
title('Force')
subplot(2,1,2)
plot([1:length(FM{k}.BasedApp{i})],FM{k}.CPComboCurve{i},[1:length(FM{k}.BasedApp{i})],FM{k}.RoV{i},[1:length(FM{k}.BasedApp{i})],FM{k}.GoF{i})
legend('RoV*GoF','RoV','GoF')
title('RoV*GoF')

% Old
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
subplot(2,1,1)
plot(FM{k}.BasedApp{i})
drawpoint('Position',[250 0],'Color','y')
title('Force')
subplot(2,1,2)
plot(FM{k}.RoV{i})
title('RoV')

% Force curve with labels
f = figure('Name','Look at it!','Color','w');
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
plot(FM{k}.THApp{i},FM{k}.BasedApp{i},FM{k}.THRet{i},FM{k}.BasedRet{i})
xlabel('Cantilever Tip Height in [m]');
ylabel('Force in [N]');
legend('Approach','Retract')
title('Typical AFM Force Curves')
drawpoint('Position',[-1.3e-6 0],'Color','red','Label','Contact Point???')

%%%%PRESENTATION%%%%%%

f = figure('Name','Look at it!','Color','w');
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
plot(FM{k}.THApp{i},FM{k}.BasedApp{i})
title('green-Manual  blue-CNN  black-Der.STD  yellow-RoV  white-GoF  red-Combined')
drawpoint('Position',[FM{k}.Man_CP(i,1) FM{k}.Man_CP(i,2)],'Color','green');
drawpoint('Position',[FM{k}.CP(i,1) FM{k}.CP(i,2)],'Color','blue');
drawpoint('Position',[FM{k}.CP_old(i,1) FM{k}.CP_old(i,2)],'Color','black');
drawpoint('Position',[FM{k}.CP_RoV(i,1) FM{k}.CP_RoV(i,2)],'Color','yellow');
drawpoint('Position',[FM{k}.CP_GoF(i,1) FM{k}.CP_GoF(i,2)],'Color','white');
drawpoint('Position',[FM{k}.CP_Combo(i,1) FM{k}.CP_Combo(i,2)],'Color','red');

%%%%%%%%%%%%%%%%%%%%%%%%

errorbar([1:6],Mean(k,:),STD(k,:),'O')
axis([0 7 (min(Mean(k,:)-1.2*max(STD(k,:)))) (max(Mean(k,:)+1.2*max(STD(k,:))))])
ylabel('Mean E-Modulus [MPa]')
xticks([1:6])
title(FM{k}.Name)
xticklabels({'CNN','Derivative STD','RoV','GoF','Combo','Manual'})
k = k + 1;
m = 1;
for i=1:21
    while Val_idx(m) < FM{k}.NCurves
        ValCurves(m,1) = Val_idx(i);
    end
end
