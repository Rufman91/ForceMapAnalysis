h = waitbar(0,'Yikes');
NMaps = 21;
for k=1:Nmaps
    waitbar(k/NMaps,h,sprintf('%i/%i processed',k,Nmaps))
%     FM{k}.estimate_cp_rov(20);
%     FM{k}.CP_CNN_predict('Fast')
%     FM{k}.estimate_cp_cnn('Dropout',100);
%     FM{k}.old_CP
%     FM{k}.DropoutNet = DropoutNet;
%     FM{k}.estimate_cp_gof;
%     FM{k}.estimate_cp_combined;
%     [E_CNN{k},Hertz_CNN{k}] = FM{k}.calculate_e_mod('CNN','parabolic');
%     [E_old{k},Hertz_old{k}] = FM{k}.calculate_e_mod('OLD','parabolic');
%     [E_RoV{k},Hertz_RoV{k}] = FM{k}.calculate_e_mod('RoV','parabolic');
%     [E_GoF{k},Hertz_GoF{k}] = FM{k}.calculate_e_mod('GoF','parabolic');
%     [E_Combo{k},Hertz_Combo{k}] = FM{k}.calculate_e_mod('Combo','parabolic');
%     [E_Manual{k},Hertz_Manual{k}] = FM{k}.calculate_e_mod('manual','parabolic');
    Mean(k,1) = mean(E_CNN{k}(find(FM{k}.SelectedCurves)));
    Mean(k,2) = mean(E_old{k}(find(FM{k}.SelectedCurves)));
    Mean(k,3) = mean(E_RoV{k}(find(FM{k}.SelectedCurves)));
    Mean(k,4) = mean(E_GoF{k}(find(FM{k}.SelectedCurves)));
    Mean(k,5) = mean(E_Combo{k}(find(FM{k}.SelectedCurves)));
    Mean(k,6) = mean(E_Manual{k}(find(FM{k}.SelectedCurves)));
    STD(k,1) = std(E_CNN{k}(find(FM{k}.SelectedCurves)));
    STD(k,2) = std(E_old{k}(find(FM{k}.SelectedCurves)));
    STD(k,3) = std(E_RoV{k}(find(FM{k}.SelectedCurves)));
    STD(k,4) = std(E_GoF{k}(find(FM{k}.SelectedCurves)));
    STD(k,5) = std(E_Combo{k}(find(FM{k}.SelectedCurves)));
    STD(k,6) = std(E_Manual{k}(find(FM{k}.SelectedCurves)));
end
close(h);
f = figure('Name','Look at it!');
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
plot(FM{k}.THApp{i},FM{k}.BasedApp{i})
drawpoint('Position',[FM{k}.Man_CP(i,1) FM{k}.Man_CP(i,2)],'Color','green');
drawpoint('Position',[FM{k}.CP(i,1) FM{k}.CP(i,2)],'Color','blue');
drawpoint('Position',[FM{k}.CP_old(i,1) FM{k}.CP_old(i,2)],'Color','black');
drawpoint('Position',[FM{k}.CP_RoV(i,1) FM{k}.CP_RoV(i,2)],'Color','yellow');
drawpoint('Position',[FM{k}.CP_GoF(i,1) FM{k}.CP_GoF(i,2)],'Color','white');
drawpoint('Position',[FM{k}.CP_Combo(i,1) FM{k}.CP_Combo(i,2)],'Color','red');

scatter(find(FM{k}.SelectedCurves),E_CNN{k}(find(FM{k}.SelectedCurves)));
errorbar([1:6],Mean(k,:),STD(k,:))
xticklabels({'A','B','C','D'})
k = k + 1;