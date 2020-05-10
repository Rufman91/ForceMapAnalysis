for k=1:21
%     FM{k}.CP_RoV
%     FM{k}.CP_CNN_predict('Fast')
%     FM{k}.CP_CNN_predict('Dropout',100)
%     FM{k}.old_CP
%     FM{k}.DropoutNet = DropoutNet;
    
end
k = randi(21);
Curves = find(FM{k}.SelectedCurves);
i = Curves(randi(length(Curves)));
f = figure('Name','Look at it!');
plot(FM{k}.THApp{i},FM{k}.BasedApp{i})
drawpoint('Position',[FM{k}.Man_CP(i,1) FM{k}.Man_CP(i,2)],'Color','green')
drawpoint('Position',[FM{k}.CP(i,1) FM{k}.CP(i,2)],'Color','blue')
drawpoint('Position',[FM{k}.CP_old(i,1) FM{k}.CP_old(i,2)],'Color','black')