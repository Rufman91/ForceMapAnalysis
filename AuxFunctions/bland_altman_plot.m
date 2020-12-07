function bland_altman_plot(X,Y)

Mean = (Y+X)/2;
Diff = Y-X;
figure('Color','w','Name','Bland-Altman_Plot','Units','normalized','Position',[0.2 0.6 0.3 0.3])

plot(Mean , Diff,'rO')
hold on
plot(Mean,mean(Diff)*ones(1,10),'b-')
plot(Mean, (mean(Diff)+1.96*std(Diff))*ones(1,10),'r-')
plot(Mean,(mean(Diff)-1.96*std(Diff))*ones(1,10),'r-')
title('Bland-Altman Absolute')
legend({'Diff vs. Mean','Mean of Diff','MoD +- 1.96*std(MoD)'})
hold off


Mean = sqrt(Y.*X);
Diff = log(Y./X);
figure('Color','w','Name','Bland-Altman_Plot','Units','normalized','Position',[0.5 0.6 0.3 0.3])

plot(Mean , Diff,'rO')
hold on
plot(Mean,mean(Diff)*ones(1,10),'b-')
plot(Mean, (mean(Diff)+1.96*std(Diff))*ones(1,10),'r-')
plot(Mean,(mean(Diff)-1.96*std(Diff))*ones(1,10),'r-')
title('Bland-Altman relative Log.')
legend({'LogDiff vs. GeoMean','Mean of LogDiff','MoLD +- 1.96*std(MoLD)'})
hold off

end