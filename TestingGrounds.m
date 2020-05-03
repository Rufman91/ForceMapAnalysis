
Range = find(FM{6}.selected_curves);
plot(FM{6}.thapp{Range(i)},FM{6}.basedapp{Range(i)})
drawpoint('Position',[FM{6}.CP_old(Range(i),1) FM{6}.CP_old(Range(i),2)])
i = i + 1;