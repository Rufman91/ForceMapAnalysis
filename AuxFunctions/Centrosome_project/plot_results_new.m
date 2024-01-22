
for i = 1:E.NumForceMaps
    if ismember(i, datavalues)
        % skip evaluation round
    else
        cd(E.ForceMapFolders{i,1})
        load(strcat('Processed',s2,'.mat'))
        CsEModHertz_data{i} = CsEModHertz(:);
        CsEModHertz_mean(i) = mean(CsEModHertz(:), 'omitnan');
        CsHeight_mean(i) = mean(CsFlatHeight(:),'omitnan');
        CsInden_mean(i) = mean(CsFlatInden(:),'omitnan');
        CsRadius_data(i) = CsRadius;
        CsAspectRatio(i) = mean(CsFlatHeight(:),'omitnan')/(CsRadius*2);
        AngleThr(i) = T2;
        CsFlatArea_data(i) = CsFlatArea;
        CsVolume_data(i) = CsVolume;
    end
end

CsHeight_mean(CsHeight_mean == 0) = NaN; 
CsVolume_data(CsVolume_data == 0) = NaN; 
CsEModHertz_mean(CsEModHertz_mean == 0) = NaN; 
figure('name', 'Height dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(CsHeight_mean(i).*1e9, CsEModHertz_mean(i).*1e-3, 50, [0 0.4470 0.7410], "filled");
end
ylabel('Indentation modulus [kPa]'); 
xlabel('Centrosome height [nm]'); 

figure('name', 'Volume dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(CsVolume_data(i)*1e+18, CsEModHertz_mean(i).*1e-3, 50, [0 0.4470 0.7410], "filled");
end
ylabel('Indentation modulus [kPa]'); 
xlabel('Volume [\mum^3]');