function CP_Seq_predict(obj,Network)
ImgSize = Network.Layers(1).InputSize;
objcell{1,1} = obj;
X = CP_Sequential_batchprep(objcell,ImgSize(1));
Ypredicted = predict(Network,X);
iRange = find(obj.selected_curves);
k = 1;
for i=iRange'
    obj.CP(i,1) = Ypredicted(k,1)*range(obj.thapp{i})+min(obj.thapp{i});
    obj.CP(i,2) = Ypredicted(k,2)*range(obj.basedapp{i})+min(obj.basedapp{i});
    k = k + 1;
end
