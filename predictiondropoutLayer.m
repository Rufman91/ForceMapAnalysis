classdef predictiondropoutLayer < nnet.layer.Layer
    % This custom layer was created for estimating prediction uncertainty
    % during inference. 
    % Don't use this layer during training as it doesn't save the
    % parameters needed for backpropagation. As is, this layer will also
    % backpropagate through channels that have been set to zero, which is
    % not the usual procedure. This will lead to untested behaviour.
    
    properties
        p   % Dropout probability
    end
    
    methods
        function layer = predictiondropoutLayer(p,Name)
            layer.p = p;
            layer.Name = Name;
        end
        
        function Z = predict(layer,X)
            Dim(1) = size(X,1);
            Dim(2) = size(X,2);
            Dim(3) = size(X,3);
            Dim(4) = size(X,4);
            Z = zeros(Dim);
            for i=1:Dim(1)
                for j=1:Dim(2)
                    for k=1:Dim(3)
                        for l=1:Dim(4)
                            if layer.p < rand
                                Z(i,j,k,l) = 1/(1-layer.p);
                            end
                        end
                    end
                end
            end
            Z = X.*Z;
        end
    end
end