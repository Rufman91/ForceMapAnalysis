classdef predictiondropoutLayer < nnet.layer.Layer
    % This custom layer was created for estimating prediction uncertainty
    % during inference. 
    % Don't use this layer during training as it doesn't save the
    % parameters needed for backpropagation. As is, this layer will also
    % backpropagate through channels that have been set to zero, which is
    % not the usual procedure. This will lead to untested behaviour.
    
    properties
        Name % Layer Name
        p   % Dropout probability
    end
    
    methods
        function layer = predictiondropoutLayer(p,name)
            layer.p = p;
            layer.Name = name;
        end
        
        function Z = predict(layer,X)
            Dim = size(X);
            X = reshape(X,[],1);
            Z = zeros(size(X));
            for i=1:length(X)
                if layer.p < rand
                    Z(i) = X(i)*1/(1-layer.p);
                else
                    Z(i) = 0;
                end
            end
            Z = reshape(Z,Dim);
        end
    end
end