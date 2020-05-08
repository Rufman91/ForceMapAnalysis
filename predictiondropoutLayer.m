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
            if ~isa(X, 'dlarray')
                superfloatOfX = superiorfloat(X);
            else
                superfloatOfX = superiorfloat(extractdata(X));
            end
            dropoutScaleFactor = cast( 1 - layer.p, superfloatOfX );
            dropoutMask = ( rand(size(X), 'like', X) > layer.p ) / dropoutScaleFactor;
            Z = X.*dropoutMask;
        end
    end
end