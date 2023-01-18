function [MaxIdx, MinIdx, MaxValue, NPeaks] = PC2HM_findPeakVectors(vectors, varargin)
    p = inputParser;
    defaultThreshold = 0;
    defaultMinPeakDistance = 1;
    defaultMinPeakHeight = -inf;
    addParameter(p,'Threshold',defaultThreshold);
    addParameter(p,'MinPeakDistance',defaultMinPeakDistance);
    addParameter(p,'MinPeakHeight',defaultMinPeakHeight);
    parse(p,varargin{:});
    
    for i=1:size(vectors,2)
        [peaks, locs] = findpeaks(vectors(:,i), 'Threshold', p.Results.Threshold, 'MinPeakDistance', p.Results.MinPeakDistance, 'MinPeakHeight', p.Results.MinPeakHeight);
        if isempty(peaks)
            [peaks, locs] = max(vectors(:,i));
        end
        [NPeaks(i), ~] = size(peaks);
        [~, maxIdx] = max(locs);
        [~, minIdx] = min(locs);
        [~, maxValueIdx] = max(peaks);
        MaxIdx(i) = locs(maxIdx(1));
        MinIdx(i) = locs(minIdx(1));
        MaxValue(i) = locs(maxValueIdx(1));
    end
end