function original_pc = PC2HM_backtransform_pointcloud(normalized_pc, scale, offset)
    original_pc = normalized_pc .* scale + offset;
end
