import numpy as np

def deconvolute_by_mathematical_morphology(in_image, eroding_geometry):
    if in_image.shape != eroding_geometry.shape:
        raise ValueError("The Image and the Eroding Geometry need to have the same dimensions")

    in_pixels_x, in_pixels_y = in_image.shape
    eg_pixels_x, eg_pixels_y = eroding_geometry.shape
    out_image = np.ones_like(in_image)

    peak_index_x, peak_index_y = np.unravel_index(np.argmax(eroding_geometry), eroding_geometry.shape)

    print("Processing deconvolution...")

    for i in range(in_pixels_x):
        for j in range(in_pixels_y):
            s_xmin = max(-peak_index_x, -i)
            s_xmax = min(eg_pixels_x - peak_index_x, in_pixels_x - i) - 1
            s_ymin = max(-j, -peak_index_y)
            s_ymax = min(eg_pixels_y - peak_index_y, in_pixels_y - j) - 1

            x_range = np.arange(i + 1 + s_xmin, i + 1 + s_xmax)
            y_range = np.arange(j + 1 + s_ymin, j + 1 + s_ymax)

            x_idx, y_idx = np.meshgrid(x_range, y_range, indexing='ij')
            ex_idx = x_idx + s_xmin + peak_index_x
            ey_idx = y_idx + s_ymin + peak_index_y

            temp_mat = in_image[x_idx, y_idx] - eroding_geometry[ex_idx, ey_idx]
            out_image[i, j] = np.min(temp_mat)

        if i % round(in_pixels_x / 10) == 0:
            print(f"{(i / in_pixels_x) * 100:.0f}% complete")

    print("Processing deconvolution... Done!")
    return out_image
