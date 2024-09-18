function DLVO_Hertz_GUI_Corrected
    % Create the figure
    hFig = figure('Position', [100, 100, 1200, 600], 'Name', 'DLVO and Hertzian Force Plotter');

    % Title for the GUI
    uicontrol('Style', 'text', 'Position', [500, 570, 200, 20], 'String', 'DLVO and Hertzian Force Model', 'FontSize', 12, 'FontWeight', 'bold');

    % Input fields for DLVO and Hertzian parameters
    addParameterInput('Tip Radius (R) [m]:', '20e-9', [20, 500]);
    addParameterInput('Debye Length (\lambda_D) [m]:', '3e-9', [20, 470]);
    addParameterInput('Surface Potential 1 (\sigma_1) [V]:', '0.1', [20, 440]);
    addParameterInput('Surface Potential 2 (\sigma_2) [V]:', '0.1', [20, 410]);
    addParameterInput('Hamaker Constant (H) [J]:', '1e-20', [20, 380]);
    addParameterInput('Effective Young''s Modulus (E^*) [Pa]:', '1e9', [20, 350]);

    % Input fields for plot limits
    addParameterInput('X-axis Min (nm):', '-5', [20, 320]);
    addParameterInput('X-axis Max (nm):', '50', [20, 290]);
    addParameterInput('Y-axis Min (N):', '-1e-9', [20, 260]);
    addParameterInput('Y-axis Max (N):', '1e-9', [20, 230]);

    % Button to calculate and plot
    uicontrol('Style', 'pushbutton', 'Position', [20, 190, 260, 30], 'String', 'Plot DLVO and Hertzian Forces', ...
              'Callback', @plotForces);

    % Axes for the plot
    hAxes = axes('Parent', hFig, 'Position', [0.4, 0.15, 0.55, 0.75]);
    xlabel(hAxes, 'Separation Distance D (nm)');
    ylabel(hAxes, 'Force F(D) (N)');
    title(hAxes, 'DLVO and Hertzian Forces');
    grid(hAxes, 'on');

    % Callback function to plot forces
    function plotForces(~, ~)
        % Retrieve input values
        R = str2double(get(findobj('Tag', 'Tip Radius (R) [m]:'), 'String'));
        lambda_D = str2double(get(findobj('Tag', 'Debye Length (\lambda_D) [m]:'), 'String'));
        sigma_1 = str2double(get(findobj('Tag', 'Surface Potential 1 (\sigma_1) [V]:'), 'String'));
        sigma_2 = str2double(get(findobj('Tag', 'Surface Potential 2 (\sigma_2) [V]:'), 'String'));
        H = str2double(get(findobj('Tag', 'Hamaker Constant (H) [J]:'), 'String'));
        E_star = str2double(get(findobj('Tag', 'Effective Young''s Modulus (E^*) [Pa]:'), 'String'));

        xMin = str2double(get(findobj('Tag', 'X-axis Min (nm):'), 'String'));
        xMax = str2double(get(findobj('Tag', 'X-axis Max (nm):'), 'String'));
        yMin = str2double(get(findobj('Tag', 'Y-axis Min (N):'), 'String'));
        yMax = str2double(get(findobj('Tag', 'Y-axis Max (N):'), 'String'));

        % Permittivity of water at room temperature (F/m)
        epsilon = 80 * 8.854e-12; 

        % Separation distance range
        D = linspace(xMin * 1e-9, xMax * 1e-9, 500);

        % Calculate DLVO forces (for D > 0)
        DLVO_indices = D > 0;
        D_DLVO = D(DLVO_indices);
        DLVO_part1 = -H * R ./ (6 * D_DLVO.^2);
        DLVO_part2 = (2 * pi * lambda_D * R / epsilon) * ...
                     ((sigma_1^2 + sigma_2^2) .* exp(-2 * D_DLVO / lambda_D) + ...
                      2 * sigma_1 * sigma_2 .* exp(-D_DLVO / lambda_D));
        F_DLVO = DLVO_part1 + DLVO_part2;

        % Hertzian contact force (for D <= 0, i.e., indentation)
        Hertz_indices = D <= 0;
        D_Hertz = abs(D(Hertz_indices));  % Indentation depth is positive
        F_Hertz = (4/3) * E_star * sqrt(R) .* (D_Hertz).^(3/2);

        % Plotting the DLVO and Hertzian forces
        cla(hAxes);
        hold(hAxes, 'on');
        plot(hAxes, D_DLVO * 1e9, F_DLVO, 'b-', 'LineWidth', 1.5, 'DisplayName', 'DLVO Force');
        plot(hAxes, D(Hertz_indices) * 1e9, F_Hertz, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Hertzian Contact Force');
        hold(hAxes, 'off');
        
        % Set plot limits
        xlim(hAxes, [xMin, xMax]);
        ylim(hAxes, [yMin, yMax]);

        % Add legend
        legend(hAxes, 'show');
    end

    % Function to add input fields for parameters
    function addParameterInput(label, defaultVal, position)
        uicontrol('Style', 'text', 'Position', [position(1), position(2), 200, 20], 'String', label);
        uicontrol('Style', 'edit', 'Position', [position(1) + 210, position(2), 100, 20], 'String', defaultVal, 'Tag', label);
    end
end
