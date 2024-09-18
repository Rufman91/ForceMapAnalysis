function visualize_DLVO_force
% Create the figure
hFig = figure('Position', [100, 100, 1200, 600], 'Name', 'DLVO Force Calculator');

% Display DLVO formula using annotation (which supports LaTeX)
annotation('textbox', [0.5, 0.5, 0.5, 0.2], 'String', ...
    '$F_{DLVO}(D) = \frac{-HR}{6D^2} + \frac{2\pi \lambda_D R}{\epsilon} \left[ (\sigma_1^2 + \sigma_2^2) e^{-\frac{2D}{\lambda_D}} + 2\sigma_1\sigma_2 e^{-\frac{D}{\lambda_D}} \right]$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% Button to display DLVO theory explanation
uicontrol('Style', 'pushbutton', 'Position', [20, 500, 150, 30], 'String', 'What is DLVO?', ...
    'Callback', @(src, event) showDLVOExplanation());

% Create UI elements for input parameters with labels mapping to formula
uicontrol('Style', 'text', 'Position', [20, 430, 200, 20], 'String', 'Tip Radius (R) [m]:');
hRadius = uicontrol('Style', 'edit', 'Position', [230, 430, 100, 20], 'String', '20e-9', ...
    'TooltipString', 'Radius of the AFM tip (R)');

uicontrol('Style', 'text', 'Position', [20, 400, 200, 20], 'String', 'Debye Length (\lambda_D) [m]:');
hDebyeLength = uicontrol('Style', 'edit', 'Position', [230, 400, 100, 20], 'String', '3e-9', ...
    'TooltipString', 'Debye length (\lambda_D) relates to ionic strength of the solution');

uicontrol('Style', 'text', 'Position', [20, 370, 200, 20], 'String', 'Surface Potential 1 (\sigma_1) [V]:');
hSigma1 = uicontrol('Style', 'edit', 'Position', [230, 370, 100, 20], 'String', '0.1', ...
    'TooltipString', 'Surface potential of the AFM tip (\sigma_1)');

uicontrol('Style', 'text', 'Position', [20, 340, 200, 20], 'String', 'Surface Potential 2 (\sigma_2) [V]:');
hSigma2 = uicontrol('Style', 'edit', 'Position', [230, 340, 100, 20], 'String', '0.1', ...
    'TooltipString', 'Surface potential of the substrate (\sigma_2)');

uicontrol('Style', 'text', 'Position', [20, 310, 200, 20], 'String', 'Hamaker Constant (H) [J]:');
hHamaker = uicontrol('Style', 'edit', 'Position', [230, 310, 100, 20], 'String', '1e-20', ...
    'TooltipString', 'Hamaker constant (H) determines van der Waals attraction');

% Button to calculate and plot
uicontrol('Style', 'pushbutton', 'Position', [20, 270, 310, 30], 'String', 'Plot DLVO Force', ...
    'Callback', @(src, event) plotDLVO(hRadius, hDebyeLength, hSigma1, hSigma2, hHamaker));

% Axes for the plot
hAxes = axes('Parent', hFig, 'Position', [0.4, 0.1, 0.55, 0.8]);
xlabel(hAxes, 'Separation Distance D (nm)');
ylabel(hAxes, 'Force F(D) (N)');
title(hAxes, 'DLVO Force as a function of Separation Distance');
grid(hAxes, 'on');
end

function plotDLVO(hRadius, hDebyeLength, hSigma1, hSigma2, hHamaker)
% Retrieve user inputs
R = str2double(get(hRadius, 'String'));
lambda_D = str2double(get(hDebyeLength, 'String'));
sigma_1 = str2double(get(hSigma1, 'String'));
sigma_2 = str2double(get(hSigma2, 'String'));
H = str2double(get(hHamaker, 'String'));
epsilon_water = 80 * 8.854e-12;  % Permittivity of water at room temperature (F/m)

% Distance range (D)
D = linspace(1e-10, 5e-8, 2000);  % D in meters

% Calculate van der Waals force component
VdW_force = -H * R ./ (6 * D.^2);

% Calculate electrostatic force component
Electrostatic_force = (2 * pi * lambda_D * R / epsilon_water) * ...
    ((sigma_1^2 + sigma_2^2) .* exp(-2 * D / lambda_D) + ...
    2 * sigma_1 * sigma_2 .* exp(-D / lambda_D));

% Total DLVO force
F_DLVO = VdW_force + Electrostatic_force;

% Plotting the DLVO force and components
hFig = gcf;
hAxes = findobj(hFig, 'Type', 'axes');

% Clear previous plot
cla(hAxes);

hold(hAxes, 'on');
plot(hAxes, D * 1e9, F_DLVO * 1e9, 'b', 'LineWidth', 1.5);
plot(hAxes, D * 1e9, VdW_force * 1e9, '--r', 'LineWidth', 1.5);
plot(hAxes, D * 1e9, Electrostatic_force * 1e9, '-.g', 'LineWidth', 1.5);
hold(hAxes, 'off');

xlabel(hAxes, 'Separation Distance D (nm)');
ylabel(hAxes, 'Force F(D) (nN)');
title(hAxes, 'DLVO Force as a function of Separation Distance');
legend(hAxes, 'Total DLVO Force', 'van der Waals Force', 'Electrostatic Force');
grid(hAxes, 'on');
end

function showDLVOExplanation()
msgbox({'The DLVO theory describes the interaction between two surfaces in a liquid medium before they come into bulk contact.', ...
    'It combines the effects of:', ...
    '- van der Waals forces (attractive)', ...
    '- Electrostatic forces (repulsive or attractive depending on charge)', ...
    'The total force is a sum of these contributions, and the theory helps explain phenomena like colloidal stability.',...
    '',...
    'Cit.:',...
    'Butt, H.-J. (1991)',...
    'Electrostatic interaction in atomic force microscopy.',...
    'Biophys. J.,',...
    '60, 777–785.'}, ...
    'What is DLVO Theory?');
end
