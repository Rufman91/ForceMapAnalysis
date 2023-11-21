function CantileverDimensioningTool()
%CREATEAFMGUI creates a visually appealing GUI for AFM cantilever optimization

% Define figure and set layout
fig = figure('Name','AFM Cantilever Optimization','NumberTitle','off','Position',[100 100 950 500]);
set(fig,'MenuBar','none','ToolBar','none');

stiffness = nan;
depth = nan;
radius = nan;
force = nan;
modulus = nan;
ind2bendratio = nan;
sensitivity = nan;
volt = nan;
AllowedNumError = 1e-3;

ErrorColor = [0.980392156862745,0.431372549019608,0.431372549019608];

BackPanelDetector = uicontrol('Style','text','String','','Position',[50 365 300 100]);
BackPanelCanti = uicontrol('Style','text','String','','Position',[50 265 300 100]);
DefaultColor = BackPanelCanti.BackgroundColor;
BackPanelMaterial = uicontrol('Style','text','String','','Position',[50 65 300 200]);

h = axes('Visible','off','Units','pixels','Position',[15 0 50 500]);
MaterialSection = text('String','Material','Units','pixels','Position',[10 165]);
set(MaterialSection, 'Rotation', 90);
CantiSection = text('String','Cantilever','Units','pixels','Position',[10 315]);
set(CantiSection, 'Rotation', 90);
DetectorSection = text('String','Detector','Units','pixels','Position',[10 415]);
set(DetectorSection, 'Rotation', 90);

% Define input fields and labels
voltLabel = uicontrol('Style','text','String','Sensor Voltage Difference (V)','Position',[50 450 150 30]);
voltInput = uicontrol('Style','edit','String','','Position',[200 450 150 30],'Callback',@valueChangedCallback);

sensitivityLabel = uicontrol('Style','text','String','Expected Optical Lever Sensitivity (m/V)','Position',[50 400 150 30]);
sensitivityInput = uicontrol('Style','edit','String','','Position',[200 400 150 30],'Callback',@valueChangedCallback);

stiffnessLabel = uicontrol('Style','text','String','Cantilever Stiffness k (N/m)','Position',[50 350 150 30]);
stiffnessInput = uicontrol('Style','edit','String','','Position',[200 350 150 30],'Callback',@valueChangedCallback);

ind2bendratioLabel = uicontrol('Style','text','String','Indent-to-Bend Ratio D/B (-)','Position',[50 300 150 30]);
ind2bendratioInput = uicontrol('Style','edit','String',num2str(1),'Position',[200 300 150 30],'Callback',@valueChangedCallback);

radiusLabel = uicontrol('Style','text','String','Tip Radius R (m)','Position',[50 250 150 30]);
radiusInput = uicontrol('Style','edit','String','','Position',[200 250 150 30],'Callback',@valueChangedCallback);

depthLabel = uicontrol('Style','text','String','Maximum Indentation Depth D (m)','Position',[50 200 150 30]);
depthInput = uicontrol('Style','edit','String','','Position',[200 200 150 30],'Callback',@valueChangedCallback);


forceLabel = uicontrol('Style','text','String','Maximum Force F (N)','Position',[50 150 150 30]);
forceInput = uicontrol('Style','edit','String','','Position',[200 150 150 30],'Callback',@valueChangedCallback);

modulusLabel = uicontrol('Style','text','String','Expected Indentation Modulus E (Pa)','Position',[50 100 150 30]);
modulusInput = uicontrol('Style','edit','String','','Position',[200 100 150 30],'Callback',@valueChangedCallback);

ErrorText1 = {'Error: Input values do not satisfy system of equations';...
    'Please erase one field that will then be recalculated to fit the model!'};
ErrorText2 = {'Not enough assumptions to fill in the model!';'Fill out 7 of the 8 fields!'};

errorLabel = uicontrol('Style','text','String',...
    ErrorText1,...
    'Position',[50 50 300 30],'ForegroundColor','red');
errorLabel.Visible = 0;

ax = axes('Position',[0.4 0.2 .55 .6]);
imshow('CantileverDimensioning.png', 'Parent', ax);

% Define callback functions
    function valueChangedCallback(~,~)
        % Update stiffness value and run determine_missing_values function
        volt = str2double(get(voltInput,'String'));
        if ~isnumeric(volt)
            volt = nan;
        end
        sensitivity = str2double(get(sensitivityInput,'String'));
        if ~isnumeric(sensitivity)
            sensitivity = nan;
        end
        stiffness = str2double(get(stiffnessInput,'String'));
        if ~isnumeric(stiffness)
            stiffness = nan;
        end
        depth = str2double(get(depthInput,'String'));
        if ~isnumeric(depth)
            depth = nan;
        end
        radius = str2double(get(radiusInput,'String'));
        if ~isnumeric(radius)
            radius = nan;
        end
        force = str2double(get(forceInput,'String'));
        if ~isnumeric(force)
            force = nan;
        end
        modulus = str2double(get(modulusInput,'String'));
        if ~isnumeric(modulus)
            modulus = nan;
        end
        ind2bendratio = str2double(get(ind2bendratioInput,'String'));
        if ~isnumeric(ind2bendratio)
            ind2bendratio = nan;
        end
        determine_missing_values(volt,sensitivity,stiffness,depth,radius,force,modulus,ind2bendratio);
    end


    function determine_missing_values(volt,sensitivity,stiffness,depth,radius,force,modulus,ind2bendratio)
        %DETERMINEMISSINGVALUES determines missing values for AFM cantilever optimization
        
        % Check if any of the input values are missing
        if isnan(stiffness) || isnan(depth) || isnan(radius) ||...
                isnan(force) || isnan(modulus) || isnan(ind2bendratio) ||...
                isnan(volt) || isnan(sensitivity)
            % Check which values are missing
            if isnan(volt)
                % Calculate volt using the given equation
                volt = depth/(ind2bendratio*sensitivity);
                set(voltInput,'String',num2str(volt));
            end
            if isnan(sensitivity)
                % Calculate sensitivity using the given equation
                sensitivity = depth/(ind2bendratio*volt);
                set(sensitivityInput,'String',num2str(sensitivity));
            end
            if isnan(stiffness)
                % Calculate stiffness using the given equation
                stiffness = ind2bendratio * force / depth;
                set(stiffnessInput,'String',num2str(stiffness));
            end
            if isnan(ind2bendratio)
                % Calculate ind2bendratio using the given equation
                ind2bendratio = depth * stiffness / force;
                set(ind2bendratioInput,'String',num2str(ind2bendratio));
            end
            if isnan(depth)
                % Calculate depth using the given equation
                depth = ((3*force*(1-0.5^2))/(4*modulus*sqrt(radius)))^(2/3);
                set(depthInput,'String',num2str(depth));
            end
            if isnan(radius)
                % Calculate radius using the given equation
                radius = ((3*force*(1-0.5^2))/(4*modulus*depth^(3/2)))^2;
                set(radiusInput,'String',num2str(radius));
            end
            if isnan(force)
                % Calculate force using the given equation
                force = modulus/(1-0.5^2)*(4*sqrt(radius)*depth^(3/2)/3);
                set(forceInput,'String',num2str(force));
            end
            if isnan(modulus)
                % Calculate modulus using the given equation
                modulus = 3*force*(1-0.5^2)/(4*sqrt(radius)*depth^(3/2));
                set(modulusInput,'String',num2str(modulus));
            end
            AnyMissing = 0;
            if  isnan(depth) || isnan(radius) || isnan(force) || isnan(modulus)
                errorLabel.String = ErrorText2;
                errorLabel.Visible = 1;
                AnyMissing = 1;
                set_material_error_background(true);
            else
                set_material_error_background(false);
            end
            if isnan(stiffness) || isnan(ind2bendratio)
                errorLabel.String = ErrorText2;
                errorLabel.Visible = 1;
                AnyMissing = 1;
                set_canti_error_background(true);
            else
                set_canti_error_background(false);
            end
            if isnan(volt) || isnan(sensitivity)
                errorLabel.String = ErrorText2;
                errorLabel.Visible = 1;
                AnyMissing = 1;
                set_detector_error_background(true);
            else
                set_detector_error_background(false);
            end
            if ~AnyMissing
                errorLabel.Visible = 0;
                set_material_error_background(false);
                set_canti_error_background(false);
                set_detector_error_background(false);
            end
            % If all input values are provided, check if the system of equations is solvable
            if abs(force - modulus/(1-0.5^2)*(4*sqrt(radius)*depth^(3/2)/3)) > force*AllowedNumError
                % If system is not solvable, display error message
                set_material_error_background(true);
                errorLabel.String = ErrorText1;
                errorLabel.Visible = 1;
            end
            if abs(ind2bendratio - depth * stiffness / force) > ind2bendratio*AllowedNumError
                set_canti_error_background(true);
                % If system is not solvable, display error message
                errorLabel.String = ErrorText1;
                errorLabel.Visible = 1;
            end
            if abs(volt - depth/(ind2bendratio*sensitivity)) > volt*AllowedNumError
                set_detector_error_background(true);
                % If system is not solvable, display error message
                errorLabel.String = ErrorText1;
                errorLabel.Visible = 1;
            end
        else
            % If all input values are provided, check if the system of equations is solvable
            if abs(force - modulus/(1-0.5^2)*(4*sqrt(radius)*depth^(3/2)/3)) > force*AllowedNumError
                % If system is not solvable, display error message
                set_material_error_background(true);
                errorLabel.String = ErrorText1;
                errorLabel.Visible = 1;
            end
            if abs(ind2bendratio - depth * stiffness / force) > ind2bendratio*AllowedNumError
                set_canti_error_background(true);
                % If system is not solvable, display error message
                errorLabel.String = ErrorText1;
                errorLabel.Visible = 1;
            end
            if abs(volt - depth/(ind2bendratio*sensitivity)) > volt*AllowedNumError
                set_detector_error_background(true);
                % If system is not solvable, display error message
                errorLabel.String = ErrorText1;
                errorLabel.Visible = 1;
            end
        end
    end

    function set_material_error_background(Switch)
        if Switch
            BackPanelMaterial.BackgroundColor = ErrorColor;
            depthLabel.BackgroundColor = ErrorColor;
            modulusLabel.BackgroundColor = ErrorColor;
            radiusLabel.BackgroundColor = ErrorColor;
            forceLabel.BackgroundColor = ErrorColor;
        else
            BackPanelMaterial.BackgroundColor = DefaultColor;
            depthLabel.BackgroundColor = DefaultColor;
            modulusLabel.BackgroundColor = DefaultColor;
            radiusLabel.BackgroundColor = DefaultColor;
            forceLabel.BackgroundColor = DefaultColor;
        end
    end

    function set_canti_error_background(Switch)
        if Switch
            BackPanelCanti.BackgroundColor = ErrorColor;
            stiffnessLabel.BackgroundColor = ErrorColor;
            ind2bendratioLabel.BackgroundColor = ErrorColor;
        else
            BackPanelCanti.BackgroundColor = DefaultColor;
            stiffnessLabel.BackgroundColor = DefaultColor;
            ind2bendratioLabel.BackgroundColor = DefaultColor;
        end
    end

    function set_detector_error_background(Switch)
        if Switch
            BackPanelDetector.BackgroundColor = ErrorColor;
            voltLabel.BackgroundColor = ErrorColor;
            sensitivityLabel.BackgroundColor = ErrorColor;
        else
            BackPanelDetector.BackgroundColor = DefaultColor;
            voltLabel.BackgroundColor = DefaultColor;
            sensitivityLabel.BackgroundColor = DefaultColor;
        end
    end
end

