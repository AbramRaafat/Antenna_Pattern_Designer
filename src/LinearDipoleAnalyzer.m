classdef LinearDipoleAnalyzer < matlab.apps.AppBase
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        MainGridLayout          matlab.ui.container.GridLayout
        ControlPanel            matlab.ui.container.Panel
        VisualTabGroup          matlab.ui.container.TabGroup
        PatternTab              matlab.ui.container.Tab
        CurrentTab              matlab.ui.container.Tab
        InfoTab                 matlab.ui.container.Tab
        PolarAxes               matlab.graphics.axis.PolarAxes 
        ThreeDAxes              matlab.ui.control.UIAxes
        CurrentAxes             matlab.ui.control.UIAxes
        LengthLabel             matlab.ui.control.Label
        LengthEditField         matlab.ui.control.NumericEditField
        ModelLabel              matlab.ui.control.Label
        ModelDropDown           matlab.ui.control.DropDown
        ShowHPBWCheckBox        matlab.ui.control.CheckBox
        CalculateButton         matlab.ui.control.Button
        ResultsPanel            matlab.ui.container.Panel
        DirectivityLabel        matlab.ui.control.Label
        RradLabel               matlab.ui.control.Label
        HPBWLabel               matlab.ui.control.Label
        WarningLabel            matlab.ui.control.Label
        FormulaTextArea         matlab.ui.control.TextArea
    end
    methods (Access = private)
        
        function UpdatePlots(app)
            disp('Starting Calculation...'); 
            
            L_lam = app.LengthEditField.Value; 
            model = app.ModelDropDown.Value;   
            
            % Physics Constants
            k = 2*pi;          % Wavenumber (normalized to lambda=1)
            H = L_lam / 2;     % Half-length
            
            z = linspace(-H, H, 500);
            
           % current distribution
            switch model
                case 'General (Sinusoidal)'
                    % standard standing wave assumption
                    % I(z) = sin(k(H - |z|)) note: this assume that current
                    % goes for zero at the end of the wire (Hence: using a infinitizm approx from general is needed)
                    Iz = sin(k * (H - abs(z)));
                    
                    % Normalization (Crucial for visualization of small L)
                    peak_I = max(abs(Iz));
                    if peak_I > 1e-9 
                        Iz = Iz / peak_I;
                    end
                    
                    % Dynamic Explanation
                    if L_lam < 0.1
                        model_desc = sprintf('Physics Note:\nAt L=%.2f (Small), sin(x) ~ x.\nThe sinusoidal wave approximates a TRIANGLE.\nIt is NOT uniform because current must be 0 at the ends.', L_lam);
                    else
                        model_desc = 'Physics Note:\nStandard sinusoidal standing wave distribution for finite wires.';
                    end
                    
                case 'Short (Triangular)'
                    % Approximation for L < lambda/10
                    % I(z) = I0 * (1 - |z|/H)
                    Iz = 1 - abs(z)/H;
                    model_desc = 'Approximation:\nLinear drop-off assumed for short wires.';
                    
                case 'Infinitesimal (Uniform)'
                    % I(z) = I0 (Constant)
                    Iz = ones(size(z));
                    model_desc = 'Idealization (Hertzian):\nAssumes capacitive plates at ends.\nImpossible for a simple wire.';
            end
            
            % Plot Current
            plot(app.CurrentAxes, z, Iz, 'LineWidth', 2, 'Color', 'b');
            yline(app.CurrentAxes, 0, 'k--');
            title(app.CurrentAxes, sprintf('Current Distribution (%s)', model));
            xlabel(app.CurrentAxes, 'z (\lambda)');
            ylabel(app.CurrentAxes, 'Magnitude (Normalized to Peak)');
            grid(app.CurrentAxes, 'on');
            ylim(app.CurrentAxes, [0 1.1]);

            hold(app.CurrentAxes, 'on');
            line(app.CurrentAxes, [-H H], [0 0], 'Color', 'k', 'LineWidth', 4);
            hold(app.CurrentAxes, 'off');
            
            % ====================== pattern =================
            theta = linspace(0, 2*pi, 1000); 
            
            switch model
                case 'General (Sinusoidal)'
                    % General Formula
                    num = cos(k * H * cos(theta)) - cos(k * H);
                    den = sin(theta);
                    F_theta = num ./ den;
                    
                    % Handle Singularities at 0, pi, 2pi
                    F_theta(abs(den) < 1e-10) = 0;
                    
                case {'Short (Triangular)', 'Infinitesimal (Uniform)'}
                    % Both Short and Infinitesimal have sin(theta) patterns
                    F_theta = sin(theta);
            end
            
            % Magnitude & Normalization
            E_mag = abs(F_theta);
            E_max = max(E_mag);
            if E_max < 1e-6; E_max = 1; end 
            E_norm = E_mag / E_max;

            % ==========================================================================
            
            % METRICS 
  
            theta_phys = linspace(eps, pi-eps, 500);
            
            if strcmp(model, 'General (Sinusoidal)')
                num_p = cos(k * H * cos(theta_phys)) - cos(k * H);
                den_p = sin(theta_phys);
                E_mag_phys = abs(num_p ./ den_p);
            else
                E_mag_phys = abs(sin(theta_phys));
            end
            E_norm_phys = E_mag_phys / max(E_mag_phys);
            
            % R_rad (Approximate formulas for display)
            integrand = (E_mag_phys.^2) .* sin(theta_phys);
            integral_val = trapz(theta_phys, integrand);
            
            if strcmp(model, 'General (Sinusoidal)')
                 if L_lam < 0.1
                     Rrad = 20 * pi^2 * (L_lam)^2; % Matches short dipole 
                 else
                     % This constant 60 factor is a simplification for the integral relation P = I^2 R
                     Rrad = 60 * integral_val; 
                 end
            elseif strcmp(model, 'Short (Triangular)')
                 Rrad = 20 * pi^2 * (L_lam)^2;
            else % Infinitesimal
                 Rrad = 80 * pi^2 * (L_lam)^2;
            end
            
            % Directivity
            integrand_norm = (E_norm_phys.^2) .* sin(theta_phys);
            beam_solid_angle = 2*pi * trapz(theta_phys, integrand_norm);
            D = 4*pi / beam_solid_angle;
            D_dBi = 10*log10(D);
            
            % HPBW
            [~, max_idx_p] = max(E_norm_phys);
            idx_left = max_idx_p;
            while idx_left > 1 && E_norm_phys(idx_left) > 0.707; idx_left=idx_left-1; end
            idx_right = max_idx_p;
            while idx_right < length(theta_phys) && E_norm_phys(idx_right) > 0.707; idx_right=idx_right+1; end
            
            width_rad = theta_phys(idx_right) - theta_phys(idx_left);
            HPBW_deg = rad2deg(width_rad);
            
            % Update Dashboard
            app.RradLabel.Text = sprintf('R_rad: %.2f \\Omega', Rrad);
            app.DirectivityLabel.Text = sprintf('Directivity: %.2f (%.2f dBi)', D, D_dBi);
            app.HPBWLabel.Text = sprintf('HPBW: %.1f^o', HPBW_deg);
            
            % Validation Warning
            if (strcmp(model, 'Short (Triangular)') || strcmp(model, 'Infinitesimal (Uniform)')) && L_lam > 0.1
                app.WarningLabel.Text = 'Warning: Model inaccurate for L > 0.1\lambda';
                app.WarningLabel.FontColor = [0.8 0 0]; % Red
            else
                app.WarningLabel.Text = sprintf('Valid %s Analysis', model);
                app.WarningLabel.FontColor = [0 0.5 0]; % Green
            end
            
            % ===== plots ====
            disp('Updating Visuals...'); 
            
            % 2D Polar
            cla(app.PolarAxes);
            polarplot(app.PolarAxes, theta, E_norm, 'LineWidth', 2, 'Color', 'r');
            app.PolarAxes.ThetaZeroLocation = 'top'; 
            app.PolarAxes.ThetaDir = 'clockwise';
            title(app.PolarAxes, sprintf('Pattern: %s', model));
            
            if app.ShowHPBWCheckBox.Value
                hold(app.PolarAxes, 'on');
                polarplot(app.PolarAxes, [theta_phys(idx_left) theta_phys(idx_left)], [0 1], '--b');
                polarplot(app.PolarAxes, [theta_phys(idx_right) theta_phys(idx_right)], [0 1], '--b');
                hold(app.PolarAxes, 'off');
            end
            
            % 3D Surface
            phi_3d = linspace(0, 2*pi, 120);
            [TH, PH] = meshgrid(theta_phys, phi_3d);
            R = repmat(E_norm_phys, 120, 1); 
            X = R .* sin(TH) .* cos(PH);
            Y = R .* sin(TH) .* sin(PH);
            Z = R .* cos(TH);
            
            cla(app.ThreeDAxes); 
            surf(app.ThreeDAxes, X, Y, Z, R, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
            colormap(app.ThreeDAxes, 'jet');
            shading(app.ThreeDAxes, 'interp');
            light(app.ThreeDAxes, 'Position', [1 1 1], 'Style', 'infinite');
            view(app.ThreeDAxes, 3);
            axis(app.ThreeDAxes, 'equal');
            title(app.ThreeDAxes, '3D Radiation Pattern');
            
            % Info tap
            app.FormulaTextArea.Value = {
                sprintf('Analysis Model: %s', model);
                sprintf('Length L = %.4f lambda', L_lam);
                '-----------------------------------';
                sprintf('Peak Angle: %.1f deg', rad2deg(theta_phys(max_idx_p)));
                sprintf('Rad Resistance: %.4f Ohms', Rrad);
                '';
                'Theory Notes:';
                model_desc;
            };
            
            disp('Update Complete.');
        end
        
        function CalculateButtonPushed(app, ~)
            disp('Button Click Detected');
            app.CalculateButton.Text = 'Calculating...';
            app.CalculateButton.Enable = 'off';
            drawnow;
            try
                app.UpdatePlots();
            catch ME
                disp(['Error: ', ME.message]);
                app.WarningLabel.Text = 'Error occurred.';
            end
            app.CalculateButton.Text = 'Calculate Pattern';
            app.CalculateButton.Enable = 'on';
            drawnow;
        end
    end
    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 950 600];
            app.UIFigure.Name = 'App 1: Linear Dipole Analyzer';
            app.MainGridLayout = uigridlayout(app.UIFigure);
            app.MainGridLayout.ColumnWidth = {280, '1x'}; 
            app.MainGridLayout.RowHeight = {'1x'};
            
            app.ControlPanel = uipanel(app.MainGridLayout);
            app.ControlPanel.Title = 'Configuration';
            app.ControlPanel.Layout.Column = 1;
            app.ControlPanel.Layout.Row = 1;
            
            ctrlGrid = uigridlayout(app.ControlPanel);
            ctrlGrid.ColumnWidth = {'1x', '1x'};
            ctrlGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', '2x', 'fit'}; 
            
            app.LengthLabel = uilabel(ctrlGrid);
            app.LengthLabel.Layout.Row = 1;
            app.LengthLabel.Layout.Column = 1;
            app.LengthLabel.Text = 'Length (L / \lambda):';
            app.LengthLabel.FontWeight = 'bold';
            
            app.LengthEditField = uieditfield(ctrlGrid, 'numeric');
            app.LengthEditField.Layout.Row = 1;
            app.LengthEditField.Layout.Column = 2;
            app.LengthEditField.Value = 0.5;
            app.LengthEditField.Limits = [0.01 10];
            app.ModelLabel = uilabel(ctrlGrid);
            app.ModelLabel.Layout.Row = 2;
            app.ModelLabel.Layout.Column = 1;
            app.ModelLabel.Text = 'Current Model:';
            app.ModelLabel.FontWeight = 'bold';
            
            app.ModelDropDown = uidropdown(ctrlGrid);
            app.ModelDropDown.Layout.Row = 2;
            app.ModelDropDown.Layout.Column = 2;
            app.ModelDropDown.Items = {'General (Sinusoidal)', 'Short (Triangular)', 'Infinitesimal (Uniform)'};
            app.ModelDropDown.Value = 'General (Sinusoidal)';
            app.ShowHPBWCheckBox = uicheckbox(ctrlGrid);
            app.ShowHPBWCheckBox.Layout.Row = 3;
            app.ShowHPBWCheckBox.Layout.Column = [1 2]; 
            app.ShowHPBWCheckBox.Text = 'Show -3dB Beamwidth Markers';
            app.ShowHPBWCheckBox.Value = true;
            app.CalculateButton = uibutton(ctrlGrid, 'push');
            app.CalculateButton.Layout.Row = 4;
            app.CalculateButton.Layout.Column = [1 2];
            app.CalculateButton.Text = 'Calculate Pattern';
            app.CalculateButton.FontWeight = 'bold';
            app.CalculateButton.ButtonPushedFcn = @(btn, event) CalculateButtonPushed(app, event);
            app.ResultsPanel = uipanel(ctrlGrid);
            app.ResultsPanel.Layout.Row = 6;
            app.ResultsPanel.Layout.Column = [1 2];
            app.ResultsPanel.Title = 'Calculated Metrics';
            
            resGrid = uigridlayout(app.ResultsPanel);
            resGrid.ColumnWidth = {'1x'};
            resGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
            
            app.DirectivityLabel = uilabel(resGrid);
            app.DirectivityLabel.Text = 'Directivity: --';
            app.DirectivityLabel.FontWeight = 'bold';
            
            app.RradLabel = uilabel(resGrid);
            app.RradLabel.Text = 'R_rad: --';
            
            app.HPBWLabel = uilabel(resGrid);
            app.HPBWLabel.Text = 'HPBW: --';
            
            app.WarningLabel = uilabel(resGrid);
            app.WarningLabel.Text = '';
            app.WarningLabel.FontAngle = 'italic';
            
            app.VisualTabGroup = uitabgroup(app.MainGridLayout);
            app.VisualTabGroup.Layout.Column = 2;
            app.VisualTabGroup.Layout.Row = 1;
            
            app.PatternTab = uitab(app.VisualTabGroup);
            app.PatternTab.Title = 'Radiation Pattern (2D/3D)';
            patGrid = uigridlayout(app.PatternTab);
            patGrid.ColumnWidth = {'1x', '1x'};
            patGrid.RowHeight = {'1x'};
            app.PolarAxes = polaraxes('Parent', patGrid);
            app.PolarAxes.Layout.Column = 1;
            app.ThreeDAxes = uiaxes(patGrid);
            app.ThreeDAxes.Layout.Column = 2;
            app.CurrentTab = uitab(app.VisualTabGroup);
            app.CurrentTab.Title = 'Current Distribution';
            app.CurrentAxes = uiaxes(app.CurrentTab);
            app.CurrentAxes.Position = [50 50 500 400];
            app.InfoTab = uitab(app.VisualTabGroup);
            app.InfoTab.Title = 'Analysis Details';
            app.FormulaTextArea = uitextarea(app.InfoTab);
            app.FormulaTextArea.Position = [20 20 600 500];
            app.FormulaTextArea.Editable = 'off';
            app.FormulaTextArea.FontName = 'Monospaced';
            app.CalculateButtonPushed(0);
            app.UIFigure.Visible = 'on';
        end
    end
    methods (Access = public)
        function app = LinearDipoleAnalyzer
            createComponents(app)
            registerApp(app, app.UIFigure)
            if nargout == 0, clear app; end
        end
        function delete(app)
            delete(app.UIFigure);
        end
    end
end