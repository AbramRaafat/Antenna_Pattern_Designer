classdef LinearArrayDesigner < matlab.apps.AppBase

    properties (Access = public)
        UIFigure              matlab.ui.Figure
        MainGridLayout        matlab.ui.container.GridLayout
        ControlPanel          matlab.ui.container.Panel
        VisualTabGroup        matlab.ui.container.TabGroup
        CartesianTab          matlab.ui.container.Tab
        PolarTab              matlab.ui.container.Tab
        ThreeDTab             matlab.ui.container.Tab
        MetricsPanel          matlab.ui.container.Panel
        CartesianAxes         matlab.ui.control.UIAxes
        PolarAxes             matlab.graphics.axis.PolarAxes
        ThreeDAxes            matlab.ui.control.UIAxes
        TypeLabel             matlab.ui.control.Label
        TypeDropDown          matlab.ui.control.DropDown
        NLabel                matlab.ui.control.Label
        NEditField            matlab.ui.control.NumericEditField
        DLabel                matlab.ui.control.Label
        DEditField            matlab.ui.control.EditField
        ScanLabel             matlab.ui.control.Label
        ScanSlider            matlab.ui.control.Slider
        ScanValueLabel        matlab.ui.control.Label
        AlphaLabel            matlab.ui.control.Label
        AlphaEditField        matlab.ui.control.EditField
        SLLLabel              matlab.ui.control.Label
        SLLEditField          matlab.ui.control.NumericEditField
        CalcButton            matlab.ui.control.Button
        DirLabel              matlab.ui.control.Label
        FNBWLabel             matlab.ui.control.Label
        MeasSLLLabel          matlab.ui.control.Label
        TargetSLLResLabel     matlab.ui.control.Label
        WarningLabel          matlab.ui.control.Label
        NoteLabel             matlab.ui.control.Label
    end
    

    methods (Access = private)
        
        function UpdateAnalysis(app, sourceOfChange)
            % sourceOfChange: 'Type', 'N', 'd', 'Slider', 'AlphaBox', 'SLL'
            
            app.CalcButton.Text = 'Calculating...';
            app.CalcButton.Enable = 'off';
            drawnow;
            
            try
                % Gather Inputs
                N = app.NEditField.Value;
                type = app.TypeDropDown.Value;
                
                d_str = app.DEditField.Value;
                d_lam = str2num(d_str); %#ok<ST2NM>
                if isempty(d_lam), d_lam = 0.5; end 
                
                % Handle SLL Visibility
                if strcmp(type, 'Dolph-Chebyshev')
                    app.SLLEditField.Enable = 'on';
                    SLL_target = app.SLLEditField.Value;
                    app.TargetSLLResLabel.Text = sprintf('Target SLL: %.1f dB', SLL_target);
                else
                    app.SLLEditField.Enable = 'off';
                    SLL_target = -20;
                    app.TargetSLLResLabel.Text = 'Target SLL: N/A';
                end
                
                % Constants
                k = 2*pi;
                bd = k * d_lam;
                
                % ALPHA & SCAN LOGIC
                alpha = 0; 
                
                if strcmp(type, 'Hansen-Woodward')
                    % Hansen overrides everything -> Endfire
                    alpha = -(bd + pi/N);
                    app.AlphaEditField.Value = num2str(alpha, '%.4f');
                    app.ScanSlider.Value = 0;
                    app.ScanValueLabel.Text = '0 deg (Forced)';
                else
                    if strcmp(sourceOfChange, 'Slider')
                        theta_deg = app.ScanSlider.Value;
                        alpha = -bd * cosd(theta_deg);
                        app.AlphaEditField.Value = num2str(alpha, '%.4f');
                        app.ScanValueLabel.Text = sprintf('%d deg', round(theta_deg));
                        
                    elseif strcmp(sourceOfChange, 'AlphaBox')
                        alpha_str = app.AlphaEditField.Value;
                        alpha = str2num(alpha_str); %#ok<ST2NM>
                        if isempty(alpha), alpha = 0; end
                        
                        % Update Slider based on Alpha
                        if bd > 1e-9
                            cos_th = -alpha / bd;
                            if abs(cos_th) <= 1.0001
                                if cos_th > 1, cos_th = 1; end
                                if cos_th < -1, cos_th = -1; end
                                theta_deg = acosd(cos_th);
                                app.ScanSlider.Value = theta_deg;
                                app.ScanValueLabel.Text = sprintf('%d deg', round(theta_deg));
                            else
                                app.ScanValueLabel.Text = 'Invisible';
                            end
                        end
                        
                    else
                        % Logic to prevent "Stuck at 0" when switching from Hansen
                        if app.ScanSlider.Value < 0.1 && ~strcmp(sourceOfChange, 'Slider')
                            app.ScanSlider.Value = 90;
                            app.ScanValueLabel.Text = '90 deg';
                        end
                        theta_deg = app.ScanSlider.Value;
                        alpha = -bd * cosd(theta_deg);
                        app.AlphaEditField.Value = num2str(alpha, '%.4f');
                    end
                end
                
                % Weights Calculation
                switch type
                    case 'Uniform'
                        weights = ones(1, N);
                    case 'Binomial'
                        if N < 20
                            weights = zeros(1,N);
                            for i=0:N-1
                                weights(i+1) = nchoosek(N-1, i);
                            end
                        else
                            n = 0:N-1;
                            weights = exp(-(n - (N-1)/2).^2 / (0.5*(N-1)));
                        end
                    case 'Hansen-Woodward'
                        weights = ones(1, N);
                    case 'Dolph-Chebyshev'
                        weights = app.ComputeChebyshevWeights(N, SLL_target);
                end
                
                % Compute Array Factor & Metrics
                
                % A. Physical Grid (0 to 180 deg) for Metrics & Cartesian
                theta_phys = linspace(0, pi, 1000); 
                psi_phys = bd * cos(theta_phys) + alpha;
                
                AF_phys = zeros(size(theta_phys));
                for n = 1:N
                    AF_phys = AF_phys + weights(n) * exp(1i * (n-1) * psi_phys);
                end
                
                % Normalize
                max_AF = max(abs(AF_phys));
                if max_AF > 1e-9
                    AF_phys_norm = abs(AF_phys) / max_AF;
                else
                    AF_phys_norm = abs(AF_phys);
                end

                % ================== metric =====================  
                
                % Directivity
                integrand = (AF_phys_norm.^2) .* sin(theta_phys);
                denom = trapz(theta_phys, integrand);
                if denom > 0
                    D = 2 / denom;
                    app.DirLabel.Text = sprintf('Directivity: %.2f dBi', 10*log10(D));
                else
                    app.DirLabel.Text = 'Directivity: N/A';
                end
                
                % Measured SLL
                [pks, ~] = findpeaks(AF_phys_norm);
                if AF_phys_norm(1) > AF_phys_norm(2), pks = [pks, AF_phys_norm(1)]; end
                if AF_phys_norm(end) > AF_phys_norm(end-1), pks = [pks, AF_phys_norm(end)]; end
                
                pks = sort(unique(round(pks, 4)), 'descend');
                if length(pks) >= 2
                    meas_sll = 20*log10(pks(2));
                    if meas_sll < -100, meas_sll = -Inf; end
                    app.MeasSLLLabel.Text = sprintf('Meas. SLL: %.1f dB', meas_sll);
                else
                    app.MeasSLLLabel.Text = 'Meas. SLL: None';
                end
                
                % FNBW 
                [~, max_idx] = max(AF_phys_norm);
                idx_L = max_idx;
                while idx_L > 1 && AF_phys_norm(idx_L) > 0.05, idx_L = idx_L - 1; end
                idx_R = max_idx;
                while idx_R < length(theta_phys) && AF_phys_norm(idx_R) > 0.05, idx_R = idx_R + 1; end
                
                fnbw = rad2deg(theta_phys(idx_R) - theta_phys(idx_L));
                app.FNBWLabel.Text = sprintf('FNBW (Approx): %.1f deg', fnbw);
                
                % Visible region: psi is in [alpha - bd, alpha + bd]
                % We check if any integer multiple of 2*pi (other than 0) falls in this range.
                
                psi_min = alpha - bd;
                psi_max = alpha + bd;
                
                % Normalize by 2*pi to find integers
                min_k = ceil(psi_min / (2*pi));
                max_k = floor(psi_max / (2*pi));
                
                has_GL = false;
                
                % If the range includes any integer, we must check if it's ONLY zero.
                if max_k >= min_k
                    % We have integers in range Check if any is non-zero.
                    % If min_k > 0 or max_k < 0, then 0 is definitely not the only one.
                    % If min_k <= 0 <= max_k, we have 0 plus potentially others.
                    
                    if min_k > 0 || max_k < 0
                        has_GL = true;
                    elseif min_k < 0 || max_k > 0
                        has_GL = true; % Zero is included, but so are others (-1 or 1)
                    end
                end
                
                if has_GL
                    app.WarningLabel.Text = 'WARNING: Grating Lobes!';
                    app.WarningLabel.FontColor = [1 0 0]; % Red
                else
                    app.WarningLabel.Text = 'Status: OK (No GL)';
                    app.WarningLabel.FontColor = [0 0.5 0]; % Green
                end

                % Visualization
                
                cla(app.CartesianAxes);
                hold(app.CartesianAxes, 'on');
                grid(app.CartesianAxes, 'on');
                
                if strcmp(type, 'Dolph-Chebyshev')
                    

                    psi_math = linspace(-2*pi, 2*pi, 2000);
                    
                    % Array Factor
                    AF_math = zeros(size(psi_math));
                    for n = 1:N
                        AF_math = AF_math + weights(n) * exp(1i * (n-1) * psi_math);
                    end
                    
                    % Normalize 
                    max_AF_math = max(abs(AF_math));
                    AF_math_norm = abs(AF_math) / max_AF_math;
                    
                    % Plot the Invisible/Mathematic
                    x_math = cos(psi_math / 2);
                    [x_math_sorted, idx_m] = sort(x_math);
                    y_math_sorted = AF_math_norm(idx_m);
                    
                    plot(app.CartesianAxes, x_math_sorted, y_math_sorted, ...
                        'LineWidth', 1.5, 'Color', [0.7 0.7 0.7], 'LineStyle', '--', ...
                        'DisplayName', 'Invisible Region');
                        
                    % Re-normalize physical AF to the same global max
                    x_vals = cos(psi_phys / 2);
                    AF_phys_plot = abs(AF_phys) / max_AF_math; % Use global max [VIP]
                    
                    [x_sorted, idx] = sort(x_vals);
                    y_sorted = AF_phys_plot(idx);
                    
                    plot(app.CartesianAxes, x_sorted, y_sorted, ...
                        'LineWidth', 2, 'Color', [0 0.45 0.74], ...
                        'DisplayName', 'Visible Region');
                    
                    xlabel(app.CartesianAxes, 'Z / Z_0  ( \approx cos(\psi/2) )');
                    ylabel(app.CartesianAxes, 'Normalized |AF|');
                    title(app.CartesianAxes, 'Chebyshev Polynomial View');
                    xlim(app.CartesianAxes, [-1 1]); 
                    yline(app.CartesianAxes, 10^(SLL_target/20), ':r', 'Target SLL');
                    legend(app.CartesianAxes, 'Location', 'northeast');
                    
                else

                    plot(app.CartesianAxes, rad2deg(theta_phys), AF_phys_norm, 'LineWidth', 2, 'Color', [0 0.45 0.74]);
                    xlabel(app.CartesianAxes, '\theta (Degrees)');
                    ylabel(app.CartesianAxes, 'Normalized |AF|');
                    title(app.CartesianAxes, 'Array Factor (Cartesian)');
                    xlim(app.CartesianAxes, [0 180]);
                end
                
                ylim(app.CartesianAxes, [0 1.05]);
                hold(app.CartesianAxes, 'off');
                
                % Polar Plot 
                theta_full = linspace(0, 2*pi, 1000);
                psi_full = bd * cos(theta_full) + alpha;
                AF_full = zeros(size(theta_full));
                for n = 1:N
                    AF_full = AF_full + weights(n) * exp(1i * (n-1) * psi_full);
                end
                AF_full_norm = abs(AF_full) / max(abs(AF_full));
                
                cla(app.PolarAxes);
                polarplot(app.PolarAxes, theta_full, AF_full_norm, 'LineWidth', 2, 'Color', 'r');
                app.PolarAxes.ThetaZeroLocation = 'top';
                app.PolarAxes.ThetaDir = 'clockwise';
                
                % 3D Plot
                phi_3d = linspace(0, 2*pi, 60);
                [TH, PH] = meshgrid(theta_phys, phi_3d);
                R_3d = repmat(AF_phys_norm, 60, 1);
                X = R_3d .* sin(TH) .* cos(PH);
                Y = R_3d .* sin(TH) .* sin(PH);
                Z = R_3d .* cos(TH);
                
                cla(app.ThreeDAxes);
                surf(app.ThreeDAxes, X, Y, Z, R_3d, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
                view(app.ThreeDAxes, 3);
                axis(app.ThreeDAxes, 'equal');
                colormap(app.ThreeDAxes, 'jet');
                title(app.ThreeDAxes, '3D Pattern');
                
            catch ME
                disp(['ERROR: ' ME.message]);
                app.WarningLabel.Text = 'Error in Calculation';
                app.WarningLabel.FontColor = [1 0 0];
            end
            
            app.CalcButton.Text = 'Update Array';
            app.CalcButton.Enable = 'on';
            drawnow;
        end
        
        function w = ComputeChebyshevWeights(~, N, SLL_dB)
            try 
                if exist('chebwin', 'file')
                    w = chebwin(N, abs(SLL_dB)).'; 
                    w = w / max(w);
                    return;
                end
            catch
            end
            
            % Fallback (my old method that work but for some reason the app freeze when I use it)
            try
                R0 = 10^(abs(SLL_dB)/20);
                M = N - 1;
                val_acosh = log(R0 + sqrt(R0^2 - 1));
                z0 = cosh(val_acosh / M);
                k_idx = 0:N-1;
                psi_k = 2*pi*k_idx / N;
                x_k = z0 * cos(psi_k / 2);
                AF_sample = zeros(size(x_k));
                for i = 1:length(x_k)
                    x = x_k(i);
                    if abs(x) <= 1
                        AF_sample(i) = cos(M * acos(x));
                    else
                        AF_sample(i) = cosh(M * log(abs(x) + sqrt(x^2 - 1))) * sign(x)^M;
                    end
                end
                w = real(ifft(AF_sample));
                w = fftshift(w);
                w = w / max(abs(w));
            catch
                w = ones(1,N); 
            end
        end
        
    end
    

    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1000 680];
            app.UIFigure.Name = 'App 2: Linear Array Designer';
            app.MainGridLayout = uigridlayout(app.UIFigure);
            app.MainGridLayout.ColumnWidth = {280, '1x'};
            app.MainGridLayout.RowHeight = {'1x'};

            app.ControlPanel = uipanel(app.MainGridLayout);
            app.ControlPanel.Title = 'Array Synthesis';
            app.ControlPanel.Layout.Column = 1;
            cGrid = uigridlayout(app.ControlPanel);
            cGrid.ColumnWidth = {'1x', '1x'};
            cGrid.RowHeight = {'fit','fit', 'fit','fit', 'fit','fit', 'fit','fit', 'fit','fit', 'fit','fit', '2x', 'fit'};

            app.TypeLabel = uilabel(cGrid); app.TypeLabel.Text = 'Array Type:'; app.TypeLabel.FontWeight = 'bold';
            app.TypeDropDown = uidropdown(cGrid);
            app.TypeDropDown.Items = {'Uniform', 'Hansen-Woodward', 'Dolph-Chebyshev', 'Binomial'};
            app.TypeDropDown.ValueChangedFcn = @(src, e) app.UpdateAnalysis('Type');

            app.NLabel = uilabel(cGrid); app.NLabel.Text = 'Elements (N):';
            app.NEditField = uieditfield(cGrid, 'numeric');
            app.NEditField.Value = 10;
            app.NEditField.Limits = [2 200];
            app.NEditField.ValueChangedFcn = @(src, e) app.UpdateAnalysis('N');

            app.DLabel = uilabel(cGrid); app.DLabel.Text = 'Spacing (d/\lambda):';
            app.DEditField = uieditfield(cGrid, 'text'); 
            app.DEditField.Value = '0.5';
            app.DEditField.ValueChangedFcn = @(src, e) app.UpdateAnalysis('d');

            app.ScanLabel = uilabel(cGrid); app.ScanLabel.Text = 'Scan Angle:';
            app.ScanValueLabel = uilabel(cGrid); app.ScanValueLabel.Text = '90 deg'; app.ScanValueLabel.HorizontalAlignment = 'right';
            app.ScanSlider = uislider(cGrid);
            app.ScanSlider.Layout.Column = [1 2];
            app.ScanSlider.Limits = [0 180];
            app.ScanSlider.Value = 90;
            app.ScanSlider.ValueChangedFcn = @(src, e) app.UpdateAnalysis('Slider');

            app.AlphaLabel = uilabel(cGrid); app.AlphaLabel.Text = 'Phase Shift (\alpha):';
            app.AlphaEditField = uieditfield(cGrid, 'text'); 
            app.AlphaEditField.Value = '0';
            app.AlphaEditField.ValueChangedFcn = @(src, e) app.UpdateAnalysis('AlphaBox');

            app.SLLLabel = uilabel(cGrid); app.SLLLabel.Text = 'Target SLL (dB):';
            app.SLLEditField = uieditfield(cGrid, 'numeric');
            app.SLLEditField.Value = -20;
            app.SLLEditField.Limits = [-100 -5];
            app.SLLEditField.ValueChangedFcn = @(src, e) app.UpdateAnalysis('SLL');

            app.CalcButton = uibutton(cGrid, 'push');
            app.CalcButton.Layout.Column = [1 2];
            app.CalcButton.Text = 'Update Array';
            app.CalcButton.FontWeight = 'bold';
            app.CalcButton.ButtonPushedFcn = @(src, e) app.UpdateAnalysis('Button');

            app.MetricsPanel = uipanel(cGrid);
            app.MetricsPanel.Layout.Row = 14;
            app.MetricsPanel.Layout.Column = [1 2];
            app.MetricsPanel.Title = 'Design Metrics';
            mGrid = uigridlayout(app.MetricsPanel);
            mGrid.ColumnWidth = {'1x'};
            mGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            
            app.WarningLabel = uilabel(mGrid); app.WarningLabel.Text = 'Status: OK'; app.WarningLabel.FontWeight = 'bold';
            app.DirLabel = uilabel(mGrid); app.DirLabel.Text = 'Directivity: --';
            app.FNBWLabel = uilabel(mGrid); app.FNBWLabel.Text = 'FNBW: --';
            app.MeasSLLLabel = uilabel(mGrid); app.MeasSLLLabel.Text = 'Meas. SLL: --';
            app.TargetSLLResLabel = uilabel(mGrid); app.TargetSLLResLabel.Text = 'Target SLL: --'; app.TargetSLLResLabel.FontAngle = 'italic';
            app.NoteLabel = uilabel(mGrid); app.NoteLabel.Text = 'Note: --'; app.NoteLabel.FontAngle = 'italic'; app.NoteLabel.FontSize = 10;

            app.VisualTabGroup = uitabgroup(app.MainGridLayout);
            
            app.CartesianTab = uitab(app.VisualTabGroup); app.CartesianTab.Title = 'Cartesian';
            cg = uigridlayout(app.CartesianTab); cg.ColumnWidth = {'1x'}; cg.RowHeight = {'1x'};
            app.CartesianAxes = uiaxes(cg);
            
            app.PolarTab = uitab(app.VisualTabGroup); app.PolarTab.Title = 'Polar Pattern';
            pg = uigridlayout(app.PolarTab); pg.ColumnWidth = {'1x'}; pg.RowHeight = {'1x'};
            app.PolarAxes = polaraxes('Parent', pg);
            
            app.ThreeDTab = uitab(app.VisualTabGroup); app.ThreeDTab.Title = '3D Revolved';
            tg = uigridlayout(app.ThreeDTab); tg.ColumnWidth = {'1x'}; tg.RowHeight = {'1x'};
            app.ThreeDAxes = uiaxes(tg);
            
            app.UIFigure.Visible = 'on';
            drawnow;
            app.UpdateAnalysis('Init');
        end
    end
    
    methods (Access = public)
        function app = LinearArrayDesigner
            createComponents(app)
            registerApp(app, app.UIFigure)
            if nargout == 0, clear app; end
        end
        function delete(app)
            delete(app.UIFigure);
        end
    end
end