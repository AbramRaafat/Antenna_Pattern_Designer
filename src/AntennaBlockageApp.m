classdef AntennaBlockageApp < matlab.apps.AppBase
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        GridLayout              matlab.ui.container.GridLayout
        VisualsPanel            matlab.ui.container.Panel
        TabGroup                matlab.ui.container.TabGroup
        GeoTab                  matlab.ui.container.Tab
        GeoAxes                 matlab.ui.control.UIAxes
        HPlaneTab               matlab.ui.container.Tab
        HPlaneAxes              matlab.ui.control.UIAxes
        EPlaneTab               matlab.ui.container.Tab
        EPlaneAxes              matlab.ui.control.UIAxes
        PolarTab                matlab.ui.container.Tab
        PolarAxes               matlab.ui.control.UIAxes
        DpatternTab             matlab.ui.container.Tab
        ThreeDAxes              matlab.ui.control.UIAxes
        parametersPanel         matlab.ui.container.Panel
        EffAreaLabel            matlab.ui.control.Label
        DirectivityLabel        matlab.ui.control.Label
        OffsetYField            matlab.ui.control.NumericEditField
        OffsetYdyLabel          matlab.ui.control.Label
        OffsetXField            matlab.ui.control.NumericEditField
        OffsetXdxLabel          matlab.ui.control.Label
        ResultsLabel            matlab.ui.control.Label
        SLLLabel                matlab.ui.control.Label
        GainLossLabel           matlab.ui.control.Label
        CalculateButton         matlab.ui.control.Button
        BlockHeightField        matlab.ui.control.NumericEditField
        BlockHcyEditFieldLabel  matlab.ui.control.Label
        BlockWidthField         matlab.ui.control.NumericEditField
        BlockWcxEditFieldLabel  matlab.ui.control.Label
        ApertureHeightField     matlab.ui.control.NumericEditField
        HieghtbEditFieldLabel   matlab.ui.control.Label
        ApertureWidthField      matlab.ui.control.NumericEditField
        WidthaEditFieldLabel    matlab.ui.control.Label
        ScanAngleField          matlab.ui.control.NumericEditField
        ScanAngleLabel          matlab.ui.control.Label
    end
    
    methods (Access = private)
        function CalculateButtonPushed(app, event)
            a  = app.ApertureWidthField.Value;
            b  = app.ApertureHeightField.Value;
            cx = app.BlockWidthField.Value;
            cy = app.BlockHeightField.Value;
            dx = app.OffsetXField.Value; 
            dy = app.OffsetYField.Value; 
            scan_deg = app.ScanAngleField.Value;
            
            % Scan Calculation [Lecture Page 26, Eq 67]
            % gamma = beta * sin(theta_0)
            u_scan = sin(deg2rad(scan_deg));
            
            Area_Ideal = a * b;
            
            % Field Calculations
            theta = linspace(-pi/2, pi/2, 1000); 
            u_plot = sin(theta); % u = sin(theta)
            
            % Obliquity Factor [Lecture Page 8, Eq 21]
            % E_field ~ (1 + cos(theta))
            obliquity = (1 + cos(theta)) / 2;
            
            % H-PLANE CUT (Scanning Plane)
            % [Lecture Page 26, Eq 66]-> epsilon(u) = a * Sa((u - gamma)a/2)

            E_main_H = (a*b) * sinc(a * (u_plot - u_scan)); 
            
            % Blockage Term [Lecture Page 22, Eq 59]
            phase_block_H = exp(1i * 2 * pi * dx * (u_plot - u_scan));
            E_block_H = (cx*cy) * sinc(cx * (u_plot - u_scan)) .* phase_block_H;
            
            % Total Pattern
            pat_lin_H = abs(E_main_H - E_block_H) / Area_Ideal;
            pat_lin_H = pat_lin_H .* obliquity; 
            
            % --- E-PLANE CUT (Transverse Plane) ---
            % note for fixing the E-PLANE when scanning by shifting the
            % plan to the scane angle u = u_scan not 0
            
            E_main_E = (a*b) * sinc(b * u_plot); 
            
            phase_block_E = exp(1i * 2 * pi * dy * u_plot);
            E_block_E = (cx*cy) * sinc(cy * u_plot) .* phase_block_E;
            
            pat_lin_E = abs(E_main_E - E_block_E) / Area_Ideal;
            pat_lin_E = pat_lin_E .* obliquity;
            
            %%% PLotting 
            
            % Rectangular Plots
            plot(app.HPlaneAxes, u_plot, pat_lin_H, 'LineWidth', 2, 'Color', 'r');
            title(app.HPlaneAxes, sprintf('H-Plane Cut (Through Scan Peak %.1f^o)', scan_deg));
            xlabel(app.HPlaneAxes, 'sin(\theta)'); grid(app.HPlaneAxes, 'on'); ylim(app.HPlaneAxes, [0 1]);
            xline(app.HPlaneAxes, u_scan, '--k', 'Scan Angle'); 
            
            plot(app.EPlaneAxes, u_plot, pat_lin_E, 'LineWidth', 2, 'Color', 'b');
            title(app.EPlaneAxes, 'E-Plane Cut (Orthogonal through Peak)');
            xlabel(app.EPlaneAxes, 'sin(\theta)'); grid(app.EPlaneAxes, 'on'); ylim(app.EPlaneAxes, [0 1]);

            % Polar Plot 
            ang_H = pi/2 - theta; 
            xh = pat_lin_H .* cos(ang_H); 
            yh = pat_lin_H .* sin(ang_H);
            
            
            ang_E = pi/2 - theta;
            xe = pat_lin_E .* cos(ang_E); 
            ye = pat_lin_E .* sin(ang_E);
            
            cla(app.PolarAxes); hold(app.PolarAxes, 'on');
            
            plot(app.PolarAxes, xh, yh, 'r', 'LineWidth', 2);
            % note for the actual 3D this will be rotated to, but for clear
            % side loop to main loop we project it across the scanning
            % angle 
            
            plot(app.PolarAxes, xe, ye, 'b', 'LineWidth', 2);
            
            legend(app.PolarAxes, 'H-Plane (Scanned)', 'E-Plane (Profile)');
            title(app.PolarAxes, 'Polar Cuts (Comparison)');
            axis(app.PolarAxes, 'equal'); 
            xlim(app.PolarAxes,[-1 1]); ylim(app.PolarAxes,[0 1]); 
            grid(app.PolarAxes, 'on'); 
            
            % 3D Plot
            [Theta_3D, Phi_3D] = meshgrid(linspace(0, pi/2, 60), linspace(0, 2*pi, 60));
            U_3D = sin(Theta_3D) .* cos(Phi_3D);
            V_3D = sin(Theta_3D) .* sin(Phi_3D);
            obliquity_3D = (1 + cos(Theta_3D))/2;
            
            % Apply Scan Shift to U [Lecture Eq 66 Logic]
            U_Prime = U_3D - u_scan; 
            
            E_main_3D = (a*b) * sinc(a * U_Prime) .* sinc(b * V_3D);
            Phase_3D = 2*pi*(dx*U_Prime + dy*V_3D);
            E_block_3D = (cx*cy) * sinc(cx * U_Prime) .* sinc(cy * V_3D) .* exp(1i*Phase_3D);
            
            Mag_3D = (abs(E_main_3D - E_block_3D) / Area_Ideal) .* obliquity_3D;
            
            R = Mag_3D;
            X3 = R .* sin(Theta_3D) .* cos(Phi_3D);
            Y3 = R .* sin(Theta_3D) .* sin(Phi_3D);
            Z3 = R .* cos(Theta_3D);
            
            cla(app.ThreeDAxes);
            surf(app.ThreeDAxes, X3, Y3, Z3, Mag_3D, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            hold(app.ThreeDAxes, 'on');
            
            theta_cut = linspace(-pi/2, pi/2, 200);
            u_cut = sin(theta_cut);
            obl_cut = (1+cos(theta_cut))/2;
            E_H_val = (abs( ((a*b)*sinc(a*(u_cut - u_scan))) - ...
                       ((cx*cy)*sinc(cx*(u_cut - u_scan)).*exp(1i*2*pi*dx*(u_cut - u_scan))) ) / Area_Ideal) .* obl_cut;
            
            XH_line = E_H_val .* sin(theta_cut);
            ZH_line = E_H_val .* cos(theta_cut);
            plot3(app.ThreeDAxes, XH_line, zeros(size(XH_line)), ZH_line, 'r', 'LineWidth', 3);
            

            E_E_val = (abs( ((a*b)*sinc(b*u_cut)) - ...
                       ((cx*cy)*sinc(cy*u_cut).*exp(1i*2*pi*dy*u_cut)) ) / Area_Ideal) .* obl_cut;
            
            
            u_line = u_scan * ones(size(theta_cut));
            v_line = u_cut; 
            
            xE_flat = zeros(size(E_E_val));
            yE_flat = E_E_val .* sin(theta_cut);
            zE_flat = E_E_val .* cos(theta_cut);
            
            % Rotation Matrix for rotation around Y 
            % this took me to time to figure so :_ )

            % [cos th  0  sin th]
            % [   0    1    0   ]
            % [-sin th 0  cos th]
            th_rot = deg2rad(scan_deg);
            xE_rot = xE_flat * cos(th_rot) + zE_flat * sin(th_rot);
            yE_rot = yE_flat;
            zE_rot = -xE_flat * sin(th_rot) + zE_flat * cos(th_rot);
            
            plot3(app.ThreeDAxes, xE_rot, yE_rot, zE_rot, 'b', 'LineWidth', 3);
            
            view(app.ThreeDAxes, 3); axis(app.ThreeDAxes, 'equal'); 
            
            %%% Results
            Area_Net = (a*b - cx*cy);
            if Area_Net > 0
                gain_loss_db = 20 * log10(Area_Net / Area_Ideal);
            else
                gain_loss_db = -inf;
            end
            app.GainLossLabel.Text = sprintf("Gain Degradation: %.2f dB", gain_loss_db);
            
            % SLL
            [pks, ~] = findpeaks(pat_lin_H);
            pks_sorted = sort(pks, 'descend');
            if length(pks_sorted) >= 2
                sll_val = 20 * log10(pks_sorted(2) / pks_sorted(1));
                app.SLLLabel.Text = sprintf("Max SLL (H-Plane): %.2f dB", sll_val);
            else
                app.SLLLabel.Text = "Max SLL: N/A";
            end
            
            app.EffAreaLabel.Text = sprintf("Effective Area: %.2f sq.lam", Area_Net);
            if Area_Net > 0
                app.DirectivityLabel.Text = sprintf("Directivity: %.2f dBi", 10 * log10(4 * pi * Area_Net));
            else
                app.DirectivityLabel.Text = "Directivity: N/A";
            end

            % antenna ploting 
            cla(app.GeoAxes); hold(app.GeoAxes, 'on');
            rectangle(app.GeoAxes, 'Position', [-a/2, -b/2, a, b], 'FaceColor', [0.8 0.9 1], 'EdgeColor', 'b');
            rectangle(app.GeoAxes, 'Position', [dx-cx/2, dy-cy/2, cx, cy], 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'r');
            if abs(u_scan) > 0.01
                quiver(app.GeoAxes, 0, 0, u_scan*a/2, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);
            end
            axis(app.GeoAxes, 'equal'); grid(app.GeoAxes, 'on');
            limit_val = max([a, b]) / 2 + 1;
            xlim(app.GeoAxes, [-limit_val, limit_val]); ylim(app.GeoAxes, [-limit_val, limit_val]);
        end
    end
    
    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'Aperture Blockage App';
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1.4x', '3.6x'};
            app.GridLayout.RowHeight = {'3x'};
            app.parametersPanel = uipanel(app.GridLayout);
            app.parametersPanel.Layout.Row = 1;
            app.parametersPanel.Layout.Column = 1;
            
            % Parameters
            lbl_pos = [2 411 54 22; 0 375 58 22; 0 335 72 22; 0 295 70 22; 4 254 72 22; 4 219 72 22; 4 185 80 22];
            field_pos = [70 411 100 22; 68 375 100 22; 68 335 100 22; 68 295 100 22; 72 254 100 22; 72 219 100 22; 90 185 82 22];
            labels = {'Width (a)', 'Height (b)', 'Block W (cx)', 'Block H (cy)', 'Offset X (dx)', 'Offset Y (dy)', 'Scan Ang'};
            
            for i = 1:7
               l = uilabel(app.parametersPanel); l.Position = lbl_pos(i,:); l.Text = labels{i};
               f = uieditfield(app.parametersPanel, 'numeric'); f.Position = field_pos(i,:); 
               if i==1, app.ApertureWidthField=f; f.Value=3;
               elseif i==2, app.ApertureHeightField=f; f.Value=4;
               elseif i==3, app.BlockWidthField=f; f.Value=1;
               elseif i==4, app.BlockHeightField=f; f.Value=4;
               elseif i==5, app.OffsetXField=f;
               elseif i==6, app.OffsetYField=f;
               elseif i==7, app.ScanAngleField=f; f.Value=0; end
            end
            
            app.CalculateButton = uibutton(app.parametersPanel, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Position = [1 150 170 25];
            app.CalculateButton.Text = 'Calculate';
            
            app.GainLossLabel = uilabel(app.parametersPanel); app.GainLossLabel.Position = [2 100 158 22];
            app.SLLLabel = uilabel(app.parametersPanel); app.SLLLabel.Position = [2 75 158 22];
            app.DirectivityLabel = uilabel(app.parametersPanel); app.DirectivityLabel.Position = [4 50 147 22];
            app.EffAreaLabel = uilabel(app.parametersPanel); app.EffAreaLabel.Position = [6 25 145 22];

            app.VisualsPanel = uipanel(app.GridLayout);
            app.VisualsPanel.Layout.Row = 1; app.VisualsPanel.Layout.Column = 2;
            app.TabGroup = uitabgroup(app.VisualsPanel);
            app.TabGroup.Position = [13 29 410 396];
            
            tabs = {'Geometry', 'H-Plane', 'E-Plane', 'Polar', '3D Pattern'};
            axes_props = {'GeoAxes', 'HPlaneAxes', 'EPlaneAxes', 'PolarAxes', 'ThreeDAxes'};
            for i=1:5
                t = uitab(app.TabGroup); t.Title = tabs{i};
                ax = uiaxes(t); ax.Position = [16 11 378 336];
                if i==5, ax.Position = [16 37 378 310]; end
                app.(axes_props{i}) = ax;
            end
            
            app.UIFigure.Visible = 'on';
        end
    end
    
    methods (Access = public)
        function app = AntennaBlockageApp
            createComponents(app)
            registerApp(app, app.UIFigure)
            if nargout == 0, clear app; end
        end
        function delete(app), delete(app.UIFigure); end
    end
end