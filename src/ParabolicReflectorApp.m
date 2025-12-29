classdef ParabolicReflectorApp < matlab.apps.AppBase

    
    properties (Access = public)
        UIFigure              matlab.ui.Figure
        MainGridLayout        matlab.ui.container.GridLayout
        ControlPanel          matlab.ui.container.Panel
        VisualTabGroup        matlab.ui.container.TabGroup
        MetricsPanel          matlab.ui.container.Panel
        GeoTab                matlab.ui.container.Tab
        ApertureTab           matlab.ui.container.Tab
        EffTab                matlab.ui.container.Tab
        TwoDTab               matlab.ui.container.Tab 
        ThreeDTab             matlab.ui.container.Tab
        GeoAxes               matlab.ui.control.UIAxes
        ApAxes                matlab.ui.control.UIAxes
        EffAxes               matlab.ui.control.UIAxes
        TwoDAxes              matlab.ui.control.UIAxes
        ThreeDAxes            matlab.ui.control.UIAxes
        FreqEdit              matlab.ui.control.NumericEditField
        DiaEdit               matlab.ui.control.NumericEditField
        FDRatioEdit           matlab.ui.control.NumericEditField
        Theta0Edit            matlab.ui.control.NumericEditField
        ConfigDrop            matlab.ui.control.DropDown 
        FeedModelDrop         matlab.ui.control.DropDown 
        HornQEdit             matlab.ui.control.NumericEditField
        SubRefDiaEdit         matlab.ui.control.NumericEditField
        DirLabel              matlab.ui.control.Label
        EffLabel              matlab.ui.control.Label
        SpillLabel            matlab.ui.control.Label
        TaperLabel            matlab.ui.control.Label
        BlockLabel            matlab.ui.control.Label
        WarningLabel          matlab.ui.control.Label
        IsUpdating            logical = false;
    end
    
    methods (Access = private)
        
        function UpdateAnalysis(app, source)
            if app.IsUpdating, return; end
            app.IsUpdating = true;
            
            try
                %INPUTS
                f_ghz = app.FreqEdit.Value;
                D = app.DiaEdit.Value;
                config = app.ConfigDrop.Value;   
                feedModel = app.FeedModelDrop.Value; 
                q_pow = app.HornQEdit.Value;     
                d_sub = app.SubRefDiaEdit.Value; 
                lambda = 0.3 / f_ghz;
                
                % SYNC (this code i used to ensure the sync between diff parameters)
                if strcmp(source, 'FD')
                    fd = app.FDRatioEdit.Value;
                    theta0_rad = 2 * acot(4 * fd);
                    app.Theta0Edit.Value = rad2deg(theta0_rad);
                elseif strcmp(source, 'Theta')
                    theta0_deg = app.Theta0Edit.Value;
                    theta0_rad = deg2rad(theta0_deg);
                    fd = 0.25 * cot(theta0_rad/2);
                    app.FDRatioEdit.Value = fd;
                else
                    fd = app.FDRatioEdit.Value;
                    theta0_rad = 2 * acot(4 * fd);
                end
                
                theta0_deg = rad2deg(theta0_rad); 
                f = fd * D; 
                
                % BLOCKAGE CHECK
                if strcmp(config, 'Cassegrain')
                    app.SubRefDiaEdit.Enable = 'on';
                    D_block = d_sub;
                    app.WarningLabel.Text = 'Cassegrain Mode';
                    app.WarningLabel.FontColor = [0 0.5 0];
                else
                    app.SubRefDiaEdit.Enable = 'off';
                    D_block = d_sub; 
                    app.WarningLabel.Text = 'Front-Fed Mode';
                    app.WarningLabel.FontColor = [0 0.5 0];
                end
                
                if D_block >= D
                     error('Blockage > Dish Diameter');
                end
                
                %  EFFICIENCY CALCS
                th_max = pi/2; 
                th = linspace(0, th_max, 500);
                
                if strcmp(feedModel, 'Dipole (Eq.16)')
                    G_f = (cos(th)).^2; 
                    q_used = 2;
                else
                    G_f = (cos(th)).^q_pow;
                    q_used = q_pow;
                end
                
                % Integrals
                idx_cap = th <= theta0_rad;
                P_cap = trapz(th(idx_cap), G_f(idx_cap) .* sin(th(idx_cap)));
                P_tot = trapz(th, G_f .* sin(th));
                eff_spill = P_cap / P_tot;
                
                integrand = sqrt(G_f(idx_cap)) .* tan(th(idx_cap)/2);
                integral_val = trapz(th(idx_cap), integrand);
                eta_illum = 2 * (cot(theta0_rad/2)^2) * (abs(integral_val)^2) / P_tot;
                eff_taper = eta_illum / eff_spill;
                if isnan(eff_taper), eff_taper = 0; end
                
                eff_block = (1 - (D_block/D)^2)^2;
                eta_total = eta_illum * eff_block;
                D_iso = (pi * D / lambda)^2 * eta_total;
                D_db = 10*log10(D_iso);
                
                % Update Metrics
                app.DirLabel.Text = sprintf('Directivity: %.2f dBi', D_db);
                app.EffLabel.Text = sprintf('Total Eff: %.1f%%', eta_total*100);
                app.SpillLabel.Text = sprintf('Spillover: %.1f%%', eff_spill*100);
                app.TaperLabel.Text = sprintf('Taper: %.1f%%', eff_taper*100);
                app.BlockLabel.Text = sprintf('Blockage: %.1f%%', eff_block*100);
                
                % ============== visuals =================
                
                % RAY TRACING
                cla(app.GeoAxes); hold(app.GeoAxes, 'on');
                y_main = linspace(-D/2, D/2, 200);
                x_main = (y_main.^2) / (4*f);
                plot(app.GeoAxes, x_main, y_main, 'k', 'LineWidth', 2.5);
                
                num_rays = 11; 
                ray_colors_in = [0.8 0 0];   
                ray_colors_out = [0 0.4 0.8]; 
                
                if strcmp(config, 'Cassegrain')
                    z_sub = f*0.7; 
                    plot(app.GeoAxes, [z_sub z_sub], [-d_sub/2 d_sub/2], 'm', 'LineWidth', 3);
                    plot(app.GeoAxes, 0, 0, 'rx', 'MarkerSize', 8, 'LineWidth', 2);
                    y_rays = linspace(0.01, d_sub/2, num_rays); 
                    for y_r = y_rays
                        scale = (D/2) / (d_sub/2); y_m = y_r * scale; x_m = (y_m^2)/(4*f);
                        plot(app.GeoAxes, [0 z_sub], [0 y_r], 'Color', ray_colors_in);
                        plot(app.GeoAxes, [0 z_sub], [0 -y_r], 'Color', ray_colors_in);
                        plot(app.GeoAxes, [z_sub x_m], [y_r y_m], 'Color', ray_colors_in);
                        plot(app.GeoAxes, [z_sub x_m], [-y_r -y_m], 'Color', ray_colors_in);
                        plot(app.GeoAxes, [x_m f*1.6], [y_m y_m], 'Color', ray_colors_out);
                        plot(app.GeoAxes, [x_m f*1.6], [-y_m -y_m], 'Color', ray_colors_out);
                    end
                else
                    plot(app.GeoAxes, f, 0, 'rx', 'MarkerSize', 10);
                    if d_sub>0, plot(app.GeoAxes, [f f], [-d_sub/2 d_sub/2], 'k-', 'LineWidth', 3); end
                    y_rays = linspace(d_sub/2 + 0.05, D/2, num_rays);
                    for y_r = y_rays
                        x_r = (y_r^2)/(4*f);
                        plot(app.GeoAxes, [f x_r], [0 y_r], 'Color', ray_colors_in);
                        plot(app.GeoAxes, [f x_r], [0 -y_r], 'Color', ray_colors_in);
                        plot(app.GeoAxes, [x_r f*1.6], [y_r y_r], 'Color', ray_colors_out);
                        plot(app.GeoAxes, [x_r f*1.6], [-y_r -y_r], 'Color', ray_colors_out);
                    end
                end
                axis(app.GeoAxes, 'equal'); grid(app.GeoAxes, 'on');
                xlabel(app.GeoAxes, 'Z (m)'); ylabel(app.GeoAxes, 'Y (m)');
                
                % CROSS POLARIZATION 
                cla(app.ApAxes); hold(app.ApAxes, 'on');
                [rho, phi] = meshgrid(linspace(0, D/2, 16), linspace(0, 2*pi, 24));
                [Y, Z] = pol2cart(phi, rho);
                theta_grid = 2 * atan(rho ./ (2*f));
                
                if strcmp(feedModel, 'Horn (Ideal)')
                    Ex = zeros(size(Y)); 
                    Ey = -ones(size(Y)) .* (cos(th_max * rho/(D/2)).^q_pow); 
                    t_str = 'Ideal Horn (Huygens Source)'; col = 'b';
                else
                    term1 = sin(phi).*cos(phi).*(1 - cos(theta_grid));
                    term2 = sin(phi).^2 .* cos(theta_grid) + cos(phi).^2;
                    mag = cos(theta_grid).^2; 
                    Ex = mag .* term1;
                    Ey = -mag .* term2;
                    t_str = 'Dipole Feed (Eq. 16)'; col = 'r';
                end
                quiver(app.ApAxes, Y, Z, Ex, Ey, 0.6, 'Color', col, 'LineWidth', 1.5);
                title(app.ApAxes, t_str);
                axis(app.ApAxes, 'equal'); xlim(app.ApAxes, [-D/2 D/2]); ylim(app.ApAxes, [-D/2 D/2]);
                
                % SPILLOVER
                cla(app.EffAxes); hold(app.EffAxes, 'on');
                yyaxis(app.EffAxes, 'left');
                plot(app.EffAxes, rad2deg(th), G_f, 'LineWidth', 2);
                xline(app.EffAxes, theta0_deg, 'r--', 'Dish Edge');
                fill(app.EffAxes, [theta0_deg, 90, 90, theta0_deg], [0, 0, 1, 1], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                title(app.EffAxes, sprintf('Feed Illumination (n=%.1f)', q_used));
                xlim(app.EffAxes, [0 90]); grid(app.EffAxes, 'on');
                
                % 2D PATTERN 
                cla(app.TwoDAxes); hold(app.TwoDAxes, 'on');
                theta_2d = linspace(-10, 10, 500); % +/- 10 degrees zoom
                u_2d = (pi*D/lambda) * sind(theta_2d);
                % J1(u)/u
                val = ones(size(u_2d));
                mask = abs(u_2d) > 1e-6;
                val(mask) = 2*besselj(1,u_2d(mask))./u_2d(mask);
                Pat_2d_dB = 20*log10(abs(val)) + D_db; % Normalize to Gain
                
                plot(app.TwoDAxes, theta_2d, Pat_2d_dB, 'LineWidth', 2, 'Color', [0 0.45 0.74]);
                yline(app.TwoDAxes, D_db - 3, 'k--', 'HPBW Level');
                yline(app.TwoDAxes, D_db - 17.6, 'r:', 'First Sidelobe (-17.6dB)');
                
                ylim(app.TwoDAxes, [D_db-40, D_db+2]);
                xlabel(app.TwoDAxes, 'Theta (deg)');
                ylabel(app.TwoDAxes, 'Gain (dBi)');
                title(app.TwoDAxes, '2D Cartesian Cut (Zoomed)');
                grid(app.TwoDAxes, 'on');
                
                % 3D PATTERN 
                cla(app.ThreeDAxes);
                

                theta_p = linspace(0, deg2rad(15), 100); 
                phi_p = linspace(0, 2*pi, 80);
                [TH, PH] = meshgrid(theta_p, phi_p);
                
                u_3d = (pi*D/lambda) * sin(TH);
                val = ones(size(u_3d));
                mask = u_3d > 1e-6;
                val(mask) = 2*besselj(1,u_3d(mask))./u_3d(mask);
                Pat_Lin = val.^2;
                Pat_dB = 10*log10(Pat_Lin);
                

                % Floor at -40dB, Range 40dB
                Range_dB = 40;
                R_vis = Pat_dB + Range_dB; 
                R_vis(R_vis < 0) = 0;
                
                X = R_vis .* sin(TH) .* cos(PH);
                Y = R_vis .* sin(TH) .* sin(PH);
                Z = R_vis .* cos(TH);
                
                s = surf(app.ThreeDAxes, X, Y, Z, Pat_dB);
                s.EdgeColor = 'none';
                s.FaceColor = 'interp';
                
                view(app.ThreeDAxes, 3);
                axis(app.ThreeDAxes, 'equal');
                camlight(app.ThreeDAxes, 'headlight');
                lighting(app.ThreeDAxes, 'gouraud');
                material(app.ThreeDAxes, 'dull');
                colormap(app.ThreeDAxes, 'jet');
                caxis(app.ThreeDAxes, [-30 0]);
                
                title(app.ThreeDAxes, '3D Pattern (Radius = dB Magnitude)');
                xlabel(app.ThreeDAxes, 'U (Relative)'); ylabel(app.ThreeDAxes, 'V (Relative)');
                
            catch ME
                app.WarningLabel.Text = 'Error Calculation';
                disp(ME.message);
            end
            app.IsUpdating = false;
        end
        
        % UI CONSTRUCTION
        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Parabolic Reflector Designer', 'WindowState', 'maximized'); 
            app.MainGridLayout = uigridlayout(app.UIFigure, [1 2]);
            app.MainGridLayout.ColumnWidth = {320, '1x'};
            
            app.ControlPanel = uipanel(app.MainGridLayout, 'Title', 'Controls');
            cgrid = uigridlayout(app.ControlPanel, [14 2]);
            cgrid.RowHeight = {'fit','fit','fit','fit','fit','fit','fit','fit','fit','fit','fit','fit','1x','fit'};

            uilabel(cgrid, 'Text', 'Freq (GHz):');
            app.FreqEdit = uieditfield(cgrid, 'numeric', 'Value', 12, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('Freq'));
            
            uilabel(cgrid, 'Text', 'Config:');
            app.ConfigDrop = uidropdown(cgrid, 'Items', {'Front-Fed', 'Cassegrain'}, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('Config'));
            
            uilabel(cgrid, 'Text', 'Dish Dia (m):');
            app.DiaEdit = uieditfield(cgrid, 'numeric', 'Value', 1.5, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('Dia'));
            
            uilabel(cgrid, 'Text', 'Sub-Ref Dia (m):');
            app.SubRefDiaEdit = uieditfield(cgrid, 'numeric', 'Value', 0.2, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('Blk'));
            
            uilabel(cgrid, 'Text', 'f/D Ratio:');
            app.FDRatioEdit = uieditfield(cgrid, 'numeric', 'Value', 0.6, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('FD'));
            
            uilabel(cgrid, 'Text', 'Theta0 (deg):');
            app.Theta0Edit = uieditfield(cgrid, 'numeric', 'Value', 45.2, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('Theta'));
            
            uilabel(cgrid, 'Text', 'Feed Model:');
            app.FeedModelDrop = uidropdown(cgrid, 'Items', {'Horn (Ideal)', 'Dipole (Eq.16)'}, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('Feed'));
            
            uilabel(cgrid, 'Text', 'Horn Power (q):');
            app.HornQEdit = uieditfield(cgrid, 'numeric', 'Value', 4, 'ValueChangedFcn', @(s,e)app.UpdateAnalysis('Q'));
            
            app.MetricsPanel = uipanel(cgrid);
            app.MetricsPanel.Layout.Row = 13; app.MetricsPanel.Layout.Column = [1 2];
            app.MetricsPanel.Title = 'Design Metrics';
            mgrid = uigridlayout(app.MetricsPanel, [6 1]);
            mgrid.RowHeight = {'fit','fit','fit','fit','fit','fit'};
            
            app.DirLabel = uilabel(mgrid, 'Text', 'Directivity: --', 'FontWeight', 'bold', 'FontColor', 'blue');
            app.EffLabel = uilabel(mgrid, 'Text', 'Total Eff: --');
            app.SpillLabel = uilabel(mgrid, 'Text', 'Spillover: --');
            app.TaperLabel = uilabel(mgrid, 'Text', 'Taper: --');
            app.BlockLabel = uilabel(mgrid, 'Text', 'Blockage: --');
            app.WarningLabel = uilabel(mgrid, 'Text', 'Status: OK', 'FontColor', [0 0.6 0]);
            
            app.VisualTabGroup = uitabgroup(app.MainGridLayout);
            
            app.GeoTab = uitab(app.VisualTabGroup, 'Title', 'Ray Tracing');
            g1 = uigridlayout(app.GeoTab); g1.ColumnWidth = {'1x'}; g1.RowHeight = {'1x'};
            app.GeoAxes = uiaxes(g1);
            
            app.ApertureTab = uitab(app.VisualTabGroup, 'Title', 'Polarization');
            g2 = uigridlayout(app.ApertureTab); g2.ColumnWidth = {'1x'}; g2.RowHeight = {'1x'};
            app.ApAxes = uiaxes(g2);
            
            app.EffTab = uitab(app.VisualTabGroup, 'Title', 'Efficiency');
            g3 = uigridlayout(app.EffTab); g3.ColumnWidth = {'1x'}; g3.RowHeight = {'1x'};
            app.EffAxes = uiaxes(g3);
            
            app.TwoDTab = uitab(app.VisualTabGroup, 'Title', '2D Cut (Cartesian)');
            g4 = uigridlayout(app.TwoDTab); g4.ColumnWidth = {'1x'}; g4.RowHeight = {'1x'};
            app.TwoDAxes = uiaxes(g4);
            

            app.ThreeDTab = uitab(app.VisualTabGroup, 'Title', '3D Pattern (dB)');
            g5 = uigridlayout(app.ThreeDTab); g5.ColumnWidth = {'1x'}; g5.RowHeight = {'1x'};
            app.ThreeDAxes = uiaxes(g5);
            
            app.UIFigure.Visible = 'on';
            app.UpdateAnalysis('Init');
        end
    end
    
    methods (Access = public)
        function app = ParabolicReflectorApp
            createComponents(app)
            registerApp(app, app.UIFigure)
        end
        function delete(app)
            delete(app.UIFigure)
        end
    end
end