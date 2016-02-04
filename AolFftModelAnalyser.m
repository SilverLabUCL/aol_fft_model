classdef AolFftModelAnalyser < handle

    properties
        waves
        time_range
        z_range
        wavelengths
        wavelength_weightings
        plot_intensity_squared
        colormap
    end
    
    methods
        function obj = AolFftModelAnalyser(waves)
            obj.waves = waves;
            obj.time_range = 0;
            obj.z_range = 0;
            obj.wavelengths = [800, 800-1.75, 800+1.75, 800-2.5, 800+2.5] * 1e-9; 
            % obj.wavelengths = [920, 920-2.5, 920+2.5, 920-3.5, 920+3.5] * 1e-9;
            obj.wavelength_weightings = [1, 1/sqrt(2), 1/sqrt(2), 1/2, 1/2];
            obj.plot_intensity_squared = false;
            obj.colormap = 'jet';
        end
        
        function res = analyse_and_plot_focus(obj)
            intensity_3d = calculate_psf_through_aol(obj);
            res = obj.get_psf_dimensions(intensity_3d);
            obj.plot_psf_xy(propagated_wave_2d)
            obj.plot_psf_xz(propagated_wave_2d)
        end
        
        function aol_fft_model = get_aol_fft_model(obj)
            aol_fft_model = AolFftModel();
            aol_fft_model.z_range = obj.z_range;
        end
        
        function intensity_sum = calculate_psf_through_aol(obj)
            aol_fft_model = obj.get_aol_fft_model();
            intensity_sum = 0;
            
            for n = 1:numel(obj.wavelengths)    
            for time = obj.time_range 
                aol_fft_model.wavevector = 2*pi/obj.wavelengths(n);
                aol_fft_model.time = time;
                focal_region_wave_3d = aol_fft_model.calculate_psf_3d(obj.waves);
                intensity_sum = intensity_sum + abs(focal_region_wave_3d).^2 * obj.wavelength_weightings(n);
            end
            end
        end

        function res = get_psf_dimensions(~, field_3d)
            m = get_aol_fft_model();
            [~, x, ~] = m.get_focal_region_xy();
            z = m.z_range;
            
            max_intensity_sqr = max(abs(field_3d(:)).^2);

            half_or_more_x = squeeze(max(abs(field_3d), [], 3)).^2 >= max_intensity_sqr/2;
            x_res = max(x(half_or_more_x)) - min(x(half_or_more_x));
            half_or_more_z = squeeze(max(max(abs(field_3d)))).^2 >= max_intensity_sqr/2;
            z_res = max(z(half_or_more_z)) - min(z(half_or_more_z));
            
            max_val = max(abs(field_3d(:)));
            x_pos = median(x(max(abs(field_3d), [], 3) == max_val));
            z_pos = median(z(max(max(abs(field_3d), [], 1), [], 2) == max_val));

            max_intensity_sqr = max(max(abs(field_3d(:,:,ceil(numel(z)/2))).^2));
            total_fluores = sum(abs(field_3d(:)).^2);
            res = [[x_pos, z_pos] * 1e6, [x_res, z_res] * 1e6, log10(max_intensity_sqr), log10(total_fluores)];
        end
        
        function plot_psf_xy(obj, propagated_wave_2d)
            [~, x, y] = get_aol_fft_model().aol_region_xy();
            [~,z_plane_index] = max(squeeze(max(max(propagated_wave_2d))));
            obj.plot_psf(x, y, abs(propagated_wave_2d(:,:,z_plane_index)).^obj.pwr);
        end          
        function plot_psf_xz(obj, propagated_wave_2d)
            [~, x, ~] = get_aol_fft_model().aol_region_xy();
            [zz, xx] = meshgrid(obj.z_list, max(x,[],1));
            obj.plot_psf(xx, zz, abs(squeeze(max(propagated_wave_2d.^obj.pwr, [], 1))));
        end            
        function plot_psf_yz(obj, propagated_wave_2d)
            [~, ~, y] = get_aol_fft_model().aol_region_xy();
            [zz, yy] = meshgrid(obj.z_list, max(y,[],2));
            obj.plot_psf(yy, zz, abs(squeeze(max(propagated_wave_2d.^obj.pwr, [], 2))));
        end 
        function plot_psf(obj, a, b, c)
            h = pcolor(a, b, c);
            set(h,'EdgeColor','none')
            colormap(obj.colormap)
            %caxis([0,2000^obj.pwr])
            axis equal
            %axis tight
            %axis off
            %set(gcf, 'Position', [0,0,800,800]);
            %set(gca,'position',[0 0 1 1],'units','normalized')
        end
        
        function plot_focal_plane(obj, waves)
            focal_plane_wave_2d = obj.get_aol_fft_model().get_wave_in_focal_plane(waves);
            labels = {'focal plane', 'x', 'y'};
            [~, x, y] = get_aol_fft_model().focal_region_xy();
            obj.plot_wave_plane(focal_plane_wave_2d, x, y, labels)
        end
        function plot_leaving_aol(obj, waves)
            sampled_wave = obj.get_aol_fft_model().get_wave_leaving_aol(waves);
            figure; hold on; 
            plot(unwrap(ifftshift(angle(sampled_wave(ceil(size(sampled_wave,2)/2),:)))), 'k'); 
            plot(unwrap(ifftshift(angle(sampled_wave(:,ceil(size(sampled_wave,2)/2))))), 'r--');
            labels = {'input to relay and obj', 'x', 'y'};
            [~, x, y] = get_aol_fft_model().aol_region_xy();
            obj.plot_wave_plane(sampled_wave, x, y, labels)
        end
        function plot_aod_phase_shift(obj, n, waves)
            directions = linspace(0, 2 * pi, num_aods + 1);
            sampled_wave = obj.get_aol_fft_model().aod_phase_shift(directions(n), waves{n});
            figure; hold on; 
            plot(unwrap(ifftshift(angle(sampled_wave(ceil(size(sampled_wave,2)/2),:)))), 'k'); 
            plot(unwrap(ifftshift(angle(sampled_wave(:,ceil(size(sampled_wave,2)/2))))), 'r--');
            labels = {'aod', 'x', 'y'};
            [~, x, y] = get_aol_fft_model().aol_region_xy();
            obj.plot_wave_plane(sampled_wave, x, y, labels)
        end 
        function plot_wave_plane(obj, wave_function_2d, x, y, labels)
            figure()
            obj.subplot_wave_plane(abs(wave_function_2d), x, y, labels, 1)
            obj.subplot_wave_plane(angle(wave_function_2d), x, y, labels, 2)
        end
        function subplot_wave_plane(~, z, x, y, labels, n)
            subplot(1,2,n)
            hh = pcolor(x,y,z);
            title(labels{1})
            xlabel(labels{2})
            ylabel(labels{3})
            set(hh,'EdgeColor','none')
        end
    end
end

