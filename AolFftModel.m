classdef AolFftModel < handle
    
    properties
        fft_adjustment
        fft_number_of_samples

        wavevector
        acoustic_velocity
        aod_spacing
        aod_half_aperture_width
        beam_width
        aol_to_objective_scaling
        objective_magnification
        tube_lens_focal_length

        time
        z_range
    end
    
    methods
        function obj = AolFftModel()
            obj.fft_adjustment = 1e2; % adjustments as necessary for the FFT - sets ratio of k to max(kx)
            obj.fft_number_of_samples = 2.^8 - 1; % computational speed vs accuracy
            
            obj.wavevector = 2*pi/800e-9;
            obj.acoustic_velocity = 613; % speed of sound in TeO2
            obj.aod_spacing = 0.04554;
            obj.aod_half_aperture_width = 15e-3; % the aperture width
            obj.beam_width = 3.75e-3; % 68% of field within +- beam_sigma, 95% of field within +- 2 beam_sigma, equiv. beam intensity falls off to 1/e after beam_sigma
            obj.aol_to_objective_scaling = 0.8; % scaling by the relay between AOL and objective
            obj.objective_magnification = 20;
            obj.tube_lens_focal_length = 160e-3;
            
            obj.z_range = linspace(-100,100,199)*1e-6; % distances from the nominal focal plane of the objective to evaluate the field at             
            obj.time = 0;
        end
        
        function focal_region_wave_3d = calculate_psf_3d(obj, waves)
            focal_plane_wave_2d = obj.get_wave_in_focal_plane(waves);
            space_width_focal_plane = obj.focal_region_xy();
            ref_ind_water = 4/3;
            focal_region_wave_3d = obj.propagate_wave(focal_plane_wave_2d, space_width_focal_plane, ref_ind_water);
        end
        
        function focal_plane_wave_2d = get_wave_in_focal_plane(obj, waves)
            wave_leaving_aol = obj.get_wave_leaving_aol(waves);
            focal_plane_wave_2d = fftshift(fft2(ifftshift(wave_leaving_aol))); % transform to focal plane
        end
        
        function sampled_wave = get_wave_leaving_aol(obj, waves)     
            % take a number of waves of phase shift for each AOD and a time to calculate the wave out of the last AOD
            num_aods = numel(waves.aods);
            distances = [ones(1,num_aods-1)*obj.aod_spacing, 0];
            directions = linspace(0, 2 * pi, num_aods + 1);
            ref_ind_air = 1; 
            
            sampled_wave = obj.input_wave();
            
            for n = 1:num_aods
                phase_shift = obj.aod_phase_shift(directions(n), waves{n});
                sampled_wave = obj.propagate_wave(sampled_wave, distances(n), space_width, ref_ind_air) .* phase_shift;
            end
            
            sampled_wave = sampled_wave .* obj.aperture(obj.aod_half_aperture_width ./ obj.aol_to_objective_scaling);
            sampled_wave = sampled_wave .* exp(1i * 2*pi * waves.focus     * ((x.^2 + y.^2)    ./ obj.aod_half_aperture_width.^2));
            sampled_wave = sampled_wave .* exp(1i * 2*pi * waves.spherical * ((x.^2 + y.^2).^2 ./ obj.aod_half_aperture_width.^4)); 
        end
        
        function input = input_wave(obj)
            [x, y, ~] = obj.aol_region_xy();
            gaussian = exp(-(x.^2 + y.^2)/2/(obj.beam_width).^2);
            input = gaussian .* obj.aperture(obj.aod_half_aperture_width);
        end
        
        function apertured = aperture(obj, half_width)
            [x, y, ~] = obj.aol_region_xy();
            apertured = (sqrt(x.^2 + y.^2) < half_width);
        end

        function propagated_wave = propagate_wave(obj, wave_2d, distance, space_width, ref_ind)
            % use FFT to propagate the field backwards and forwards
            count = [0:floor((obj.fft_number_of_samples - 1)/2) floor(-(obj.fft_number_of_samples - 1)/2):-1];
            
            [k_x, k_y] = meshgrid(count * 2 * pi / space_width);
            k_z = sqrt(obj.wavevector.^2 * ref_ind.^2 - k_x.^2 - k_y.^2);
            mask = ~(abs(imag(k_z)) > 0);

            ft_wave_2d = fft2(ifftshift(wave_2d)); % recentre wave with fftshift
            
            propagated_wave = zeros([size(k_z), numel(distance)]);
            for n = 1:numel(distance) % to handle multiple distances in one hit
                propagated_wave(:,:,n) = fftshift(ifft2((ft_wave_2d .* exp(1i * mask .* k_z .* distance(n)) .* mask))); % shift back 
            end
            propagated_wave = squeeze(propagated_wave); % in case only one distance required
        end
        
        function phase_shift = aod_phase_shift(obj, orientation, aod_waves)
            r = x.*cos(orientation) + y.*sin(orientation) - obj.acoustic_velocity * t;
            aod_phase = 2*pi * aod_waves ./ obj.aod_half_aperture_width .^ (1:5);
            phase_shift = 0;
            for m = 1:5
                phase_shift = phase_shift * exp(1i * (aod_phase(m) .* r.^m));
            end
        end
        
        function [width, x, y] = focal_region_xy(obj)
            space_width = pi * obj.fft_number_of_samples / obj.wavevector * obj.fft_adjustment;
            integer_list = floor(-obj.fft_number_of_samples/2+0.5):floor(obj.fft_number_of_samples/2-0.5);
            k_x_discrete = integer_list * 2 * pi / (space_width * obj.aol_to_objective_scaling);
            plane_wave_angle_air = k_x_discrete / obj.wavevector;
            [x, y] = meshgrid(plane_wave_angle_air * obj.tube_lens_focal_length / obj.objective_magnification); % new positions in the focal plane  
            width = (max(x(:)) - min(x(:)));
            %width = obj.fft_number_of_samples * 2 * pi / (space_width * obj.aol_to_objective_scaling) / obj.wavevector * obj.tube_lens_focal_length / obj.objective_magnification;
        end        
        
        function [width, x, y] = aol_region_xy(obj)
            samples = linspace(-1/2, 1/2, obj.fft_number_of_samples) * space_width;
            [x, y] = meshgrid(samples);
            width = pi * obj.fft_number_of_samples / obj.wavevector * obj.fft_adjustment;
        end
    end
end

