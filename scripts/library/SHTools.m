%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spherical Harmonics Synthesis and Analysis Tools (SHSAT)
% - Spherical coordinate Azimuth convention: 0 - 2pi, anti-clockwise from x = +1
% - Spherical coordinate Elevation convention: 0 - pi, from z = +1 to -1
% (c) Abdullah Fahim, abdullah.fahim@hotmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef SHTools < handle
    methods (Static)
        %% Spherical harmonics
        function y = getRealSHPerMode(n, m, el, az)
            % y = getRealSHPerMode(n, m, az, el)
            % n, m = order and mode, scalar
            % az, el = Azimuth and elevation in radian, vector Qx1
            % y = Real Spherical Harmonics, Qx1

            L = legendre(n, cos(el));
            L = L(abs(m) + 1, :);

            Norm = (-1).^m * sqrt( (2*n+1)*factorial(n-abs(m)) ./ (4*pi*factorial(n+abs(m))) );

            if m < 0
                K = sqrt(2)*sin(abs(m) * az);
            elseif m > 0
                K = sqrt(2)*cos(m * az);
            else
                K = 1;
            end


            y = Norm .* L(:) .* K(:);
        end
        
        function y = getRealSH(order, el, az)
            % y = getRealSH(order, el, az)
            % order = SH order
            % az, el = Azimuth and elevation in radian, vector Qx1
            % y = Real Spherical Harmonics, (order+1)^2 x Q
            
            
            el = el(:); % because Legndre always return in column format for each theta
            az = az(:);
            
            if numel(el) ~= numel(az)
                if numel(el) == 1
                    el = repmat(el, numel(az), 1);
                elseif numel(az) == 1
                    az = repmat(az, numel(el), 1);
                else
                    error('Number mismatch between azimuths and elevations.');
                end
            end
            
            NN = (order + 1)^2;
            Q = length(el);
            
            y = zeros(NN, Q);

            y(1, :) = 1 / sqrt(4*pi);
            
            if order > 0
                ind = 2;
                for n = 1:order
                    % (n+1) x Q
                    L = legendre(n, cos(el));
                    L = [flipud(L(2:end, :)); L];
                    
                    m=(0:n)';
                    
                    K = zeros(2*n+1, Q);
                    
                    % m = 0
                    K(n+1,:) = 1;
                    
                    % Positive m startes from n+1+1th index onwards
                    K(n + 1 + m(2:end),:) = sqrt(2) * cos(m(2:end) * az');
                    
                    % Negative m startes from n+1-1th index backwards
                    K(n + 1 - m(end:-1:2),:) = sqrt(2)*sin(m(end:-1:2) * az');
                    
                    m = (-n:n)';
                    Norm = (-1).^m .* sqrt( (2*n+1)*factorial(n-abs(m)) ./ (4*pi*factorial(n+abs(m))) );
                    
                    nextind = (n+1)^2;
                    % ( 1 - 2 * (m < 0 & rem(m,2) ~=0) )
                    y(ind:nextind, :) = bsxfun(@times, Norm, L.*K );
                    
                    ind = nextind+1;
                end
            end
        end

        function y = getComplexSHPerMode(n, m, el, az)
            % y = getComplexSHPerMode(n, m, el, az)
            % n = order, scaler
            % m = mode (-n <= m <= n), scaler
            % az, el = Azimuth and elevation in radian, vector Qx1
            % y = Complex Spherical Harmonics, Qx1
            % Wolfram Definition, http://keisan.casio.com/exec/system/1180573409

            el = el(:); % because Legndre always return in column format for each theta
            az = az(:);

            M = abs(m);

            if (M > n)
                y = 0;
                return;
            end

            m_index = M + 1;
            a1=((2*n+1)/(4*pi));

            a2=factorial(n-M) ./ factorial(n+M);
            C=sqrt(a1*a2);
            L=legendre(n,cos(el));
            % transpose L to convert L(m_index, :) to column
            % Y_{n |m|}(.)
            y=C .* transpose(L(m_index, :)).* exp(1i .* M .* az);
            if(m < 0)
                y = (-1) .^ M .* conj(y);
            end
        end
        
        function y = getComplexSH(order, el, az)
            % y = getComplexSH(order, el, az)
            % order = Highest order to calculate
            % az, el = Azimuth and elevation in radian, vector Qx1
            % y = Complex Spherical Harmonics, (order+1)^2 x Q
            % Wolfram Definition, http://keisan.casio.com/exec/system/1180573409

            el = el(:); % because Legndre always return in column format for each theta
            az = az(:);
            
            if numel(el) ~= numel(az)
                if numel(el) == 1
                    el = repmat(el, numel(az), 1);
                elseif numel(az) == 1
                    az = repmat(az, numel(el), 1);
                else
                    error('Number mismatch between azimuths and elevations.');
                end
            end
            
            NN = (order + 1)^2;
            Q = length(el);
            
            y = zeros(NN, Q);

            y(1, :) = 1 / sqrt(4*pi);
            
            if order > 0
                ind = 2;
                for n = 1:order
                    % (n+1) x Q
                    L = legendre(n, cos(el));
                    L = [flipud(L(2:end, :)); L];
                    
                    m = (-n:n)';
                    K = exp(1i * m * az');
                    
                    M = abs(m);
                    
                    Norm = (2*n+1) / (4*pi) * factorial(n-M) ./ factorial(n+M);
                    
                    nextind = (n+1)^2;
                    % ( 1 - 2 * (m < 0 & rem(m,2) ~=0) )
                    y(ind:nextind, :) = (-1).^( (m < 0) .* m ) .* bsxfun(@times, sqrt(Norm), L.*K );
                    
                    ind = nextind+1;
                end
            end
        end
        
        function y = getComplexEvenSH(order, el, az)
            % y = getComplexEvenSH(order, el, az)
            % order = Highest order to calculate
            % az, el = Azimuth and elevation in radian, vector Qx1
            % y = Even Complex Spherical Harmonics, (order + 1) * (order + 2) / 2 x Q
            % Wolfram Definition, http://keisan.casio.com/exec/system/1180573409

            el = el(:); % because Legndre always return in column format for each theta
            az = az(:);
            
            if numel(el) ~= numel(az)
                if numel(el) == 1
                    el = repmat(el, numel(az), 1);
                elseif numel(az) == 1
                    az = repmat(az, numel(el), 1);
                else
                    error('Number mismatch between azimuths and elevations.');
                end
            end
            
            NN = (order + 1) * (order + 2) / 2;
            Q = length(el);
            
            y = zeros(NN, Q);

            y(1, :) = 1 / sqrt(4*pi);
            
            if order > 0
                ind = 2;
                for n = 1:order
                    % (n+1) x Q
                    L = legendre(n, cos(el));
                    L = [flipud(L(2:end, :)); L];
                    L = L(1:2:end, :);
                    
                    m = (-n:2:n)';
                    K = exp(1i * m * az');
                    
                    M = abs(m);
                    
                    Norm = (2*n+1) / (4*pi) * factorial(n-M) ./ factorial(n+M);
                    
                    nextind = ind + length(Norm) - 1;
                    % ( 1 - 2 * (m < 0 & rem(m,2) ~=0) )
                    y(ind:nextind, :) = (-1).^( (m < 0) .* m ) .* bsxfun(@times, sqrt(Norm), L.*K );
                    
                    ind = nextind+1;
                end
            end
        end
        
        %% Spherical harmonics decomposition
        function Alpha = alphaFarField3DPerMode(n, m, el_i, az_i, amp_i, polarity)
            % Alpha = alphaFarField3DPerMode(n, m, el_i, az_i, amp_i, polarity)
            % n = order, scaler
            % m = mode, scaler
            % el_i => incident elevations, K x 1
            % az_i => incident azimuths, K x 1
            % amp_i = incident amplitudes, scalar or K(source) x T (time snapshot)
            % polarity = 'sink' or 'source'. Eigenmik uses sink
            % Alpha = Spherical harmonics coefficients for plane wave, 1 x T



            if(strcmp(polarity, 'sink'))
                const_i = 1i;
            else
                const_i = -1i;
            end

            if(abs(m) > n)
                error('Mode can not be greater than order');
            end


            el_i = el_i(:);
            az_i = az_i(:);

            Alpha = 4 * pi * const_i.^n * ctranspose(SHTools.getComplexSHPerMode(n, m, el_i, az_i)) * amp_i;

            % We used source at infinity, exp(-1i .* k' * x) => alpha = (-1i)^n
            % For sink at infinity, we should use exp(1i .* k' * x) => alpha = (1i)^n
        end
        
        function Alpha = alphaFarField3D(order, el_i, az_i, amp_i, polarity)
            % Alpha = alphaFarField3D(order, el_i, az_i, amp_i, polarity)
            % order = scaler
            % el_i => incident elevations, K x 1
            % az_i => incident azimuths, K x 1
            % amp_i = incident amplitudes, scalar or K(source) x T (time snapshot)
            % polarity = 'sink' or 'source'. Practical mics use sink
            % Alpha = Spherical harmonics coefficients for plane wave, (order+1)^2 x T



            if(strcmp(polarity, 'sink'))
                const_i = 1i;
            else
                const_i = -1i;
            end

            el_i = el_i(:);
            az_i = az_i(:);
            
            if size(amp_i, 1) == 1
                amp_i = repmat(amp_i, length(el_i), 1);
            end

            n_arr = SHTools.getACNOrderArray(order);
            
            Alpha = 4 * pi * bsxfun(@times, const_i.^n_arr, conj(SHTools.getComplexSH(order, el_i, az_i))) * amp_i;
        end
        
        function Alpha = alphaFarField3DPerSource(order, el_i, az_i, amp_i, polarity)
            % Alpha = alphaFarField3D(order, el_i, az_i, amp_i, polarity)
            % order = scaler
            % el_i => incident elevations, K x 1
            % az_i => incident azimuths, K x 1
            % amp_i = incident amplitudes, scalar or K(source) x 1
            % polarity = 'sink' or 'source'. Practical mics use sink
            % Alpha = Spherical harmonics coefficients for plane wave, (order+1)^2 x K

            if(strcmp(polarity, 'sink'))
                const_i = 1i;
            else
                const_i = -1i;
            end

            el_i = el_i(:);
            az_i = az_i(:);
            amp_i = amp_i(:);

            n_arr = SHTools.getACNOrderArray(order);
            
            Alpha = 4 * pi * bsxfun(@times, bsxfun(@times, const_i.^n_arr, conj(SHTools.getComplexSH(order, el_i, az_i))), amp_i(:).');
        end
        
        function Alpha = alphaNearField3D(order, k, r_s, el_s, az_s, amp_s)
            % Alpha = alphaNearField3D(order, k, r_s, el_s, az_s, amp_s)
            % order = scaler
            % k = Wave number, scalar
            % r_s => Source radius, sclar or K x 1
            % el_s => incident elevations, K x 1
            % az_s => incident azimuths, K x 1
            % amp_s = incident amplitudes, scalar or K(source) x T (time snapshot)
            % Alpha = Spherical harmonics coefficients for plane wave, (order+1)^2 x T


            el_s = el_s(:);
            az_s = az_s(:);
            
            if size(amp_s, 1) == 1
                amp_s = repmat(amp_s, length(el_s), 1);
            end

            n_arr = SHTools.getACNOrderArray(order);
            
            Alpha = 1i * k * bsxfun(@times, SHTools.sph_h1(n_arr, k * r_s), conj(SHTools.getComplexSH(order, el_s, az_s))) * amp_s;
        end
        
        function Alpha = alphaOmniMicArray(P, order, kr, el, az, isRigid, Mode, compensateBesselZero, applyRegularization, weight)
            % Alpha = alphaOmniMicArray(P, order, kr, el, az, isRigid, Mode, compensateBesselZero, applyRegularization, weight)
            % P = Sound pressure matrix, Q (microphones) x T (time frames)
            % order = Intended sound field order, scaler. Should be <= kr
            % kr = Wave number * array radius, scaler/vector
            % el => Micropone elevations, vector
            % az => Micropone azimuths, vector
            % isRigid = 1/0 => Array type, default 0
            % Mode = inv/mm (pseudoinverse or mode matching)
            % compensateBesselZero = Compensate Bessel zero (uused if mode ~= inv), ref: https://ieeexplore.ieee.org/document/8357900
            % applyRegularization = Apply regularization (used if mode == inv)
            % weight= Array Weight (used if mode ~= inv)
            % Alpha = Sound field coefficients, (N+1)^2 x T

            % Important if true order of the sound field > N
            N_true = max(order, max(ceil(kr)));
            
            n_arr = SHTools.getACNOrderArray(N_true);
            
            % NN x Q (or NN x 1)
            bn = SHTools.sph_bn(n_arr, kr, isRigid);
            
            if strcmp(Mode, 'inv')
                % We must consider N_true, not N, for pinv
                
                % NN x Q
                Y_mat = bsxfun(@times, SHTools.getComplexSH(N_true, el, az), bn);
            
                if applyRegularization
                    Alpha = SHTools.pinvRegularized(Y_mat.', P);
                else
                    Alpha = pinv(Y_mat.') * P;
                end
                
                if order < N_true
                    Alpha = Alpha(1:(order+1)^2, :);
                end
            else
                % It's enough to use N for mode-matching case
                bn = bn(1:(order+1)^2, :);
                
                % Fixed flooring on bn
                if compensateBesselZero
                    bn = max(abs(bn), .05) .* exp(1i*angle(bn));
                end
            
                % NN x Q
                Y_mat = bsxfun(@rdivide, conj(SHTools.getComplexSH(order, el, az)), bn);
                
                if weight ~= 1
                    Y_mat = bsxfun(@times, Y_mat, weight(:).');
                end
                
                Alpha = Y_mat * P;
            end
            
            if (max(ceil(kr)) + 1)^2 > length(az)
                warning('(true_order + 1)^2 > #mics, spatial alising might occur.');
            end
        end

        function evenAlpha = evenAlphaPlanarArrayWithCenterMic(P0, P, order, kr, az, compensateBesselZero, applyRegularization)
            % evenAlpha = evenAlphaPlanarArrayWithCenterMic(P0, P, order, kr, az, compensateBesselZero, applyRegularization)
            % el is implicitly pi/2
            % ref: https://ieeexplore.ieee.org/abstract/document/8547121
            % P0 = Pressure at center mic, 1 x T
            % P = Sound pressure matrix, Q (microphones) x T (time frames)
            % order = Sound field order, scaler
            % kr = Wave number * array radius, scaler
            % az => Micropone azimuths, vector
            % compensateBesselZero = Compensate Bessel zero, 0/1. , ref: https://ieeexplore.ieee.org/document/8357900
            % applyRegularization = Apply regularization
            % evenAlpha = Even sound field coefficients, (N+1)(N+2)/2 x T

            % Importent if true order of the sound field > N
            N_true = max(order, ceil(kr));
            
            % Q x T
            P = bsxfun(@minus, P, P0 * SHTools.sph_besselj(0, kr));

            Q = size(P, 1);
            NN = (N_true + 1) * (N_true + 2) / 2;

            if Q < NN - 1
                warning('Under-determined set of equations');
            end

            az = az(:);

            [~, n_even_arr] = SHTools.getACNOrderArray(N_true);
            
            bn = SHTools.sph_besselj(n_even_arr(2:end), kr);
            
            % Fixed flooring on bn
            if compensateBesselZero
                bn = max(abs(bn), .05) .* exp(1i*angle(bn));
            end
            
            % NN x Q
            Y_mat = SHTools.getComplexEvenSH(N_true, pi/2, az);
            
            % (NN-1) x Q
            Y_mat = bsxfun(@times, bn, Y_mat(2:end, :));

            if applyRegularization
                evenAlpha = SHTools.pinvRegularized(Y_mat.', P);
            else
                % (NN-1) x L
                evenAlpha = pinv(Y_mat.') * P;
            end
            
            % NN x L
            evenAlpha = [sqrt(4*pi) * P0; evenAlpha];
            
            if order < N_true
                evenAlpha = evenAlpha(1:(order + 1) * (order + 2) / 2, :);
            end
        end
        
        function Beta = betaNearField3D(order, k, r_s, el_s, az_s, amp_s)
            % Beta = betaNearField3D(order, k, r_s, el_s, az_s, amp_s)
            % order = scaler
            % k = Wave number, scalar
            % r_i => Source radius, sclar or K x 1
            % el_i => Source elevations, K x 1
            % az_i => Source azimuths, K x 1
            % amp_i = Source amplitudes, scalar or K(source) x T (time snapshot)
            % Beta = Spherical harmonics coefficients, (order+1)^2 x T


            el_s = el_s(:);
            az_s = az_s(:);
            
            if size(amp_s, 1) == 1
                amp_s = repmat(amp_s, length(el_s), 1);
            end

            n_arr = SHTools.getACNOrderArray(order);
            
            Beta = 1i * k * bsxfun(@times, SHTools.sph_bj(n_arr, k * r_s), conj(SHTools.getComplexSH(order, el_s, az_s))) * amp_s;
        end
        
        %% Spherical harmonics synthesis
        function P = soundPressureFromMeasuredAlpha(Alpha, kr, el, az, isRigid)
            % P = soundPressureFromMeasuredAlpha(Alpha, kr, el, az, isRigid)
            % Alpha = Spherical harmonics coefficients, (N+1)^2 x T
            % kr = Wave number * Micropone radius, scaler/vector Q x 1
            % el => Micropone elevations, Q x 1
            % az => Micropone azimuths, Q x 1
            % isRigid = 1/0 => Array type
            % P = Sound pressure, Q x T
            
            NN = size(Alpha, 1);
            N = sqrt(NN) - 1;
            
            if ~floor(N) == N
                error('Fractional order is not a valid option.');
            elseif N < ceil(kr)
                warning('Not enough modes in Alpha for given kr.');
            elseif NN > length(az)
                warning('Not enough mics in the array for given order.');
            end

            n_arr = SHTools.getACNOrderArray(N);
            
            % NN x Q (or NN x 1)
            bn = SHTools.sph_bn(n_arr, kr, isRigid);
            
            % NN x Q
            Y_mat = bsxfun(@times, SHTools.getComplexSH(N, el, az), bn);
            
            % Q x T
            P = Y_mat.' * Alpha;
        end
        
        function P = nearFieldNonReverbSoundPressure(order, k, r, el, az, r_s, el_s, az_s, amp_s, isRigid)
            % P = nearFieldNonReverbSoundPressure(order, k, r, el, az, r_s, el_s, az_s, amp_s, isRigid)
            % order = Sound field order, scaler
            % k = Wave number
            % r = Micropone radius, scaler/Q x 1
            % el => Micropone elevations, Q x 1
            % az => Micropone azimuths, Q x 1
            % r_s => Source radius, scalar or K x 1
            % el_s => Source elevations, K x 1
            % az_s => Source azimuths, K x 1
            % amp_s = Source amplitudes, scalar or K x T
            % isRigid = 1/0 => Array type (only for an interior field)
            % P = Sound pressure Q x T
            
            if max(r) < max(r_s)
%                 disp('r < r_s, hence, using interior field formulation.');
                
                if order < ceil(k * r)
                    order = ceil(k * r);
                    warning(['N is too low, resetting N = ceil(k * r) = ', num2str(order)]);
                end
                
                n_arr = SHTools.getACNOrderArray(order);
                
                % NN x Q (or NN x 1)
                bn = SHTools.sph_bn(n_arr, k*r, isRigid);

                % NN x Q
                Y_mat = bsxfun(@times, SHTools.getComplexSH(order, el, az), bn);

                % NN x T
                Coeff = SHTools.alphaNearField3D(order, k, r_s, el_s, az_s, amp_s);
            elseif max(r) > max(r_s)
%                 disp('r > r_s, hence, using exterior field formulation.');
                
                if order < ceil(k * r_s)
                    order = ceil(k * r_s);
                    warning(['N is too low, resetting N = ceil(k * r_s) = ', num2str(order)]);
                end
                
                n_arr = SHTools.getACNOrderArray(order);
                
                % NN x Q (or NN x 1)
                bn = SHTools.sph_h1(n_arr, k*r);

                % NN x Q
                Y_mat = bsxfun(@times, SHTools.getComplexSH(order, el, az), bn);

                % NN x T
                Coeff = SHTools.betaNearField3D(order, k, r_s, el_s, az_s, amp_s);
            else
                error('r and r_s can not have the same value');
            end
            
            % Q x T
            P = Y_mat.' * Coeff;
        end
        
        function P = farFieldNonReverbSoundPressure(order, kr, el, az, el_i, az_i, amp_i, isRigid, polarity)
            % P = farFieldNonReverbSoundPressure(order, kr, el, az, el_i, az_i, amp_i, isRigid, polarity)
            % order = Sound field order, scaler
            % kr = Wave number * Micropone radius, scaler/Q x 1
            % el => Micropone elevations, Q x 1
            % az => Micropone azimuths, Q x 1
            % el_i => incident elevations, K x 1
            % az_i => incident azimuths, K x 1
            % amp_i = incident amplitudes, scalar or K x T
            % isRigid = 1/0 => Array type
            % polarity = 'sink' or 'source'. Practical mics use sink
            % P = Sound pressure Q x T
            
            if order < ceil(kr)
                warning('N is too low, resetting N = ceil(kr);');
                order = ceil(kr);
            end
            
            n_arr = SHTools.getACNOrderArray(order);
            
            % NN x Q (or NN x 1)
            bn = SHTools.sph_bn(n_arr, kr, isRigid);
            
            % NN x Q
            Y_mat = bsxfun(@times, SHTools.getComplexSH(order, el, az), bn);
            
            % NN x T
            Alpha = SHTools.alphaFarField3D(order, el_i, az_i, amp_i, polarity);
            
            % Q x T
            P = Y_mat.' * Alpha;
        end
        
        %% Helper functions
        function x = pinvRegularized(A, b, weight)
            % x = pinvRegularize(A, b)
            % Solves Ax = b for x
            % A = [M x N]
            % b = [M x K]
            % weight = Regulrization weight
            % x = [N x K]
            
            if ~exist('weight', 'var')
                weight = .01;
            end
            
            % Regularization weight, N x N
            w = weight * eye(size(A, 2));

            % (M + N) x N
            A = [A; w];

            % (M + N) x K
            b = padarray(b, size(w, 1), 'post');
            
            x = pinv(A) * b;
        end
        
        function [n_arr, n_even_arr] = getACNOrderArray(order)
            % [n_arr, n_even_arr] = getACNOrderArray(order)
            % Get all SH orders in an array according to ACN
            n_arr = zeros((order+1)^2, 1);
            n_even_arr = zeros((order+1) * (order+2) / 2, 1);
            ind = 2;
            ind_even = 2;
            for n = 1:order
                numModes = 2*n;
                n_arr(ind:ind+numModes) = n;
                n_even_arr(ind_even:ind_even+n) = n;
                ind = ind + numModes + 1;
                ind_even = ind_even + n + 1;
            end
        end
        
        function NM = getACNOrdersModesArray(order)
            % NM = getACNOrderArray(order)
            % Get all SH orders and modes in an array according to ACN
            % NM = (order+1)^2 x 2
            NM = zeros((order+1)^2, 2);
            
            ind = 2;
            for n = 1:order
                numModes = 2*n;
                NM(ind:ind+numModes, 1) = n;
                NM(ind:ind+numModes, 2) = -n:n;
                
                ind = ind + numModes + 1;
            end
        end
        
        function [N1, N2, k] = getSHOrder(r, f)
            % [N1, N2, k] = sph_N(r, f)
            % f = frequency
            % r = Array radius
            % k = 2 pi f / c
            % N1 = kr
            % N2 = k e r/2
            k = SHTools.getK(f);
            N1 = ceil(k * r);
            N2 = ceil(k * exp(1) * r /2);
        end
        
        function k = getK(freq, c)
            if nargin < 2
                c = 343;
            end
            k = 2 * pi * freq / c;
        end
        
        function order = getOrder(freq, r, c)
            if nargin < 3
                order = ceil(SHTools.getK(freq, 343) * r);
            else
                order = ceil(SHTools.getK(freq, c) * r);
            end
        end
        
        function [X, x, y, z] = s2c(el, az, r)
            % X = s2c(el, az, r)
            % az = azimuth measured from positive x-axis on xy plane (0 to 2pi)
            % el = elevation measured from positive z axis (0 to pi)
            % X = [x, y, z] = Cartesian coordinates
            

            az = az(:);
            el = el(:);
            r = r(:);
            
            x = r .* sin(el) .* cos(az);
            y = r .* sin(el) .* sin(az);
            
            [M, N] = size(x);
            
            z = r .* cos(el) .* ones(M, N);
            
            X = [x, y, z];
        end
        
        function [el, az, r] = c2s2(varargin)
            % [el, az, r] = c2s2(varargin)
            % az = azimuth measured from positive x-axis on xy plane (0 to 2pi)
            % el = elevation measured from positive z axis (0 to pi)
            % varargin = Cartesian coordinates [x, y, z] or [X]
            % x, y, z = [Nx1]
            % X = [N x 3];

            if nargin == 1
                X = varargin{1};
                x = X(:, 1);
                y = X(:, 2);
                z = X(:, 3);
            elseif nargin == 3
                x = varargin{1};
                y = varargin{2};
                z = varargin{3};
            else
                error('Wrong input format');
            end

            r = sqrt(x .^ 2 + y .^ 2 + z .^ 2); 
            el = acos(z ./ r);
            az = atan2(y,x);

            az(az < 0) = az(az < 0) + (2 * pi);
        end
        
        function d_angle = angularDistance(p1, p2, isDegree)
            % d_angle = angularDistance(p1, p2)
            % p1 = [el1, az1], p2 = [el2, az2] (default: radian)
            % d_angle = Angular distance between p1 and p2 (default: radian)
            % isDegree = true/false
            
            if exist('isDegree', 'var') && isDegree
                d_angle = acosd( cosd(p1(1)) * cosd(p2(1)) + sind(p1(1)) * sind(p2(1)) .* ...
                    cosd(p1(2) - p2(2)) );
            else
                d_angle = acos( cos(p1(1)) * cos(p2(1)) + sin(p1(1)) * sin(p2(1)) * ...
                    cos(p1(2) - p2(2)) );
            end
        end

        %% Bessel & Hankel functions
        function bn = sph_bn(n_arr, x_arr, isRigid)
            % bn = sph_bn(n_arr, x_arr, isRigid)
            % Supports random numbers of n and x
            % n_arr = all orders
            % x_arr = all kr
            % isRigid = 0/1
            % bn = [length(n_arr) x length(x_arr)]
            n_arr = n_arr(:);
            x_arr = x_arr(:).';
            
            N = length(n_arr);
            K = length(x_arr);
            
            if N > 1
                x_arr = repmat(x_arr, N, 1);
            end
            
            if K > 1
                n_arr = repmat(n_arr, 1, K);
            end
            
            if exist('isRigid', 'var') && isRigid
                temp = SHTools.dsph_besselj(n_arr, x_arr) .* SHTools.sph_hankel1(n_arr, x_arr) ./ SHTools.dsph_hankel1(n_arr, x_arr);
                temp(isnan(temp)) = 0;
            else
                temp = 0;
            end
            
            bn = SHTools.sph_besselj(n_arr, x_arr) - temp;
        end
        
        function j = sph_bj(n_arr, x_arr)
            % j = sph_bj(n_arr, x_arr)
            % Supports random numbers of n and x
            % j = [length(n_arr) x length(x_arr)];
            
            n_arr = n_arr(:);
            x_arr = x_arr(:).';
            
            N = length(n_arr);
            K = length(x_arr);
            
            if N > 1
                x_arr = repmat(x_arr, N, 1);
            end
            
            if K > 1
                n_arr = repmat(n_arr, 1, K);
            end
            
            j = SHTools.sph_besselj(n_arr, x_arr);
        end
        
        function h1 = sph_h1(n_arr, x_arr)
            % h1 = sph_h1(n_arr, x_arr)
            % Supports random numbers of n and x
            % h1 = [length(n) x length(x)];
            
            n_arr = n_arr(:);
            x_arr = x_arr(:).';
            
            N = length(n_arr);
            K = length(x_arr);
            
            if N > 1
                x_arr = repmat(x_arr, N, 1);
            end
            
            if K > 1
                n_arr = repmat(n_arr, 1, K);
            end
            
            h1 = SHTools.sph_hankel1(n_arr, x_arr);
        end
        
        
        
        
        function j = sph_besselj(n, x)
            %SPH_BESSELJ Spherical bessel function of the first kind.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % SPH_BESSELJ.M - 15/7/2013
            % Archontis Politis, archontis.politis@aalto.fi
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % find zeros in the argument
            idx_zero = (x==0);

            j = sqrt(pi./(2*x)).*besselj(n+0.5, x);
            j(n==0 & idx_zero) = 1;
            j(n~=0 & idx_zero) = 0;

        end
        
        function j = sph_bessely(n, x)
            %SPH_BESSELY Spherical bessel function of the second kind.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % SPH_BESSELY.M - 15/7/2013
            % Archontis Politis, archontis.politis@aalto.fi
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            j = sqrt(pi./(2*x)).*bessely(n+0.5, x);
        end
        
        function j = sph_hankel1(n, x)
            %SPH_HANKEL1 Spherical hankel function of the first kind.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % SPH_HANKEL1.M - 15/7/2013
            % Archontis Politis, archontis.politis@aalto.fi
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            j = SHTools.sph_besselj(n,x) + 1i*SHTools.sph_bessely(n,x);

        end
        
        function dj = dsph_besselj(n, x)
            %DSPH_BESSELJ Spherical bessel function derivative of the first kind.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % DSPH_BESSELJ.M - 15/7/2013
            % Archontis Politis, archontis.politis@aalto.fi
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dj = 1 ./ (2 .* n + 1) .* (n .* SHTools.sph_besselj(n-1,x) - (n+1) .* SHTools.sph_besselj(n+1,x));
            % Teutch:
            % dj = sph_besselj(n-1,x) - (n+1)/x*sph_besselj(n,x);

        end
        
        function dj = dsph_bessely(n, x)
            %DSPH_BESSELY Spherical bessel function derivative of the second kind.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % DSPH_BESSELY.M - 15/7/2013
            % Archontis Politis, archontis.politis@aalto.fi
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dj = 1 ./ (2 .* n+1) .* (n .* SHTools.sph_bessely(n-1,x) - (n+1) .* SHTools.sph_bessely(n+1,x));

        end
        
        function dj = dsph_hankel1(n, x)
            %DSPH_HANKEL1 Spherical hankel function derivative of the first kind.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % DSPH_HANKEL1.M - 15/7/2013
            % Archontis Politis, archontis.politis@aalto.fi
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dj = SHTools.dsph_besselj(n,x) + 1i .* SHTools.dsph_bessely(n,x);

        end
        
        function hom = getEigenmike()
            % hom = getEigenmike()
            % hom.gain = Calibrated gains for the Eigenmic in ANU AASP Lab
            % hom.cart = mic positions in Cartesian coordinates
            % hom.dist = Distance between mics
            
            hom.rigid = 1;
            hom.r = 0.042;
            hom.el = [1.2043; 1.5708; 1.9373; 1.5708; 0.5585; 0.9599; 1.5708; 2.1817; 2.5831; 2.1817; 1.5708; 0.9599; 0.3665; 1.0123; 2.1118; 2.7751; 1.2043; 1.5708; 1.9373; 1.5708; 0.5585; 0.9599; 1.5708; 2.1817; 2.5831; 2.1817; 1.5708; 0.9599; 0.3665; 1.0123; 2.1293; 2.7751];
            hom.az = [0; 0.5585; 0; 5.7247; 0; 0.7854; 1.2043; 0.7854; 0; 5.4978; 5.0789; 5.4978; 1.5882; 1.5708; 1.5708; 1.5533; 3.1416; 3.7001; 3.1416; 2.5831; 3.1416; 3.927; 4.3459; 3.927; 3.1416; 2.3562; 1.9373; 2.3562; 4.6949; 4.7124; 4.7124; 4.7298];
            hom.cart = SHTools.s2c(hom.el, hom.az, hom.r);
            hom.dist = squareform(pdist(hom.cart));
            hom.w = 1;

            % Calibrated gain of the Eigenmike at ANU AASP lab
            hom.gain = [1;0.841783705978144;0.824422770634782;1.07805634428050;0.819368701472337;0.934155064445381;0.942818510614020;0.876545989186322;1.10264739127246;0.848483639170549;1.04461665900460;0.936951042253091;0.922256134152427;1.04479150187472;0.930715983737848;1.27539079246131;1.13043435730242;1.16837435166048;0.862191048272407;0.861692088104175;1.00844381347731;0.754351026481360;0.973011316700164;0.781515901903590;0.787410743970355;0.735444383147764;0.703536416735698;0.738175096862642;0.747236005060125;0.758104189173352;0.742033871663043;1.02480445158795];
        end
        
        function hom = reducedEigenmike(Q, hom)
            % hom = reducedEigenmike(Q, hom)
            % Pick Q microphoens from Eigenmic based on http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html
            % Require HOAGrid.m
            
            if ~exist('hom', 'var')
                hom = SHTools.getEigenmike();
            end
            
            a = hom.cart;

            b_grid = HOAGrid.getGrid(Q);

            b = b_grid.cart;

            ind = zeros(size(b, 1), 1);

            for i = 1:size(b, 1)
                Dist = pdist2(a, b(i, :));
                [~, ind(i)] = min(Dist);
            end
            
            hom.el = hom.el(ind);
            hom.az = hom.az(ind);
            hom.cart = hom.cart(ind, :);
            hom.dist = hom.dist(ind, ind);
            hom.gain = hom.gain(ind);
            hom.w = b_grid.w;
            hom.emReference = ind;
        end
    end
    
    methods
        
    end
end