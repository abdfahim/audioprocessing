%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short-time Fourier Transform Library
% (c) Abdullah Fahim, abdullah.fahim@hotmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef STFTClass < handle
    properties
        windowsize = 0;
        hopsize = 0;
        overlap = 0;
        nfft = 0;
        fs = 0;
        
        T = 0;
        F = 0;
        freqArray = [];
    end
    
    properties (SetAccess = 'private')
        pos_freq = 0;
        hop = 0.5;
        verbose = true;
    end
    
    methods
        function obj = STFTClass(fs, winlen, hop, nfft)
            % winlen = Window length
            %       : in ms if double, e.g., 256
            %       : in #samples if integer, e.g., uint32(256)
            % hop in fraction
            % nfft in #samples
            
            if ~exist('nfft','var')
                nfft = [];
            end
            if ~exist('hop','var')
                hop = 0.5;
            end
            obj.setParams(fs, winlen, hop, nfft);
        end
        
        function setVerbose(obj, v)
            obj.verbose = v;
        end
        
        function setParams(obj, varargin)
            if nargin > 1 && ~isempty(varargin{1})
                obj.fs = varargin{1};
            end
            
            if nargin > 2 && ~isempty(varargin{2})
                if isinteger(varargin{2})
                    obj.windowsize = double(varargin{2});
                else
                    obj.windowsize = floor(varargin{2} * obj.fs / 1000);
                end
            end
            
            if nargin > 3 && ~isempty(varargin{3})
                if isinteger(varargin{3})
                    obj.hopsize = varargin{3};
                    obj.hop = obj.hopsize / obj.windowsize;
                else
                    obj.hop = varargin{3};
                    obj.hopsize = floor(obj.windowsize * obj.hop);
                end
            end
            
            if nargin > 4 && ~isempty(varargin{4})
                obj.nfft = 2^nextpow2(varargin{4});
            else
                obj.nfft = 2^nextpow2(obj.windowsize);
            end
            
            obj.overlap = obj.windowsize - obj.hopsize;
            obj.pos_freq = floor(obj.nfft/2) + 1;
        end
        
        
        function X = stft(obj, x, T, plotSpectrogram)
            % X = stft(obj, x, T, plotSpectrogram)
            % x = input time domain signal
            % T = optional, number of time frames to keep
            % plotSpectrogram = optional, 0/1
            % X = STFT output
            
            if obj.verbose
                disp(['Running STFT with fs = ', num2str(obj.fs), '. If fs has changed, reinitialize the object with correct fs.']);
            end
            
            % Adding obj.overlap zero at the end to avoid distortion
            x(end+1:end+obj.overlap) = 0;
            
            % buffer adds obj.overlap zeros as initial condition.
            % using 'nodelay' would result in distortion at the beginning
            xFrames = buffer(x, obj.windowsize, obj.overlap);
            xFrames = bsxfun(@times, xFrames, hamming(obj.windowsize, 'periodic'));
            X = fft(xFrames, obj.nfft) / obj.windowsize;
            X = X(1:obj.pos_freq, :);
            
            if exist('T', 'var') && T > 0
                if T > size(X, 2)
                    X(:, end+1:T) = 0;
                else
                    X(:, T+1:end) = [];
                end
            end
            
            [obj.F, obj.T] = size(X);
            
            obj.freqArray = (0:obj.pos_freq-1).' * obj.fs / obj.nfft;
            
            
            if exist('plotSpectrogram', 'var') && plotSpectrogram
                obj.plotSpectrogram(X);
            end
        end
        
        function x = istft(obj, X, sigL)
            % x = istft(obj, X, sigL)
            % X = STFT signal
            % sigL = optional, signal length to keep
            % x = time domain signal
            
            if obj.verbose
                disp(['Running ISTFT with fs = ', num2str(obj.fs), '. If fs has changed, reinitialize the object with correct fs.']);
            end
            
            xFrames = ifft(X, obj.nfft, 'symmetric') * obj.windowsize;
            WinFunc = hamming(obj.nfft, 'periodic');
            xFrames = bsxfun(@times, xFrames, WinFunc);
            
            % Overlap-add
            x = zeros((obj.T-1) * obj.hopsize + obj.nfft, 1);
            
            for ind = 1:obj.T
                x(1 + (ind-1) * obj.hopsize : obj.nfft + (ind-1) * obj.hopsize) = ...
                    x(1 + (ind-1) * obj.hopsize : obj.nfft + (ind-1) * obj.hopsize) + ...
                    xFrames(:, ind);
            end
            
            % Normalize
            x = x * obj.hopsize / sum(WinFunc.^2);
            
            % obj.overlap zeros were added at the beginning (by buffer) and
            % end (manually) during stft
            x = x(obj.overlap+1:end-obj.overlap);
            
            if exist('sigL', 'var') && sigL > 0
                x(sigL+1:end) = [];
            end
        end
        
        function TArr = getTimeBins(obj)
            TArr = ((0:obj.T-1) * obj.hopsize + obj.windowsize / 2)/obj.fs;
        end
        
        function FArr = getFrequencyBins(obj)
            FArr = ((0:obj.F - 1) * obj.fs / obj.nfft).';
        end
        
        function plotSpectrogram(obj, X, isBinned)
            if exist('isBinned', 'var') && isBinned
                imagesc(abs(X));
            else
                imagesc(obj.getTimeBins(), obj.getFrequencyBins()/1000, abs(X));
                xlabel('Seconds');
                ylabel('KHz')
            end

            set(gca,'YDir','normal');
            
        end
        
        function [windowsize, hopsize, nfft] = getParams(obj)
            windowsize = obj.windowsize;
            hopsize = obj.hopsize;
            nfft = obj.nfft;
        end
        
        function printParams(obj)
            disp(['Window size: ', num2str(obj.windowsize)]);
            disp(['Hop size: ', num2str(obj.hopsize)]);
            disp(['NFFT: ', num2str(obj.nfft)]);
        end
        
        function freq = getFrequency(obj, freqBinMatlab)
            % freq = getFrequency(obj, freqBinMatlab)
            % Considering postive nfft
            % freqBinMatlab = 1, 2, ...., nfft/2, nfft/2 + 1, nfft/2 + 2, .... nfft
            % 1 => DC => freq = (1-1)/nfft * fs = 0
            % nfft/2 + 1 => Nyquist => freq = (nfft/2 + 1 - 1) / nfft * fs = fs / 2
            % 2 .... nfft/2 => Positive freq (total: nfft/2 - 1)
            % nfft/2 + 2 .... nfft => Negative freq (total: nfft/2 - 1)
            freq = (freqBinMatlab - 1) * obj.fs / obj.nfft;
        end
        
        
        function freqBinMatlab = getFreqBinMatlab(obj, freq)
            % freqBinMatlab = getFreqBinMatlab(obj, freq)
            if max(freq) > obj.fs/2
                error('Frequency cannot exceed fs/2 (Have a problem? Ask Mr. Nyquist).');
            end
            freqBinMatlab = round(freq / obj.fs * obj.nfft) + 1;
        end
    end
end