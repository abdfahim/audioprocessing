clear all
%% Notes
% 1. Include Path to SHTools.m class
% 2. Enter "help SHTools.functionName" to get help about that particular function
% 3. All angles are in radians

%% Input parameters

% Lets say we have 10 microphones with the following elevations and azimuths
Q = 32;
el = pi * rand(Q, 1);
az = 2 * pi * rand(Q, 1);
r = .042 * rand(Q, 1);

% And 5 sources  with the following amplitudes, radii, elevations, and azimuths
% The function can calculate for multiple timeframes at a time considering stationary sources
K = 5;
timeFrames = 500;
el_s = pi * rand(K, 1);
az_s = 2 * pi * rand(K, 1);
r_s = 2 *  rand(K, 1);
amp_s = 10 *  rand(K, timeFrames); % Let say we have 500 time frames

% Spherical harmonics order
order = 4;

% Frequency
f = 1000;

%% Get some parameters

% wavenumber
k = SHTools.getK(f);

%% Get Spherical Harmonics

% Get all real SH upto order = 4 for all azimuths and elevations
y1 = SHTools.getRealSH(order, el, az);  % [25 x 32] output matrix

% Get real SH for n = 3, m = 2 for all azimuths and elevations
y2 = SHTools.getRealSHPerMode(3, 2, el, az);  % [32 x 1] output vector

% Get all complex SH upto order = 4 for all azimuths and elevations
y3 = SHTools.getComplexSH(order, el, az);  % [25 x 32] output matrix

% Get complex SH for n = 3, m = 2 for all azimuths and elevations
y4 = SHTools.getComplexSHPerMode(3, 2, el, az);  % [32 x 1] output vector

% Get all EVEN complex SH upto order = 4 for all azimuths and elevations
y5 = SHTools.getComplexEvenSH(order, el, az);  % [15 x 32] output matrix


%% Calculate analytical alphas (free-field)

% Get alpha for n = 3, m = 2 for all sources considering far-field
% Last parameter takes 'sink'/'source', just change in polarity
% It accumulates alpha contributed by all the sources per time frame
Alpha1 = SHTools.alphaFarField3DPerMode(3, 2, el_s, az_s, amp_s, 'sink'); % [1 x 500] matrix

% And, if we want for all available modes (again, accumulating alpha for
% all sources)
Alpha2 = SHTools.alphaFarField3D(order, el_s, az_s, amp_s, 'sink'); % [25 x 500] matrix

% Instead, if we want to calculate alphas individually for each source (for
% all available modes). THIS TIME, we have to process PER TIME FRAME
Alpha3 = SHTools.alphaFarField3DPerSource(order, el_s, az_s, amp_s(:, 1), 'sink'); % [25 x 5] matrix


% For nearfield, (all available modes upto order, in all time frames)
% Alpha is accumualted for all sources
Alpha4 = SHTools.alphaNearField3D(order, k, r_s, el_s, az_s, amp_s); % [25 x 500] matrix


%% Spherical Harmonic decomposition (get alpha from measured sound pressure)

% Lets say we have measured pressure P
P = rand(Q, timeFrames) + 1i* rand(Q, timeFrames);

% Get Alpha for all all modes and time frames, [25 x 500].
% Type "help SHTools.alphaOmniMicArray" to know input parameter details
Alpha5 = SHTools.alphaOmniMicArray(P, order, k * r, el, az, true, 'inv', false, false, 1);


%% Special case: get even alpha using circular array and a center mic, based on https://ieeexplore.ieee.org/abstract/document/8547121
% Output [15 x 500]. Type "help SHTools.evenAlphaPlanarArrayWithCenterMic" for detials
P0 = rand(1, timeFrames) + 1i* rand(1, timeFrames); % Center microphone
r_circ = .042; % only circular array can be used for this particular function
Alpha6 = SHTools.evenAlphaPlanarArrayWithCenterMic(P0, P, order, k * r_circ, az, false, false);


%% Get analytical beta (free-field), [25 x 500]
Beta = SHTools.betaNearField3D(order, k, r_s, el_s, az_s, amp_s);


%% Spherical harmonics synthesis (get sound pressure from alpha/beta)
% Outputs [32 x 500] for each case (as 32 mics and 500 time frames)

% Sound pressure from measured/calculated Alpha for a RIGID array, for all micrphones and all time frames
P1 = SHTools.soundPressureFromMeasuredAlpha(Alpha2, k * r, el, az, true);

% Sound pressure from measured/calculated Alpha for an OPEN array, for all micrphones and all time frames
% The function automatically decides whether to use interior or exterior formula
P2 = SHTools.soundPressureFromMeasuredAlpha(Alpha2, k * r, el, az, false);

% Analytic (free-field) near-field sound pressure for a RIGID array, for all micrphones and all time frames
% The function automatically decides whether to use interior or exterior formula
P3 = SHTools.nearFieldNonReverbSoundPressure(order, k, r, el, az, r_s, el_s, az_s, amp_s, true);

% Analytic (free-field) near-field sound pressure for an OPEN array, for all micrphones and all time frames
P4 = SHTools.nearFieldNonReverbSoundPressure(order, k, r, el, az, r_s, el_s, az_s, amp_s, false);

% Analytic (free-field) far-field sound pressure for a RIGID array, for all micrphones and all time frames
P5 = SHTools.farFieldNonReverbSoundPressure(order, k * r, el, az, el_s, az_s, amp_s, true, 'sink');

% Analytic (free-field) far-field sound pressure for an OPEN array, for all micrphones and all time frames
P6 = SHTools.farFieldNonReverbSoundPressure(order, k * r, el, az, el_s, az_s, amp_s, false, 'sink');


%% Bessel and hankel functions

% Get bn. Any number of order and kr is supported
isRigid = true;
orders_to_get = [0, 1, 4, 6];
bn = SHTools.sph_bn(orders_to_get, k * r, isRigid); % [4 x 32]

% Get jn. Any number of order and kr is supported
bj = SHTools.sph_bj(orders_to_get, k * r); % [4 x 32]

% Get hankel function. Any number of order and kr is supported
h1 = SHTools.sph_h1(orders_to_get, k * r); % [4 x 32]


%% Some additional useful functions

% Get Eigenmic coordinates
hom = SHTools.getEigenmike();

% Pick Q best microphoens from Eigenmic based on http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html
% Requires HOAGrid.m and HOAGridData.mat
hom1 = SHTools.reducedEigenmike(Q, hom);

% Get spherical microphone array with Q microphones based on http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html
% Requires HOAGrid.m and HOAGridData.mat
hom2 = HOAGrid.getGrid(Q);

% Get all available modes (and even modes) for a particular order
NM = SHTools.getACNOrdersModesArray(order);

% Get truncation limit and k based on radius and frequency
% N1 uses kr formulat, N2 uses k e r/2 formula
[N1, N2, k] = SHTools.getSHOrder(r, f);

% Spherical to cartesian converter
[X, x, y, z] = SHTools.s2c(el, az, r);

% Cartesian to spherical converter
[el1, az1, r1] = SHTools.c2s2(X);




