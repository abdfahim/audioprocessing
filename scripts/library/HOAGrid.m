%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get HOA nodes based on http://www.mathematik.uni-dortmund.de/lsx/research/projects/fliege/nodes/nodes.html
% Requires HOAGridData.mat
% Node files from P. N. Samarasinghe
% Abdullah Fahim, abdullah.fahim@hotmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef HOAGrid < handle
    methods (Static)
        function hom = getGrid(Q)
            N = ceil(sqrt(Q) - 1);
            if N > 29
                error('Upto 900 nodes are supported');
            end
            a = load('HOAGridData.mat');
            hom = a.HOAGrid{N};
        end
    end
end