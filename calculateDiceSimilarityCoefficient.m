function [ diceSimilarityCoefficient ] = calculateDiceSimilarityCoefficient (groundTruth,sample)
% CALCULATEDICESIMILARITYCOEFFICIENT is calculating the Sorensen-Dice index / DSC
% for two logical masks.
%
%  Input:
%   - groundTruth (logical)
%   - sample (logical)
%
%  Output:
%   - diceSimilarityCoefficient (0-1)
%
%
%  Notes:
%  -----
%  In case of a size mismatch, the sample matrix gets resized to match the
%  ground truth.
%
%
%  References:
%  -----
%  Dice, L. R. (1945). Measures of the amount of ecologic association between
%  species. Ecology, 26(3), 297-302.
%
%
%  Version 2016.12.05
%  Marcel Gehrung, Werner Siemens Imaging Center (2016)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isequal(size(groundTruth),size(sample))
    warning('The size of the ground truth and the sample has to be the same for calculating the dice similarity coefficient. Trying to resize...');
    sample = imresize(sample,size(groundTruth));
    if ~isequal(size(groundTruth),size(sample))
        error('Resizing of sample failed.');
    end
end
diceSimilarityCoefficient = 2*nnz(groundTruth.*sample)/(nnz(groundTruth)+nnz(sample));
end