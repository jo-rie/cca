function [ranks] = copula_rank_transform(data)
%COPULA_RANK_TRANSFORM Rank-transforms the given data and divides it by (n+1). 
% Ties are broken at random
    n = size(data, 1);
    ranks = tiedrank(data, 1) / (n+1);
end

