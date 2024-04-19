classdef CopulaCAEmpirical < CopulaCA
    %CopulaCAEmpirical models given data as Copula Correspondence Model

    properties
        data_ranked
    end

    methods
        function obj = CopulaCAEmpirical(data, n, normalization_method)
            arguments
                data
                n
                normalization_method = 'rank'
            end
            %EmpiricalCA Construct an instance of class EmpiricalCA
            % data: 2-column matrix of the data points, that are put based
            % on a rank transformation into a frequency table
            % n: grid size in u and v direction
            m = size(data, 1);
            if normalization_method == 'rank'
                data_ranked = copula_rank_transform(data);
                discPdf = zeros(n, n);
                for i = 1:m
                    index_row = ceil(data_ranked(i, :) * n);
                    discPdf(index_row(1), index_row(2)) = ...
                        discPdf(index_row(1), index_row(2)) + n./m;
                end
            else
                error(sprintf('Normalization method %s not implemented'));
            end
            discCenteredPdf = discPdf - ones(n, n) ./ n;
            obj@CopulaCA(discCenteredPdf, n);
            obj.data_ranked = data_ranked;
        end

    end

end

