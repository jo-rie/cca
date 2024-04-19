classdef AbstractParametricCopulaCA < handle
% A class to model Copula CA for parametric copula of the Families Gauss, Gumbel, Frank and Clayton.

    properties (Access = protected)
        copFamily
        copTau
    end

    methods
        function obj = AbstractParametricCopulaCA(copFamily, copTau)
            % Constructor
            obj.copFamily = copFamily;
            obj.copTau = copTau;
        end

        function copFamily = getCopFamily(obj)
            copFamily = obj.copFamily;
        end

        function copTau = getCopTau(obj)
            copTau = obj.copTau;
        end
    end
end
