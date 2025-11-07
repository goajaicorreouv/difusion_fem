classdef ElementoFinito
    properties
        ref_nodos
        Lag
        dLag
        xx
        ww
    end
    
    methods
        function s = ElementoFinito(grado_l)
            s.ref_nodos = linspace(-1,1,grado_l+1);
            syms x
            Lag = sym(zeros(grado_l+1,1));
            for i=1:grado_l+1
                Li = sym(1);
                for k=1:grado_l+1
                    if k~=i
                        Li = Li * (x - s.ref_nodos(k)) / (s.ref_nodos(i) - s.ref_nodos(k));
                    end
                end
                Lag(i) = Li;
            end
            s.Lag = Lag;
            s.dLag = diff(Lag, x);
            [s.xx, s.ww] = GaussQuad(grado_l+1);
        end
    end
end