classdef ProblemaDifusion1D
    properties
        malla
        elemento
        A
        B
        phi_inc
        keff
    end
    
    methods
        function s = ProblemaDifusion1D(malla, elemento)
            s.malla = malla;
            s.elemento = elemento;
        end
        
        function s = ensamblar_matrices(s)
            N = s.malla.N;
            grado_l = s.malla.grado_l;
            tamano_celdas = s.malla.tamano_celdas;
            D = s.malla.D;
            sigma_a = s.malla.sigma_a;
            nu_sigma_f = s.malla.nu_sigma_f;
            
            Lag_eval = double(subs(s.elemento.Lag, sym('x'), s.elemento.xx'));
            dLag_eval = double(subs(s.elemento.dLag, sym('x'), s.elemento.xx'));
            Jac = tamano_celdas/2;
            
            
            A = sparse(N*grado_l+1, N*grado_l+1);
            B = sparse(N*grado_l+1, N*grado_l+1);
            
            for e = 1:N
                idx = s.malla.nodos(e, :);
                Ke = zeros(grado_l+1);
                Me = zeros(grado_l+1);
                Te = zeros(grado_l+1);
                for i=1:grado_l+1
                    for j=1:grado_l+1
                        for k=1:grado_l+1
                            Ke(i,j) = Ke(i,j) + s.elemento.ww(k)*D(e)*dLag_eval(i,k)*dLag_eval(j,k);
                            Me(i,j) = Me(i,j) + s.elemento.ww(k)*sigma_a(e)*Lag_eval(i,k)*Lag_eval(j,k);
                            Te(i,j) = Te(i,j) + s.elemento.ww(k)*nu_sigma_f(e)*Lag_eval(i,k)*Lag_eval(j,k);
                        end
                    end
                end
                Ke = Ke / Jac(e);
                Me = Me * Jac(e);
                Te = Te * Jac(e);
                A(idx,idx) = A(idx,idx) + Ke + Me;
                B(idx,idx) = B(idx,idx) + Te;
            end
            
            s.A = A;
            s.B = B;
        end
        
        function s = aplicar_cc(s)
            s.A(1,:) = 0; s.A(1,1)=1; s.B(1,:)=0;
            s.A(end,:) = 0; s.A(end,end)=1; s.B(end,:)=0;
        end
        
        function s = resolver_autovalor(s)
            [phi, k] = eigs(s.B, s.A, 1);
            s.phi_inc = abs(phi / max(abs(phi)));
            s.keff = k;
        end
        
        function s = graficar(s)
            figure; hold on; grid on;
            plot(s.malla.x_nodos, s.phi_inc, 'b--x', 'LineWidth', 2);
            xlabel('x (cm)'); ylabel('\phi(x)');
            title('Solución de difusión neutrónica 1D');
        end
    end
end


