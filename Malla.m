classdef Malla
    properties
        L
        N
        grado_l
        x_nodos
        nodos
        tamano_celdas
        D
        sigma_a
        nu_sigma_f
    end
    
    methods
        function s = Malla(L, N, grado_l, tamano_celdas, D, sigma_a, nu_sigma_f)
            s.L = L;
            s.N = N;
            s.grado_l = grado_l;

            s.nodos = zeros(N , grado_l+1);
            for e = 1:N
                s.nodos(e,:) = (e-1)*grado_l + (1:grado_l+1);
            end
            
            if sum(tamano_celdas) ~= L || length(tamano_celdas) ~= N || length(D) ~= N || length(sigma_a) ~= N || length(nu_sigma_f) ~= N
                error('Tamaño Invalido: Todos los vectores proporcionados deben tener longitud N = %d y la suma de los elementos del vector tamaño debe ser L = %d.', N, L);
            else
                s.tamano_celdas = tamano_celdas;
                s.D = D;
                s.sigma_a = sigma_a;
                s.nu_sigma_f = nu_sigma_f;
            end
            
            x_nodos = linspace(0, tamano_celdas(1), grado_l + 1);  
            x_e = tamano_celdas(1);       
            for e = 2:N
                x_local = linspace(0, tamano_celdas(e), grado_l + 1);        
                x_elem = x_e + x_local;   
                
                x_nodos = [x_nodos, x_elem(2:end)];
                x_e = x_e + tamano_celdas(e);
            end
            s.x_nodos = x_nodos;
        end
    end
end
