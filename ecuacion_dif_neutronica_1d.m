function [phi_inc, keff] = ecuacion_dif_neutronica_1d(L,D,sigma_a,nu_sigma_f,N,grado_l)

%%%Solución analítica
lamda = nu_sigma_f/(D*(pi/L)^2+sigma_a);
phi = @(x) sin(pi.*x/L);

%%%Matrices
h=L/N;
format long
x_nodos=linspace(0,L,N*(grado_l)+1);

% matriz de polinomios de lagrange

ref_nodos=linspace(-1,1,grado_l+1);


syms x
Lag=sym(zeros(grado_l+1,1));
for i = 1:grado_l+1
    Li=sym(1);
    for k=1:grado_l+1
        if k ~= i
            Li=Li *(x-ref_nodos(k))/(ref_nodos(i)-ref_nodos(k));
        end
    end   
    Lag(i)=Li;
end

dLag = diff(Lag, x);

%matrices locales y ensamblaje

A = sparse(N*(grado_l)+1,N*(grado_l)+1);
B = sparse(N*(grado_l)+1,N*(grado_l)+1);

[xx,ww]=GaussQuad(grado_l+1);

Lag_eval = double(subs(Lag, x, xx'));
dLag_eval = double(subs(dLag, x, xx'));

Jac=h/2;

for e=1:N
    Ke=zeros(grado_l+1,grado_l+1);
    Me=zeros(grado_l+1,grado_l+1);
    Te=zeros(grado_l+1,grado_l+1);

    idx_global = (e-1)*grado_l + (1:grado_l+1);

    for i=1:grado_l+1
        for j=1:grado_l+1
            
            Kij = 0;
            for k = 1:grado_l+1
                Kij = Kij + ww(k) * D * dLag_eval(i,k) * dLag_eval(j,k);
            end
            Ke(i,j) = Kij / Jac;
            
            
            Mij = 0;
            for k = 1:grado_l+1
                Mij = Mij + ww(k) * sigma_a * Lag_eval(i,k) * Lag_eval(j,k);
            end
            Me(i,j) = Mij * Jac;
            
            
            Tij = 0;
            for k = 1:grado_l+1
                Tij = Tij + ww(k) * nu_sigma_f * Lag_eval(i,k) * Lag_eval(j,k);
            end
            Te(i,j) = Tij * Jac;

            I = idx_global(i);
            J = idx_global(j);
            A(I,J) = A(I,J) + Ke(i,j) + Me(i,j);
            B(I,J) = B(I,J) + Te(i,j);
        end
    end
end


% Aplicar condiciones de contorno (Dirichlet homogéneas: phi(0)=phi(L)=0)
A(1,:) = 0; A(1,1) = 1; B(1,:) = 0;
A(end,:) = 0; A(end,end) = 1; B(end,:) = 0;



%%%Solución autovalores
format long
[phi_inc, keff] = eigs(B, A, 1);
phi_inc = abs(phi_inc / max(abs(phi_inc)));



%% Plot
figure
hold on
grid on
xlabel('x (cm)');
ylabel('\phi(x)');
title('Solucion DF')
plot(x_nodos, phi_inc, 'b--x', 'LineWidth', 2);
plot(x_nodos, phi(x_nodos), 'r', 'LineWidth', 2);
end