L = 350;
D = [1.446, 0.776+ zeros(1,12) ,1.446];
sigma_a = [0.0077, 0.0244+ zeros(1,12), 0.0077];
nu_sigma_f = [0, 0.0260+ zeros(1,12), 0];
tamano_celdas=[25, 25 + zeros(1, 12), 25];
N = 14;
grado_l = 3;

malla = Malla(L, N, grado_l, tamano_celdas, D, sigma_a, nu_sigma_f);
elemento = ElementoFinito(grado_l);

problema = ProblemaDifusion1D(malla, elemento);
problema = problema.ensamblar_matrices();
problema = problema.aplicar_cc();
problema = problema.resolver_autovalor();
problema.graficar();

disp(problema.keff)
