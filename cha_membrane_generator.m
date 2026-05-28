function cha_membrane_generator()
%CHA_MEMBRANE_GENERATOR  Gerador paramétrico 3D para membrana de zeólita CHA
% -------------------------------------------------------------------------
%  Disciplina : Modelagem Computacional Aplicada à Engenharia Química (05082)
%  Material   : Zeólita CHA (chabazita)  –  Chabazite-Ca
%
%  Referência cristalográfica:
%    Nature 181, 1794-1796 (1958) | COD 9016278 | AMCSD 0017695
%    Grupo espacial : R -3 m  (n.º 166, setting romboédrico)
%    Símbolo Hall   : -P 3* 2
%
%  Parâmetros da célula unitária (setting ROMBOÉDRICO):
%    a = b = c = 9.4 Å       α = β = γ = 94.3°
%    Z = 12 f.u.              V = 823.197 Å³
%    Fórmula: (Al0.32 Si0.68) O2
%
%  Sítios cristalográficos (posições fracionárias, setting R):
%    T  (Al/Si) :  x=0.34,  y=0.11,  z=−0.12  | ocp: Al=0.32, Si=0.68
%    O1         :  x=0.27,  y=−0.27, z=0.00
%    O2         :  x=0.17,  y=−0.17, z=0.50
%    O3         :  x=0.25,  y=0.25,  z=0.90
%    O4         :  x=0.03,  y=0.03,  z=0.35
%
%  Operações de simetria do grupo R -3 m (setting R, 12 ops):
%    x,y,z      | -x,-z,-y  | -z,-x,-y  | y,x,z
%    y,z,x      | -z,-y,-x  | -x,-y,-z  | x,z,y
%    z,x,y      | -y,-x,-z  | -y,-z,-x  | z,y,x
%
%  Nota: setting romboédrico NÃO usa translações de centramento.
%        A célula já é primitiva (12 operações × sítio único = máx. 12 pos.)
%
%  Tipos de átomo no LAMMPS:
%    1 = Si   (carga +2.0500 e)
%    2 = Al   (carga +1.5750 e)
%    3 = O    (carga −1.0250 e)
%
%  Imperfeições modeladas:
%    • Ocupação mista Al/Si no sítio T  (razão Si/Al ajustável)
%    • Desordem posicional           (vibração/distorção estrutural)
%    • Vacâncias aleatórias          (sítios T e O removidos)
%
%  Saída : arquivo LAMMPS  (atom_style charge)
%          Formato: atom-ID  atom-type  charge  x  y  z
%



Lx          = 100.0;       % Largura    da seção reta (Å)
Ly          = 100.0;       % Comprimento da seção reta (Å)
Lz          =  50.0;       % Espessura total da membrana (Å)

n_pores     =    6;        % Número de macroporos
r_pore      =  7.5;        % Raio de cada macroporo (Å)
dist_type   = 'random';    % Distribuição: 'regular' | 'random'
rng_seed    =   42;        % Semente aleatória (seed)

% Razão Si/Al conforme estrutura de referência:
%   Ocupação Al=0.32, Si=0.68  →  Si/Al = 0.68/0.32 ≈ 2.125
%   Use Inf para silicalita pura
si_al_ratio =  2.125;      % Razão Si/Al (SAR) baseada no CIF de referência

defect_amp  =  0.20;       % Amplitude máx. desordem posicional (Å)
vac_frac    =  0.012;      % Fração de vacâncias  (ex.: 0.012 = 1,2 %)

output_file = 'cha_membrane.lammps';   % Nome do arquivo de saída


%%  CÉLULA UNITÁRIA CHA  –  Setting ROMBOÉDRICO  (COD 9016278)
%%
%%  Parâmetros reticulares:
%%    a = b = c = 9.4 Å       α = β = γ = 94.3°
%%
%%  Matriz de vetores de rede (romboédrico → Cartesiano):
%%    Convencional para α = β = γ:
%%      cx = a * cos(α/2) * 2  (projeção)
%%    Usando a matriz M tal que M * [f1;f2;f3] = [x;y;z]:
%%      a1 = a [1, 0, 0]  (rotacionado para facilitar)
%%      Formula geral para célula romboédrica:
%%        tx = cos(α)
%%        ty = (cos(α) - cos(α)^2) / sin(α)   [= 0 para α=90°, != 0 aqui]
%%        tz = sqrt(1 - tx^2 - ty^2... ) mas usamos a decomposição padrão
%%
%%  Posições fracionárias (setting R, grupo -3m, 12 operações):
%%    T  (Al/Si): (0.34,  0.11, -0.12)  + equivalentes por simetria
%%    O1        : (0.27, -0.27,  0.00)  + equivalentes
%%    O2        : (0.17, -0.17,  0.50)  + equivalentes
%%    O3        : (0.25,  0.25,  0.90)  + equivalentes
%%    O4        : (0.03,  0.03,  0.35)  + equivalentes
%% ══════════════════════════════════════════════════════════════════════════

% ── Parâmetros reticulares (romboédrico, COD 9016278) ──────────────────
a_rh  = 9.4;       % Å  (a = b = c)
al_rh = 94.3;      % °  (α = β = γ)

% Constrói vetores de rede romboédricos no referencial Cartesiano
% Convenção padrão IUCr (a1 ao longo de x):
cosA = cosd(al_rh);
sinA = sind(al_rh);
% a2 no plano xy, a3 geral:
%   a1 = a [1, 0, 0]
%   a2 = a [cos(α), sin(α), 0]
%   a3 = a [cos(α), (cos(α)-cos²(α))/sin(α), tz]
%   tz = sqrt(1 - cos²(α) - ((cos(α)-cos²(α))/sin(α))²)
a_vec = a_rh * [1, 0, 0];
b_vec = a_rh * [cosA, sinA, 0];
t_cy  = (cosA - cosA^2) / sinA;
t_cz  = sqrt(max(0, 1 - cosA^2 - t_cy^2));
c_vec = a_rh * [cosA, t_cy, t_cz];

fprintf('━━━ Vetores de rede (romboédrico) ━━━━━━━━━━━\n');
fprintf('   a1 = [%7.4f  %7.4f  %7.4f] Å\n', a_vec);
fprintf('   a2 = [%7.4f  %7.4f  %7.4f] Å\n', b_vec);
fprintf('   a3 = [%7.4f  %7.4f  %7.4f] Å\n', c_vec);
fprintf('   |a| = %.4f  |b| = %.4f  |c| = %.4f Å\n', ...
        norm(a_vec), norm(b_vec), norm(c_vec));
fprintf('   α  = %.2f°  β = %.2f°  γ = %.2f°\n', al_rh, al_rh, al_rh);
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');

% ── Operações de simetria R -3 m (setting R, 12 ops) ──────────────────
%  Extraídas diretamente do CIF:
%   1: x,y,z      2: -x,-z,-y   3: -z,-x,-y   4: y,x,z
%   5: y,z,x      6: -z,-y,-x   7: -x,-y,-z   8: x,z,y
%   9: z,x,y     10: -y,-x,-z  11: -y,-z,-x  12: z,y,x
sym_ops = {
    @(x,y,z) [ x,  y,  z];   %  1
    @(x,y,z) [-x, -z, -y];   %  2
    @(x,y,z) [-z, -x, -y];   %  3
    @(x,y,z) [ y,  x,  z];   %  4
    @(x,y,z) [ y,  z,  x];   %  5
    @(x,y,z) [-z, -y, -x];   %  6
    @(x,y,z) [-x, -y, -z];   %  7
    @(x,y,z) [ x,  z,  y];   %  8
    @(x,y,z) [ z,  x,  y];   %  9
    @(x,y,z) [-y, -x, -z];   % 10
    @(x,y,z) [-y, -z, -x];   % 11
    @(x,y,z) [ z,  y,  x];   % 12
};

% ── Posições assimétricas (sítios Wyckoff, CIF COD 9016278) ───────────
asym_T  = [0.34,  0.11, -0.12];   % sítio T (Al/Si compartilhado)
asym_O1 = [0.27, -0.27,  0.00];
asym_O2 = [0.17, -0.17,  0.50];
asym_O3 = [0.25,  0.25,  0.90];
asym_O4 = [0.03,  0.03,  0.35];

% Aplica operações de simetria a cada sítio assimétrico
T_frac  = apply_sym(asym_T,  sym_ops);
O1_frac = apply_sym(asym_O1, sym_ops);
O2_frac = apply_sym(asym_O2, sym_ops);
O3_frac = apply_sym(asym_O3, sym_ops);
O4_frac = apply_sym(asym_O4, sym_ops);

% Conversão frac → Cartesiano
T_cart  = frac2cart(T_frac,  a_vec, b_vec, c_vec);
O1_cart = frac2cart(O1_frac, a_vec, b_vec, c_vec);
O2_cart = frac2cart(O2_frac, a_vec, b_vec, c_vec);
O3_cart = frac2cart(O3_frac, a_vec, b_vec, c_vec);
O4_cart = frac2cart(O4_frac, a_vec, b_vec, c_vec);

nT_uc = size(T_cart,1);
nO_uc = size(O1_cart,1) + size(O2_cart,1) + size(O3_cart,1) + size(O4_cart,1);

fprintf('━━━ Célula unitária CHA (COD 9016278) ━━━━━━\n');
fprintf('   T  (Al/Si): %2d posições  (asym: 0.34, 0.11,-0.12)\n', nT_uc);
fprintf('   O1         : %2d posições  (asym: 0.27,-0.27, 0.00)\n', size(O1_cart,1));
fprintf('   O2         : %2d posições  (asym: 0.17,-0.17, 0.50)\n', size(O2_cart,1));
fprintf('   O3         : %2d posições  (asym: 0.25, 0.25, 0.90)\n', size(O3_cart,1));
fprintf('   O4         : %2d posições  (asym: 0.03, 0.03, 0.35)\n', size(O4_cart,1));
fprintf('   Total      : %2d átomos/célula\n', nT_uc + nO_uc);
fprintf('   Ocupação   : Al=0.32, Si=0.68  →  Si/Al ≈ 2.125\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');

% Monta arrays da célula unitária (tipo T provisoriamente = 1/Si)
uc_xyz  = [T_cart; O1_cart; O2_cart; O3_cart; O4_cart];
uc_type = [ones(nT_uc,1); 3*ones(nO_uc,1)];   % 1=Si(T), 3=O


%%                        REPLICAÇÃO DA REDE

rng(rng_seed);

% Alcance das réplicas ao longo de cada eixo romboédrico
% Usa o comprimento máximo da projeção de cada vetor sobre x,y,z
max_proj = max(abs([a_vec; b_vec; c_vec]), [], 1);  % projeções máximas
nx_rep = ceil(Lx / max_proj(1)) + 4;
ny_rep = ceil(Ly / max_proj(2)) + 4;
nz_rep = ceil(Lz / max_proj(3)) + 2;

% Pré-alocação generosa
n_est = (nx_rep+8) * (ny_rep+8) * (nz_rep+4) * (nT_uc + nO_uc);
all_xyz  = zeros(n_est, 3);
all_type = zeros(n_est, 1);
ptr = 0;
n_uc = nT_uc + nO_uc;

for ix = -3 : nx_rep+2
    for iy = -3 : ny_rep+2
        for iz = 0 : nz_rep+1
            shift = ix*a_vec + iy*b_vec + iz*c_vec;
            idx   = ptr+1 : ptr+n_uc;
            all_xyz(idx, :) = uc_xyz + shift;
            all_type(idx)   = uc_type;
            ptr = ptr + n_uc;
        end
    end
end

all_xyz  = all_xyz(1:ptr, :);
all_type = all_type(1:ptr);

% Recorte: mantém apenas átomos dentro de [0,Lx] × [0,Ly] × [0,Lz]
tol  = 1e-4;
mask = all_xyz(:,1) >= -tol & all_xyz(:,1) <= Lx+tol & ...
       all_xyz(:,2) >= -tol & all_xyz(:,2) <= Ly+tol & ...
       all_xyz(:,3) >= -tol & all_xyz(:,3) <= Lz+tol;

all_xyz  = all_xyz(mask,:);
all_type = all_type(mask);

% Clipa para [0, L]
all_xyz(:,1) = min(max(all_xyz(:,1), 0), Lx);
all_xyz(:,2) = min(max(all_xyz(:,2), 0), Ly);
all_xyz(:,3) = min(max(all_xyz(:,3), 0), Lz);

fprintf('Átomos após replicação + recorte : %d\n', size(all_xyz,1));


%%          SUBSTITUIÇÃO ALEATÓRIA  Si → Al  (ocupação mista do sítio T)
%%
%%  Referência: Al_occ = 0.32, Si_occ = 0.68  (COD 9016278)
%%  Si/Al = 0.68/0.32 = 2.125  (valor padrão em si_al_ratio)
%%  Regra de Löwenstein aplicada de forma estatística (distribuição aleatória)

si_idx = find(all_type == 1);
n_Si   = numel(si_idx);
n_Al   = round(n_Si / (si_al_ratio + 1));
n_Al   = min(n_Al, n_Si);
al_sel = si_idx(randperm(n_Si, n_Al));
all_type(al_sel) = 2;   % 2 = Al


%%                     POSICIONAMENTO DOS MACROPOROS

pore_xy = place_pores(n_pores, r_pore, Lx, Ly, dist_type, rng_seed);
n_pores_real = size(pore_xy, 1);

% Remove átomos no interior dos macroporos (cilindros paralelos a z)
keep = true(size(all_xyz,1), 1);
for p = 1:n_pores_real
    r2 = (all_xyz(:,1)-pore_xy(p,1)).^2 + (all_xyz(:,2)-pore_xy(p,2)).^2;
    keep(r2 <= r_pore^2) = false;
end
all_xyz  = all_xyz(keep,:);
all_type = all_type(keep);


%%          IMPERFEIÇÕES ESTRUTURAIS  (membrana não-perfeita)
%%
%%  (1) Desordem posicional : desloca aleatoriamente ~30 % dos átomos
%%      → simula distorção da rede, vibrações congeladas, desordem Si/Al
%%  (2) Vacâncias           : remove aleatoriamente fração de sítios
%%      → defeitos pontuais, sítios Si/Al/O ausentes
%%  (3) Confinamento        : garante posições dentro da caixa periódica

N = size(all_xyz, 1);

% (1) Desordem posicional
n_disp = round(N * 0.30);
d_idx  = randperm(N, n_disp);
all_xyz(d_idx,:) = all_xyz(d_idx,:) + defect_amp * (2*rand(n_disp,3) - 1);

% (2) Vacâncias
n_vac = round(N * vac_frac);
v_idx = randperm(N, n_vac);
keep2 = true(N,1);
keep2(v_idx) = false;
all_xyz  = all_xyz(keep2,:);
all_type = all_type(keep2);

% (3) Confinamento periódico (wrap)
all_xyz(:,1) = mod(all_xyz(:,1), Lx);
all_xyz(:,2) = mod(all_xyz(:,2), Ly);
all_xyz(:,3) = mod(all_xyz(:,3), Lz);

N = size(all_xyz, 1);


%%            CARGAS PARCIAIS  (modelo BKS-like para aluminossilicato)
%%
%%  Nota: framework puro-Si é neutro (q_Si + 2·q_O = 0).
%%  Com substituição Al, o framework fica carregado negativamente:
%%  cátions compensadores (Ca²⁺, conforme Chabazite-Ca) devem ser
%%  adicionados em etapa posterior no LAMMPS.

charges = zeros(N,1);
charges(all_type == 1) =  2.0500;   % Si
charges(all_type == 2) =  1.5750;   % Al
charges(all_type == 3) = -1.0250;   % O

% Carga total do framework
q_total = sum(charges);
fprintf('Carga total do framework : %+.2f e\n', q_total);
if abs(q_total) > 1
    fprintf('  ⚠  Adicionar %.2f e em cátions Ca²⁺ compensadores no LAMMPS.\n\n', -q_total);
end


write_lammps(output_file, all_xyz, all_type, charges, Lx, Ly, Lz, ...
             pore_xy, r_pore, dist_type, si_al_ratio, defect_amp, vac_frac, ...
             a_rh, al_rh);


fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║         RESUMO DA GEOMETRIA  (CHA / COD 9016278) ║\n');
fprintf('╠══════════════════════════════════════════════════╣\n');
fprintf('║  Arquivo de saída  : %-26s║\n', output_file);
fprintf('║  Referência CIF    : COD 9016278 / AMCSD 0017695 ║\n');
fprintf('║  Grupo espacial    : R -3 m  (n.º 166, R setting)║\n');
fprintf('║  Célula unitária   : a=b=c=9.4 Å  α=β=γ=94.3°   ║\n');
fprintf('║  Caixa (Å)         : %5.1f × %5.1f × %5.1f          ║\n', Lx, Ly, Lz);
fprintf('║  Átomos totais     : %-26d║\n', N);
fprintf('║    ▸ Si (tipo 1)   : %-26d║\n', sum(all_type==1));
fprintf('║    ▸ Al (tipo 2)   : %-26d║\n', sum(all_type==2));
fprintf('║    ▸ O  (tipo 3)   : %-26d║\n', sum(all_type==3));
fprintf('║  Razão Si/Al       : %-26.4f║\n', si_al_ratio);
fprintf('║  Ocp. referência   : Al=0.32, Si=0.68            ║\n');
fprintf('║  Macroporos        : %d (r=%.1f Å, %-8s)      ║\n', ...
        n_pores_real, r_pore, dist_type);
fprintf('║  Vacâncias         : %.1f %%                       ║\n', vac_frac*100);
fprintf('╚══════════════════════════════════════════════════╝\n\n');

visualize_membrane(all_xyz, all_type, pore_xy, r_pore, Lx, Ly, Lz);

end   

function C = frac2cart(F, a1, a2, a3)
%FRAC2CART  Coordenadas fracionárias → Cartesianas
%   F(N,3) = [f1 f2 f3];   a1,a2,a3 = vetores de rede (1×3)
C = F(:,1)*a1 + F(:,2)*a2 + F(:,3)*a3;
end

% ─────────────────────────────────────────────────────────────────────────
function F_out = apply_sym(frac0, sym_ops)
%APPLY_SYM  Aplica operações de simetria a uma posição fracionária
%
%  Entradas:
%    frac0   : [x, y, z] posição assimétrica (1×3)
%    sym_ops : cell array de function handles {@(x,y,z) ...}
%
%  Saída:
%    F_out : posições únicas após aplicação de todas as operações (N×3)
%            as coordenadas são reduzidas para [0,1) com mod(·,1)

x0 = frac0(1);  y0 = frac0(2);  z0 = frac0(3);
n_ops = numel(sym_ops);
raw   = zeros(n_ops, 3);

for k = 1:n_ops
    v = sym_ops{k}(x0, y0, z0);
    raw(k,:) = v;
end

% Reduz para [0,1) e elimina duplicatas com tolerância
raw_mod = mod(raw, 1.0);
raw_mod = round(raw_mod * 1e7) / 1e7;
F_out   = unique(raw_mod, 'rows');
end

% ─────────────────────────────────────────────────────────────────────────
function pore_xy = place_pores(n_pores, r_pore, Lx, Ly, dist_type, seed)
%PLACE_PORES  Determina centros (x,y) dos macroporos
%
%  'regular' : grade regular com nc × nr pontos
%  'random'  : posicionamento aleatório com seed, evitando sobreposição

rng(seed);
margin = max(r_pore * 1.5, r_pore + 5.0);

if strcmp(dist_type, 'regular')
    nc = ceil(sqrt(n_pores));
    nr = ceil(n_pores / nc);
    xv = linspace(margin, Lx-margin, nc);
    yv = linspace(margin, Ly-margin, nr);
    [gx, gy] = meshgrid(xv, yv);
    pts = [gx(:), gy(:)];
    pore_xy = pts(1:min(n_pores, size(pts,1)), :);

else   % random
    pore_xy = zeros(n_pores, 2);
    k = 0;
    min_dist = 2.2 * r_pore;
    for trial = 1:200000
        if k >= n_pores, break; end
        cx = margin + (Lx - 2*margin) * rand();
        cy = margin + (Ly - 2*margin) * rand();
        if k == 0
            ok = true;
        else
            d = sqrt((pore_xy(1:k,1)-cx).^2 + (pore_xy(1:k,2)-cy).^2);
            ok = all(d > min_dist);
        end
        if ok
            k = k + 1;
            pore_xy(k,:) = [cx, cy];
        end
    end
    if k < n_pores
        warning('Apenas %d/%d poros posicionados (espaço insuficiente para r=%.1f Å).', ...
                k, n_pores, r_pore);
    end
    pore_xy = pore_xy(1:k, :);
end
end

% ─────────────────────────────────────────────────────────────────────────
function write_lammps(fname, xyz, types, charges, Lx, Ly, Lz, ...
                      pore_xy, r_pore, dtype, sar, d_amp, v_frac, a_rh, al_rh)
%WRITE_LAMMPS  Escreve arquivo de dados LAMMPS
%
%  Formato atom_style charge:
%     atom-ID   atom-type   charge    x   y   z

N = size(xyz, 1);
active_types = sort(unique(types))';
n_types = numel(active_types);

masses_all = [28.0860, 26.9820, 15.9990];   % Si, Al, O
type_names = {'Si','Al','O'};

fid = fopen(fname, 'w');
if fid < 0
    error('Não foi possível criar: %s', fname);
end

% ── Cabeçalho ──────────────────────────────────────────────────────────
fprintf(fid, '# ════════════════════════════════════════════════════════\n');
fprintf(fid, '# Membrana de Zeólita CHA (Chabazite-Ca)\n');
fprintf(fid, '# Gerado por: cha_membrane_generator.m\n');
fprintf(fid, '# ────────────────────────────────────────────────────────\n');
fprintf(fid, '# Referência : COD 9016278 / AMCSD 0017695\n');
fprintf(fid, '#              Nature 181, 1794-1796 (1958)\n');
fprintf(fid, '# Grupo esp. : R -3 m  (Hall: -P 3* 2)  n.º 166\n');
fprintf(fid, '# Setting    : Romboédrico (primitivo)\n');
fprintf(fid, '# Célula uc  : a=b=c=%.2f Å   α=β=γ=%.1f°\n', a_rh, al_rh);
fprintf(fid, '# Fórmula    : (Al0.32 Si0.68) O2  (Z=12)\n');
fprintf(fid, '# Sítios T   : (0.34, 0.11,-0.12)  ocp: Al=0.32, Si=0.68\n');
fprintf(fid, '# O1         : (0.27,-0.27, 0.00)\n');
fprintf(fid, '# O2         : (0.17,-0.17, 0.50)\n');
fprintf(fid, '# O3         : (0.25, 0.25, 0.90)\n');
fprintf(fid, '# O4         : (0.03, 0.03, 0.35)\n');
fprintf(fid, '# ────────────────────────────────────────────────────────\n');
fprintf(fid, '# Caixa      : %.2f x %.2f x %.2f Å\n', Lx, Ly, Lz);
fprintf(fid, '# Macroporos : %d  |  raio: %.2f Å  |  dist: %s\n', ...
        size(pore_xy,1), r_pore, dtype);
fprintf(fid, '# Si/Al      : %.4f  (ref. CIF: 2.125)\n', sar);
fprintf(fid, '# Imperfeições: desordem %.2f Å  |  vacâncias %.1f %%\n', d_amp, v_frac*100);
fprintf(fid, '# ════════════════════════════════════════════════════════\n\n');

% ── Contagens ──────────────────────────────────────────────────────────
fprintf(fid, '%d atoms\n',        N);
fprintf(fid, '%d atom types\n\n', n_types);

% ── Dimensões da caixa ─────────────────────────────────────────────────
fprintf(fid, '0.000000  %.6f  xlo xhi\n', Lx);
fprintf(fid, '0.000000  %.6f  ylo yhi\n', Ly);
fprintf(fid, '0.000000  %.6f  zlo zhi\n\n', Lz);

% ── Massas ─────────────────────────────────────────────────────────────
fprintf(fid, 'Masses\n\n');
for t = active_types
    fprintf(fid, '  %d  %.4f  # %s\n', t, masses_all(t), type_names{t});
end

% ── Posições atômicas ──────────────────────────────────────────────────
fprintf(fid, '\nAtoms  # charge\n\n');
for i = 1:N
    fprintf(fid, '%8d  %d  %9.4f  %13.6f  %13.6f  %13.6f\n', ...
            i, types(i), charges(i), xyz(i,1), xyz(i,2), xyz(i,3));
end

fclose(fid);
fprintf('Arquivo LAMMPS escrito com sucesso: %s\n', fname);
end

% ─────────────────────────────────────────────────────────────────────────
function visualize_membrane(xyz, types, pore_xy, r_pore, Lx, Ly, Lz)
%VISUALIZE_MEMBRANE  Gera três figuras de diagnóstico

% Cores e tamanhos por tipo de átomo
clr = [0.20 0.45 0.90;   % Si – azul
       0.85 0.35 0.10;   % Al – laranja
       0.95 0.20 0.20];  % O  – vermelho
sz_3d = [25, 30, 7];
lbl   = {'Si (tipo 1)','Al (tipo 2)','O (tipo 3)'};

% Subsample para visualização rápida (máx. 12000 átomos por tipo)
max_vis = 12000;
if size(xyz,1) > max_vis*3
    vis_idx = randperm(size(xyz,1), min(max_vis*3, size(xyz,1)));
else
    vis_idx = 1:size(xyz,1);
end
xyz_v  = xyz(vis_idx,:);
type_v = types(vis_idx);

% ── Figura 1: Vista 3D ────────────────────────────────────────────────
f1 = figure('Name','CHA Membrane – Vista 3D','Color','w','Position',[40 80 950 660]);
ax1 = axes(f1); hold(ax1,'on'); grid(ax1,'on');

for t = 1:3
    idx = type_v == t;
    if any(idx)
        scatter3(ax1, xyz_v(idx,1), xyz_v(idx,2), xyz_v(idx,3), ...
                 sz_3d(t), clr(t,:), 'filled', ...
                 'MarkerFaceAlpha', 0.45, 'DisplayName', lbl{t});
    end
end

% Contornos dos macroporos
theta = linspace(0, 2*pi, 72);
for p = 1:size(pore_xy,1)
    for zz = linspace(0, Lz, 8)
        xp = pore_xy(p,1) + r_pore*cos(theta);
        yp = pore_xy(p,2) + r_pore*sin(theta);
        plot3(ax1, xp, yp, zz*ones(1,72), 'k-', 'LineWidth',0.5, 'HandleVisibility','off');
    end
    plot3(ax1, pore_xy(p,1)*[1 1], pore_xy(p,2)*[1 1], [0 Lz], ...
          'k--', 'LineWidth',1.0, 'HandleVisibility','off');
end

xlabel(ax1,'x (Å)'); ylabel(ax1,'y (Å)'); zlabel(ax1,'z (Å)');
title(ax1, 'Zeólita CHA (COD 9016278) – Membrana 3D (subamostra visual)', ...
      'FontSize',12,'FontWeight','bold');
legend(ax1,'Location','northeast','FontSize',10);
view(ax1, 35, 22); set(ax1,'FontSize',11); axis(ax1,'tight');

% ── Figura 2: Vista superior (plano xy) ──────────────────────────────
f2 = figure('Name','CHA Membrane – Vista Superior','Color','w','Position',[40 80 700 620]);
ax2 = axes(f2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');

t_mask = (type_v==1) | (type_v==2);
scatter(ax2, xyz_v(t_mask,1), xyz_v(t_mask,2), 5, [0.78 0.78 0.84], ...
        'filled','HandleVisibility','off');

theta2 = linspace(0,2*pi,120);
for p = 1:size(pore_xy,1)
    xp = pore_xy(p,1) + r_pore*cos(theta2);
    yp = pore_xy(p,2) + r_pore*sin(theta2);
    fill(ax2, xp, yp, [0.10 0.10 0.10], 'FaceAlpha',0.30, ...
         'EdgeColor','k', 'LineWidth',1.5, 'DisplayName',sprintf('Poro %d',p));
    text(ax2, pore_xy(p,1), pore_xy(p,2), num2str(p), ...
         'HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
end

xlabel(ax2,'x (Å)'); ylabel(ax2,'y (Å)');
title(ax2,'Vista Superior – Distribuição dos Macroporos','FontSize',12,'FontWeight','bold');
legend(ax2,'Location','northeast','FontSize',9);
axis(ax2,'equal'); xlim(ax2,[0 Lx]); ylim(ax2,[0 Ly]);
set(ax2,'FontSize',11);

% ── Figura 3: Composição e distribuição z ────────────────────────────
f3 = figure('Name','CHA – Análise Composicional','Color','w','Position',[40 80 820 380]);

% Subplot 1: contagem de átomos por tipo
ax3a = subplot(1,2,1,'Parent',f3);
cnts = [sum(types==1), sum(types==2), sum(types==3)];
b = bar(ax3a, 1:3, cnts, 0.6);
b.FaceColor = 'flat';
b.CData = clr;
set(ax3a,'XTickLabel',{'Si','Al','O'},'FontSize',11);
xlabel(ax3a,'Tipo de Átomo'); ylabel(ax3a,'Número de Átomos');
title(ax3a,'Composição Atômica','FontSize',11,'FontWeight','bold');
grid(ax3a,'on');
for k=1:3
    if cnts(k) > 0
        text(ax3a, k, cnts(k)+0.015*max(cnts), num2str(cnts(k)), ...
             'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    end
end

% Subplot 2: distribuição de átomos ao longo da espessura (z)
ax3b = subplot(1,2,2,'Parent',f3);
edges = linspace(0, Lz, 31);
for t = 1:3
    idx = types == t;
    if any(idx)
        histogram(ax3b, xyz(idx,3), edges, 'FaceColor', clr(t,:), ...
                  'FaceAlpha',0.55,'EdgeColor','none','DisplayName',lbl{t});
        hold(ax3b,'on');
    end
end
xlabel(ax3b,'z (Å)'); ylabel(ax3b,'Átomos por camada');
title(ax3b,'Distribuição ao longo de z','FontSize',11,'FontWeight','bold');
legend(ax3b,'Location','northeast','FontSize',8);
grid(ax3b,'on'); set(ax3b,'FontSize',11);

sgtitle(f3,'Zeólita CHA (COD 9016278) – Análise da Membrana Gerada', ...
        'FontSize',13,'FontWeight','bold');
end
