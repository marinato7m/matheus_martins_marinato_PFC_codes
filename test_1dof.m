% =========================================================================
% TEST_1DOF.M  —  Passo 1: Verificação do Modelo 1-DOF Vertical
%
% O QUE ESTES TESTES FAZEM (e o que NÃO fazem)
% ─────────────────────────────────────────────
% Estes testes verificam que o CÓDIGO está implementado corretamente:
% sem bugs de sinal, sem erros de unidade, sem problemas numéricos.
%
% Eles NÃO provam que as equações físicas escolhidas são as corretas
% para descrever o foguete real. Essa validação física mais forte
% acontece no Passo 2 (test_3dof.m), onde comparamos nosso modelo
% contra o modelo de ponto de massa das Aulas 21+22 do professor —
% duas implementações independentes do mesmo problema.
%
% TESTE 1 — Queda livre vs. solução analítica
%   O que valida: que o ode45 integra corretamente a ODE dv/dt = −g
%   Como: compara ode45 com z(t) = z0 + vz0·t − ½g·t² (fórmula fechada)
%   Estes SÃO independentes: um integra numericamente, outro é fórmula.
%   Erro esperado: ~10⁻¹³ m (só erro de discretização do ode45)
%
% TESTE 2 — Empuxo sem gravidade vs. Tsiolkovsky
%   O que valida: que o código de empuxo e consumo de massa não tem bugs
%   Limitação: Tsiolkovsky é a solução analítica das MESMAS EDOs que o
%   ode45 integra — então 0% de erro é esperado por construção matemática,
%   não é surpresa. O teste pega bugs de implementação (sinal errado,
%   unidade errada, mdot com sinal trocado, etc.)
%
% TESTE 3 — Empuxo com gravidade vs. Tsiolkovsky + perda gravitacional
%   O que valida: que gravidade e empuxo se somam corretamente no código
%   Mesma limitação do Teste 2: a previsão teórica vem das mesmas EDOs.
%   Valida implementação, não a escolha do modelo físico.
%
% TESTE 4 — Conservação de energia em queda livre
%   O que valida: que o integrador não introduz dissipação espúria
%   Como: E = ½vz² + g·z deve ser constante (lei de conservação)
%   Valida o integrador via lei física — mais forte que Teste 1.
%
% COMO USAR
%   Coloque todos os .m na mesma pasta e execute: test_1dof
% =========================================================================

clear; close all; clc;

p = falcon9_params();

fprintf('\n=========================================================\n');
fprintf('  VERIFICACAO PASSO 1 — MODELO 1-DOF VERTICAL\n');
fprintf('  (verifica implementacao; validacao fisica no Passo 2)\n');
fprintf('=========================================================\n\n');

% ─────────────────────────────────────────────────────────────────────────
% TESTE 1 — Queda livre
% Verifica: ode45 integra dz/dt=vz, dvz/dt=−g corretamente
% Referência: solução analítica fechada z(t) = z0 + vz0·t − ½g·t²
% ─────────────────────────────────────────────────────────────────────────
fprintf('TESTE 1 — Queda livre\n');
fprintf('  Verifica: integrador numerico vs. solucao analitica\n');

z0  = p.h0;   
vz0 = p.vz0;  

% Tempo de impacto: raiz positiva de z0 + vz0·t − (1/2)g·t² = 0
t_hit = (vz0 + sqrt(vz0^2 + 2*p.g*z0)) / p.g;

% Solução fechada amostrada em 1000 pontos (referência independente).
t_vec    = linspace(0, t_hit, 1000);
z_analy  = z0 + vz0*t_vec - 0.5*p.g*t_vec.^2;
vz_analy = vz0 - p.g*t_vec;

% Integração numérica da mesma queda livre (estado [z; vz]).
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t_1, x_1] = ode45(@(t,x) [x(2); -p.g], [0, t_hit], [z0; vz0], opts);

% Reamostra a solução numérica na grade analítica e mede o erro máximo.
z_interp  = interp1(t_1, x_1(:,1), t_vec, 'spline');
err_z_abs = max(abs(z_interp - z_analy));
err_z_rel = err_z_abs / abs(z0) * 100;

fprintf('  t_hit = %.4f s (esperado ~5.27 s)\n', t_hit);
fprintf('  Erro altitude: %.2e m  (%.10f%%)\n', err_z_abs, err_z_rel);
fprintf('  Interpretacao: erro e ruido de discretizacao do ode45\n');
if err_z_rel < 0.01
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% ─────────────────────────────────────────────────────────────────────────
% TESTE 2 — Empuxo sem gravidade
% Verifica: codigo de empuxo e consumo de massa sem bugs de implementacao
% Referência: Tsiolkovsky — solução analítica das mesmas EDOs
% Limitação: 0% de erro é esperado por construção; o teste pega bugs
%            de sinal, unidade ou lógica no código, não erros de modelo
% ─────────────────────────────────────────────────────────────────────────
fprintf('TESTE 2 — Empuxo sem gravidade\n');
fprintf('  Verifica: codigo de T/m e mdot sem bugs de implementacao\n');
fprintf('  Referencia: Tsiolkovsky (solucao analitica das mesmas EDOs)\n');
fprintf('  Nota: 0%% de erro e esperado — o teste pega bugs de codigo\n');

% Para quando a massa atinge mf (fim do propelente).
opts_ev = odeset('Events',  @(t,x) ev_propelente(t, x, p.mf), ...
                 'RelTol',   1e-10, 'AbsTol',   1e-10);

% Estado [posição; velocidade; massa]; g=0 isola o efeito do empuxo.
x0_2 = [0; 0; p.m0];
[t_2, x_2, te_2, xe_2, ~] = ode45( ...
    @(t,x) ode_thrust(t, x, p.T_max, 0.0, p.Isp, p.g0), ...
    [0, 60], x0_2, opts_ev);

if isempty(te_2)
    fprintf('  FALHOU: evento nao disparou\n\n');
else
    % Compara o delta-v simulado com a previsão de Tsiolkovsky.
    mf_real  = xe_2(3);
    dv_sim   = xe_2(2) - x0_2(2);
    dv_tsiol = p.Isp * p.g0 * log(p.m0 / mf_real);
    err2     = abs(dv_sim - dv_tsiol) / dv_tsiol * 100;

    fprintf('  Duracao burn:         %.4f s\n', te_2);
    fprintf('  Massa consumida:      %.2f kg\n', p.m0 - mf_real);
    fprintf('  dv simulado:          %.4f m/s\n', dv_sim);
    fprintf('  dv Tsiolkovsky:       %.4f m/s\n', dv_tsiol);
    fprintf('  Diferenca:            %.8f%%\n', err2);
    if err2 < 0.01
        fprintf('  PASSOU (sem bugs de implementacao)\n\n');
    else
        fprintf('  FALHOU (bug de implementacao detectado)\n\n');
    end
end

% ─────────────────────────────────────────────────────────────────────────
% TESTE 3 — Empuxo com gravidade
% Verifica: gravidade e empuxo somam corretamente no codigo
% Referência: Tsiolkovsky + correção gravitacional (g·t_burn)
% Mesma limitação do Teste 2
% ─────────────────────────────────────────────────────────────────────────
fprintf('TESTE 3 — Empuxo com gravidade\n');
fprintf('  Verifica: gravidade e empuxo somam corretamente no codigo\n');
fprintf('  Referencia: Tsiolkovsky com correcao gravitacional\n');

% Mesmo burn do Teste 2, agora com g ativo (parte de 5000 m para subir).
x0_3 = [5000; 0; p.m0];
[t_3, x_3, te_3, xe_3, ~] = ode45( ...
    @(t,x) ode_thrust(t, x, p.T_max, p.g, p.Isp, p.g0), ...
    [0, 60], x0_3, opts_ev);

if isempty(te_3)
    fprintf('  FALHOU: evento nao disparou\n\n');
else
    % Previsão = Tsiolkovsky menos a perda gravitacional g·t_burn.
    t_burn3   = te_3;
    dv_sim3   = xe_3(2) - x0_3(2);
    dv_tsiol3 = p.Isp * p.g0 * log(p.m0 / xe_3(3));
    dv_perda  = p.g * t_burn3;
    dv_prev   = dv_tsiol3 - dv_perda;
    err3      = abs(dv_sim3 - dv_prev) / abs(dv_prev) * 100;

    fprintf('  Duracao burn:              %.4f s\n', t_burn3);
    fprintf('  dv simulado:               %+.4f m/s\n', dv_sim3);
    fprintf('  dv Tsiolkovsky (sem g):    %+.4f m/s\n', dv_tsiol3);
    fprintf('  dv perda grav. (g x t):    −%.4f m/s\n', dv_perda);
    fprintf('  dv previsto (Tsiol − g·t): %+.4f m/s\n', dv_prev);
    fprintf('  Diferenca:                 %.8f%%\n', err3);
    if err3 < 0.01
        fprintf('  PASSOU (sem bugs de implementacao)\n\n');
    else
        fprintf('  FALHOU (bug de implementacao detectado)\n\n');
    end
end

% ─────────────────────────────────────────────────────────────────────────
% TESTE 4 — Conservação de energia
% Verifica: integrador não introduz dissipação espúria
% Referência: lei de conservação de energia mecânica
% Mais forte que Teste 1 — valida via princípio físico externo
% ─────────────────────────────────────────────────────────────────────────
fprintf('TESTE 4 — Conservacao de energia\n');
fprintf('  Verifica: integrador nao introduz dissipacao espuria\n');
fprintf('  Referencia: E = 1/2*vz^2 + g*z = constante (lei fisica)\n');

% Energia específica ao longo da queda livre do Teste 1; deve ser constante.
e_num = 0.5*x_1(:,2).^2 + p.g*x_1(:,1);
e_var = max(e_num) - min(e_num);
e_ref = abs(0.5*vz0^2 + p.g*z0);
err_e = e_var / e_ref * 100;

fprintf('  Variacao de energia:  %.2e J/kg\n', e_var);
fprintf('  Erro relativo:        %.10f%%\n', err_e);
if err_e < 1e-4
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% ─────────────────────────────────────────────────────────────────────────
fprintf('=========================================================\n');
fprintf('  RESUMO\n');
fprintf('  Teste 1: ode45 integra gravidade corretamente\n');
fprintf('  Teste 2: codigo de empuxo e massa sem bugs\n');
fprintf('  Teste 3: gravidade + empuxo somam corretamente\n');
fprintf('  Teste 4: integrador conserva energia\n');
fprintf('\n');
fprintf('  PROXIMA ETAPA: test_3dof.m\n');
fprintf('  La comparamos nosso modelo contra o modelo de ponto\n');
fprintf('  de massa das Aulas 21+22 — validacao fisica mais forte.\n');
fprintf('=========================================================\n\n');

% ── PLOTS ─────────────────────────────────────────────────────────────────
% Diagnóstico rápido dos Testes 1 e 4 (queda livre e energia).

figure('Name','TESTE 1+4 — Queda Livre', ...
       'NumberTitle','off','Position',[50 50 900 400]);

% Altitude: analítico vs. numérico (devem coincidir).
subplot(1,3,1);
plot(t_vec, z_analy, 'r--', 'LineWidth', 2); hold on;
plot(t_1,   x_1(:,1), 'b-', 'LineWidth', 1.5); grid on;
xlabel('Tempo [s]'); ylabel('z_I [m]');
legend('Analitico','ode45');
title('Teste 1 — Altitude');

% Erro absoluto entre as duas soluções.
subplot(1,3,2);
plot(t_vec, abs(z_interp - z_analy), 'k-', 'LineWidth', 1.5); grid on;
xlabel('Tempo [s]'); ylabel('|Erro| [m]');
title('Teste 1 — Erro: ode45 vs. Analitico');

% Desvio de energia em relação ao valor inicial (deve ficar ~0).
subplot(1,3,3);
e_n = 0.5*x_1(:,2).^2 + p.g*x_1(:,1);
plot(t_1, e_n - e_n(1), 'b-', 'LineWidth', 1.5); grid on;
yline(0, 'r--', 'Esperado = 0', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('\DeltaE [J/kg]');
title('Teste 4 — Conservacao de Energia');

% Diagnóstico dos Testes 2 e 3 (burns), só se ambos dispararam o evento.
if ~isempty(te_2) && ~isempty(te_3)
    figure('Name','TESTES 2+3 — Empuxo (verifica implementacao)', ...
           'NumberTitle','off','Position',[50 530 900 400]);

    % Velocidade durante o burn: sem g (sobe sempre) vs. com g.
    subplot(1,3,1);
    plot(t_2, x_2(:,2), 'b-', 'LineWidth', 2); hold on;
    plot(t_3, x_3(:,2), 'g-', 'LineWidth', 2); grid on;
    yline(0, 'r--', 'LineWidth', 1.2);
    xlabel('Tempo [s]'); ylabel('v_z [m/s]');
    legend('Sem g (T2)','Com g (T3)');
    title('Velocidade durante o Burn');

    % Consumo de massa idêntico nos dois (mesmo empuxo).
    subplot(1,3,2);
    plot(t_2, x_2(:,3), 'b-', 'LineWidth', 2); hold on;
    plot(t_3, x_3(:,3), 'g-', 'LineWidth', 2); grid on;
    yline(p.mf, 'r--', sprintf('m_f=%.0f kg',p.mf), 'LineWidth', 1.2);
    xlabel('Tempo [s]'); ylabel('Massa [kg]');
    legend('Sem g (T2)','Com g (T3)');
    title('Consumo de Massa');

    % Altitude no Teste 3 (com g, foguete sobe).
    subplot(1,3,3);
    plot(t_3, x_3(:,1), 'g-', 'LineWidth', 2); grid on;
    xlabel('Tempo [s]'); ylabel('z_I [m]');
    title('Altitude — Teste 3 (com g, foguete sobe)');
end


% =========================================================================
% PLOTS — Publication-grade figure for the TCC (Chapter 3, V&V)
% Generates: figuras/vv_1dof.pdf  (16 x 7 cm, 1 row x 3 columns)
% =========================================================================
% Versão "limpa" das mesmas validações, formatada para entrar no texto.
% Usa o estilo padronizado de tcc_plot_style/tcc_axes/tcc_save.
S = tcc_plot_style();
 
fig = figure('Name','V&V — 1-DOF battery','NumberTitle','off','Color','w');
 
% ── Panel (a): free fall, ode45 vs analytic ─────────────────────────────
subplot(1,3,1);
plot(t_vec, z_analy, '--', 'Color', S.orange, 'LineWidth', S.lw_ref); hold on;
plot(t_1,   x_1(:,1),  '-', 'Color', S.blue,   'LineWidth', S.lw_data);
tcc_axes;
xlabel('Time [s]',        'FontSize', S.fs_label);
ylabel('Altitude z_I [m]','FontSize', S.fs_label);
title('(a) Free-fall trajectory',   'FontSize', S.fs_title, 'FontWeight', 'normal');
legend({'Analytical','Numerical (ode45)'}, ...
       'Location','southwest', 'FontSize', S.fs_legend, 'Box','off');
 
% ── Panel (b): residual ────────────────────────────────────────────────
% Erro em escala log (eps como piso para o semilogy não quebrar em zero).
subplot(1,3,2);
semilogy(t_vec, max(abs(z_interp - z_analy), eps), '-', ...
         'Color', S.blue, 'LineWidth', S.lw_data);
tcc_axes;
xlabel('Time [s]',                       'FontSize', S.fs_label);
ylabel('|z_{num} - z_{an}|  [m]',        'FontSize', S.fs_label);
title('(b) Integrator residual', 'FontSize', S.fs_title, 'FontWeight', 'normal');
 
% ── Panel (c): mechanical energy conservation ──────────────────────────
e_n = 0.5*x_1(:,2).^2 + p.g*x_1(:,1);
subplot(1,3,3);
plot(t_1, e_n - e_n(1), '-', 'Color', S.blue, 'LineWidth', S.lw_data); hold on;
yline(0, '--', 'Color', S.black, 'LineWidth', S.lw_lim);
tcc_axes;
xlabel('Time [s]',                'FontSize', S.fs_label);
ylabel('\Delta E  [J kg^{-1}]',   'FontSize', S.fs_label);
title('(c) Energy conservation',  'FontSize', S.fs_title, 'FontWeight', 'normal');
 
% Exporta o PDF na largura de página cheia definida no estilo.
tcc_save(fig, 'vv_1dof', S.w_full, 6.5);

% ── FUNÇÕES LOCAIS ─────────────────────────────────────────────────────────

function dxdt = ode_thrust(~, x, T, g_val, Isp, g0)
% ODE 1-DOF com empuxo: estado [posição; velocidade; massa].
% Aceleração = T/m - g; consumo de massa pela vazão T/(Isp*g0).
    m    = x(3);
    az   = T/m - g_val;
    mdot = -T / (Isp * g0);
    dxdt = [x(2); az; mdot];
end

function [value, isterminal, direction] = ev_propelente(~, x, mf)
% Evento de parada: massa (x(3)) cruzando mf de cima para baixo.
    value      = x(3) - mf;
    isterminal = 1;
    direction  = -1;
end