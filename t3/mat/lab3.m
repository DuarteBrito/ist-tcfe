close all
clear all

%% parametros
n = 7.29150093083668;  %----------
C = 0.0000175150173805745;    %----------
R = 45e3;     %----------
%R2 = 7e3;  %------
r_d = 65.7894736842105;  %----------
regulator = 19;   %----------
%Von = 12.0001/regulator;      %----------
Von = 0.631508421052632;


%% Vs

%n = 6.673844;  %----------
f = 50;
w = 50*2*pi;
%Von = 0.7;     %----------
A = 230/n;
rect = 2;


t= linspace(0, 1/(2*f) , 1000);

vs = abs(A*cos(w*t))-rect*Von;
vs(vs<0) = 0;

figure
plot(t,vs)
legend({'vs'},'Location','southwest')

%% vC
%comeca igual a vs
vC = vs;

%descarregamento
%C = 4.508784e-6;    %----------
%R = 2.219563e3;     %----------

R_ = 1/(R + (regulator+rect)*r_d);
%R_=1/R;
t_off = (1/w) * atan(1/w/C*R_);
v_exp = (A-rect*Von)*cos(w*t_off)*exp(-(t-t_off)/C*R_);

figure
plot(t,vs, t_off*ones(1,round(max(vs))), 1:round(max(vs)), t, v_exp)
legend({'vs','t_off','v_(exp)'},'Location','southwest')

%vc
vC(t>t_off) = max([vC(t>t_off);v_exp(t>t_off)]);
vC_control = [max(vC) min(vC)];

figure
plot(t, vC)
legend({'v_C'},'Location','southwest')

%Vc e vc
Vc = mean(vC)*ones(1, length(vC));
vc = vC-Vc;

figure
plot(t, Vc, t , vc)
legend({'V_c','v_c'},'Location','southwest')

%r_d = Von/((1e-9)*exp(17));

%% vO

%Vo
%regulator = 17;   %----------
Vo = ones(1, length(t)) * (Von*regulator);

%vo
%r_d = 80;  %----------

vo = (r_d*regulator)/(r_d*regulator+R) * vc;

vO = Vo + vo;

figure
plot(t, Vo, t, vo, t, vO)
legend({'V_o','v_o','v_O'},'Location','southwest')

%% resultados
ripple = max(vo)-min(vo)
avg = Vo(1)
cost = R*1e-3 + (rect*2+regulator)*0.1 + C*1e6
merit = 1/(cost*(ripple + (avg-12)+1e-6))

%% graficos
figure
plot(t, vC , t, vO)
title('Envelope Detector and Voltage Regulator')
xlabel('t [s]')
ylabel('V [V]')
legend({'v_C', 'v_O'},'Location','southwest')
print ("results_1.png", "-dpng");
ylim([-1 35])
print ("results_2.png", "-dpng");

figure
plot(t, vO-12)
title('v_O - 12')
xlabel('t [s]')
ylabel('v_O [V]')
legend({'v_O - 12'},'Location','southwest')
print ("results.png", "-dpng");

%% tabelas
fidCirc2 = fopen("resultados_tabela.tex","w");
fprintf(fidCirc2,"Parameter & Value \\\\\n"); 
fprintf(fidCirc2, "\\hline\n");
fprintf(fidCirc2,"Output DC level & %f V \\\\\n", Vo(1));
fprintf(fidCirc2, "\\hline\n");
fprintf(fidCirc2,"Ripple & %f V \\\\\n",ripple);
fprintf(fidCirc2, "\\hline\n");
fprintf(fidCirc2,"Cost & %f MU \\\\\n",cost);
fprintf(fidCirc2, "\\hline\n");
fprintf(fidCirc2,"M & %f \\\\\n",merit);
fprintf(fidCirc2, "\\hline\n");
fclose(fidCirc2);
close all
