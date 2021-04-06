close all
clear all
format longG
pkg load symbolic

%% getting initial data
data = importdata('../data.txt',"=",8);
data = data.data;

%% MATLAB   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% chamar variáveis

G1 = 1/data(1)*1e-3;
G2 = 1/data(2)*1e-3;
G3 = 1/data(3)*1e-3;
G4 = 1/data(4)*1e-3;
G5 = 1/data(5)*1e-3;
G6 = 1/data(6)*1e-3;
G7 = 1/data(7)*1e-3;
Vs = data(8);
C = data(9)*1e-6;
Kb = data(10)*1e-3;
Kd = data(11)*1e3;


 %% circuit 1
            %V1 V2 V3 V5 V6 V7 V8
circ1_no1 = [1,0,0,0,0,0,0];
circ1_no2 = [G1,-G1-G2-G3,G2,G3,0,0,0];
circ1_no3 = [0,G2+Kb,-G2,-Kb,0,0,0];
circ1_no5 = [0,0,0,1,0,Kd*G6,-1];
circ1_no6 = [0,-Kb,0,G5+Kb,-G5,0,0];
circ1_no7 = [0,0,0,0,0,-G6-G7,G7];
circ1_no8 = [0,G3,0,-G3-G4-G5,G5,G7,-G7];

eq_circ1 = [circ1_no1;circ1_no2;circ1_no3;circ1_no5;circ1_no6;circ1_no7;circ1_no8];
b_circ1 = [Vs;0;0;0;0;0;0];

res1 = eq_circ1\b_circ1;


fidCirc1 = fopen("resultados1.txt","w");
fprintf(fidCirc1," ,V (V)\n");
fprintf(fidCirc1,"V1,%f\n",res1(1));
fprintf(fidCirc1,"V2,%f\n",res1(2));
fprintf(fidCirc1,"V3,%f\n",res1(3));
fprintf(fidCirc1,"V5,%f\n",res1(4));
fprintf(fidCirc1,"V6,%f\n",res1(5));
fprintf(fidCirc1,"V7,%f\n",res1(6));
fprintf(fidCirc1,"V8,%f\n",res1(7));
fclose(fidCirc1);

Vxfrom1 = res1(5)-res1(7);  %solucao para o ngspice

 %% circuit 2
Vx = res1(5)-res1(7);
            %V2 V3 V5 V6 V7 V8
circ2_no2 = [-G1-G2-G3,G2,G3,0,0,0];
circ2_no3 = [G2+Kb,-G2,-Kb,0,0,0];
circ2_no5 = [0,0,1,0,Kd*G6,-1];
circ2_no6 = [0,0,0,1,0,-1];
circ3_no7 = [0,0,0,0,-G6-G7,G7];
circ2_no8 = [G3-Kb,0,-G3-G4+Kb,0,G7,-G7];

eq_circ2 = [circ2_no2;circ2_no3;circ2_no5;circ2_no6;circ3_no7;circ2_no8];
b_circ2 = [0;0;0;Vx;0;0];

res2 = eq_circ2\b_circ2;

Ix = (res2(4)-res2(3))/data(5) + data(10)*(res2(1)-res2(3));
R = (res2(4)-res2(6))/Ix;

t_circ2 = 0:0.1e-5:20e-3;
V_6n = (res2(4)-res2(6))*exp((-1/(R * data(9)*1e-3))*t_circ2);


figure
plot(t_circ2,V_6n)

title('Natural Solution')
xlabel('t [s]')
ylabel('V_6_n [V]')
legend('V_6_n')
print ("circ2.png", "-dpng");

fidCirc2 = fopen("resultados2.txt","w");
fprintf(fidCirc2," ,V (V)\n");
fprintf(fidCirc2,"V2,%f\n",res2(1));
fprintf(fidCirc2,"V3,%f\n",res2(2));
fprintf(fidCirc2,"V5,%f\n",res2(3));
fprintf(fidCirc2,"V6,%f\n",res2(4));
fprintf(fidCirc2,"V7,%f\n",res2(5));
fprintf(fidCirc2,"V8,%f\n",res2(6));
fclose(fidCirc2);

%valores para o ngspice
n_fronteira6 = res2(4);
n_fronteira8 = res2(6);



%% circuit 3
f_circ3 = 1000;
w_circ3 = 2*pi*f_circ3;
            %~V2 ~V3 ~V5 ~V6 ~V7 ~V8
circ3_no1 = [1, 0, 0, 0, 0, 0, 0];
circ3_no2 = [-G1, G1 + G3+G2, -G2, -G3, 0, 0, 0];
circ3_no3 = [0, Kb+G2, -G2, -Kb, 0, 0, 0];
circ3_no5 = [0, 0, 0, 1, 0, Kd*G6, -1];
circ3_no6 = [0, Kb, 0, -G5-Kb, G5+1i*w_circ3*C, 0, -1i*w_circ3*C];
circ3_no7 = [0, 0, 0, 0, 0, -G6-G7, G7];
circ3_no8 = [0, -G3, 0, G3+G4+G5, -1i*w_circ3*C-G5, -G7, G7+1i*w_circ3*C];

eq_circ3 = [circ3_no1; circ3_no2; circ3_no3; circ3_no5; circ3_no6; circ3_no7; circ3_no8];

b_circ3 = [exp(-1i*pi/2); 0; 0; 0; 0; 0; 0];

res3 = eq_circ3\b_circ3;

phase = angle(res3(5));
amplitude = abs(res3(5));

t_circ3 = 0:0.1e-5:20e-3;
V_6f = real(amplitude * exp(-1i*(w_circ3*t_circ3 + phase))); 

figure
plot(t_circ3,V_6f)
title('Forced Solution')
legend('V_6_f')
xlabel('t [s]')
ylabel('V_6_f [V]')
print ("circ3.png", "-dpng");

fidtot2 = fopen("resultados4.txt","w");
fprintf(fidtot2," ,V (V)\n");
fprintf(fidtot2,"V1,%f\n",real(res3(1)));
fprintf(fidtot2,"V2,%f\n",real(res3(2)));
fprintf(fidtot2,"V3,%f\n",real(res3(3)));
fprintf(fidtot2,"V5,%f\n",real(res3(4)));
fprintf(fidtot2,"V6,%f\n",real(res3(5)));
fprintf(fidtot2,"V7,%f\n",real(res3(6)));
fprintf(fidtot2,"V8,%f\n",real(res3(7)));
fclose(fidtot2);

%% ponto 5
% somar as duas solucoes anteriores

t_menor = -5e-3:0.1e-5:0-0.1e-5;
vs_menor = ones(1,length(t_menor))*Vs;
v_6_menor = ones(1,length(t_menor))*double(res1(5));

t_maior = 0:0.1e-5:20e-3;
vs_maior = sin(2*pi*f_circ3*t_maior);
v_6_maior = V_6n + V_6f;

t_p5 = [t_menor t_maior];
vs_p5 = [vs_menor vs_maior];
v_6_p5 = [v_6_menor v_6_maior];

figure
plot(t_p5,vs_p5, t_p5,v_6_p5)
title('Final Solution')
legend('V_s', 'V_6')
xlabel('t [s]')
ylabel('V [V]')
print ("tot_ponto5.png", "-dpng");


%% ponto 6
f_s = logspace(-1, 6, 100);
v6_p6 = zeros(1,100);
v6_p6 = zeros(1,100);
v_in = exp(-1i*pi/2);

for k = 1:100
    
w_p6 = 2*pi*f_s(k);

          %~V2 ~V3 ~V5 ~V6 ~V7 ~V8
p6_no1 = [1, 0, 0, 0, 0, 0, 0];
p6_no2 = [-G1, G1 + G3+G2, -G2, -G3, 0, 0, 0];
p6_no3 = [0, Kb+G2, -G2, -Kb, 0, 0, 0];
p6_no5 = [0, 0, 0, 1, 0, Kd*G6, -1];
p6_no6 = [0, Kb, 0, -G5-Kb, G5+1i*w_p6*C, 0, -1i*w_p6*C];
p6_no7 = [0, 0, 0, 0, 0, -G6-G7, G7];
p6_no8 = [0, -G3, 0, G3+G4+G5, -1i*w_p6*C-G5, -G7, G7+1i*w_p6*C];

eq_p6 = [p6_no1;p6_no2; p6_no3; p6_no5; p6_no6; p6_no7; p6_no8];

b_p6 = [exp(-1i*pi/2); 0; 0; 0; 0; 0; 0];

res6 = eq_p6\b_p6;

v6_p6(k) = res6(5);
v6_p8(k) = res6(7);

end


T6 = v6_p6/v_in;

Tc = (v6_p6-v6_p8)/v_in;

Tin = ones(1,100)*v_in/v_in;


figure
hold on
plot (log10(f_s), 20*log10(abs(T6)));
plot (log10(f_s), 20*log10(abs(Tc)));
plot (log10(f_s), 20*log10(abs(Tin)));
ylim([-10, 10])
xlabel ("log10(f) [Hz]");
ylabel ("|T| dB");
legend('V_6','V_c','V_s')
hold off
print ("T_abs_zoom.png", "-dpng");

figure
hold on
plot (log10(f_s), 20*log10(abs(T6)));
plot (log10(f_s), 20*log10(abs(Tc)));
plot (log10(f_s), 20*log10(abs(Tin)));
xlabel ("log10(f) [Hz]");
ylabel ("|T| dB");
legend('V_6','V_c','V_s')
hold off
print ("T_abs.png", "-dpng");

figure
hold on
plot (log10(f_s), ((angle(T6))*180/pi));
plot (log10(f_s), ((angle(Tc)))*180/pi);
plot (log10(f_s), (angle(Tin)*180/pi));
xlabel ("log10(f) [Hz]");
ylabel ("Phase (degrees)");
legend('V_6','V_c','V_s')
hold off
print ("T_angle.png", "-dpng");




%% NGSPICE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Circuit 1
fileID = fopen('../sim/circ1.net','w');
   
%R1
fprintf(fileID, '* R1\nR1 1 2 %.11fk\n\n',data(1));

%R2
fprintf(fileID, '* R2\nR2 2 3 %.11fk\n\n',data(2));

%R3
fprintf(fileID, '* R3\nR3 2 5 %.11fk\n\n',data(3));

%R4
fprintf(fileID, '* R4\nR4 0 5 %.11fk\n\n',data(4));

%R5
fprintf(fileID, '* R5\nR5 5 6 %.11fk\n\n',data(5));

%R6
fprintf(fileID, '* R6\nR6 0 aldr %.11fk\n\n',data(6));

% extra voltage source to measure current
fprintf(fileID, '* Extra voltage source\nValdr aldr 7 0V\n\n');

%R7
fprintf(fileID, '* R7\nR7 7 8 %.11fk\n\n',data(7));

%Ib dependent current source
fprintf(fileID, '* Ib dependent current source\nGb 6 3 (2,5) %.11fm\n\n',data(10));

%Vc dependent supply voltage
fprintf(fileID, '* Vd dependent supply voltage\nHc 5 8 Valdr %.11fk\n\n',data(11));

% Vs supply voltage
fprintf(fileID, '* Vs supply voltage\nVs 1 0 %.11fV\n\n',data(8));

% Analysis settings
fprintf(fileID, '.control\n\nop\n\necho "********************************************"\necho  "Operating point"\necho "********************************************"\n\n\necho  "circ1_TAB"\nprint all\necho  "circ1_END"\n\nquit\n.endc\n\n.end\n');
fclose(fileID); % close circ 1


%% Circ2
fileID = fopen('../sim/circ2.net','w');
   
%R1
fprintf(fileID, '* R1\nR1 1 2 %.11fk\n\n',data(1));

%R2
fprintf(fileID, '* R2\nR2 2 3 %.11fk\n\n',data(2));

%R3
fprintf(fileID, '* R3\nR3 2 5 %.11fk\n\n',data(3));

%R4
fprintf(fileID, '* R4\nR4 0 5 %.11fk\n\n',data(4));

%R5
fprintf(fileID, '* R5\nR5 5 6 %.11fk\n\n',data(5));

%R6
fprintf(fileID, '* R6\nR6 0 aldr %.11fk\n\n',data(6));

% extra voltage source to measure current
fprintf(fileID, '* Extra voltage source\nValdr aldr 7 0V\n\n');

%R7
fprintf(fileID, '* R7\nR7 7 8 %.11fk\n\n',data(7));

%Ib dependent current source
fprintf(fileID, '* Ib dependent current source\nGb 6 3 (2,5) %.11fm\n\n',data(10));

%Vc dependent supply voltage
fprintf(fileID, '* Vd dependent supply voltage\nHc 5 8 Valdr %.11fk\n\n',data(11));

% Vs supply voltage
fprintf(fileID, '* Vs supply voltage\nVs 1 0 %.11fV\n\n',0);

% Replacing capacitor with a supply voltage
fprintf(fileID, '* Va supply voltage\nVcond 6 8 %.11fV\n\n',double(Vxfrom1));

% Analysis settings
fprintf(fileID, '.control\n\nop\n\necho "********************************************"\necho  "Operating point"\necho "********************************************"\n\n\necho  "circ2_TAB"\nprint all\necho  "circ2_END"\n\nquit\n.endc\n\n.end\n');
fclose(fileID); % close circ 2

%% Circ3
fileID = fopen('../sim/circ3.net','w');
   
%R1
fprintf(fileID, '* R1\nR1 1 2 %.11fk\n\n',data(1));

%R2
fprintf(fileID, '* R2\nR2 2 3 %.11fk\n\n',data(2));

%R3
fprintf(fileID, '* R3\nR3 2 5 %.11fk\n\n',data(3));

%R4
fprintf(fileID, '* R4\nR4 0 5 %.11fk\n\n',data(4));

%R5
fprintf(fileID, '* R5\nR5 5 6 %.11fk\n\n',data(5));

%R6
fprintf(fileID, '* R6\nR6 0 aldr %.11fk\n\n',data(6));

% extra voltage source to measure current
fprintf(fileID, '* Extra voltage source\nValdr aldr 7 0V\n\n');

%R7
fprintf(fileID, '* R7\nR7 7 8 %.11fk\n\n',data(7));

%Ib dependent current source
fprintf(fileID, '* Ib dependent current source\nGb 6 3 (2,5) %.11fm\n\n',data(10));

%Vc dependent supply voltage
fprintf(fileID, '* Vd dependent supply voltage\nHc 5 8 Valdr %.11fk\n\n',data(11));

% Vs supply voltage
fprintf(fileID, '* Vs supply voltage\nVs 1 0 %.11fV\n\n',0);

% Capacitor
fprintf(fileID, '* Capacitor\nCb 6 8 %.11fuF\n\n',data(9));

% Analysis settings
fprintf(fileID, '.end\n\n.op\n\n.ic v(6)= %d v(8)= %d\n\n.end\n\n.control\n\n',double(n_fronteira6),double(n_fronteira8));
%Changing graph colours
fprintf(fileID,'*makes plots in color\nset hcopypscolor=0\nset color0=white\nset color1=black\nset color2=red\nset color3=blue\nset color4=violet\nset color5=rgb:3/8/0\nset color6=rgb:4/0/0\n\n');
% printing out the analysis
fprintf(fileID,'echo "********************************************"\n\necho  "Transient analysis"\n\necho "********************************************"\n\ntran 1e-5 20e-3\n\nhardcopy trans3.ps v(6)\n\necho trans3_FIG\n\nquit\n\n.endc\n');
fclose(fileID); % close circ 3


%% Circ4
fileID = fopen('../sim/circ4.net','w');
   
%R1
fprintf(fileID, '* R1\nR1 1 2 %.11fk\n\n',data(1));

%R2
fprintf(fileID, '* R2\nR2 2 3 %.11fk\n\n',data(2));

%R3
fprintf(fileID, '* R3\nR3 2 5 %.11fk\n\n',data(3));

%R4
fprintf(fileID, '* R4\nR4 0 5 %.11fk\n\n',data(4));

%R5
fprintf(fileID, '* R5\nR5 5 6 %.11fk\n\n',data(5));

%R6
fprintf(fileID, '* R6\nR6 0 aldr %.11fk\n\n',data(6));

% extra voltage source to measure current
fprintf(fileID, '* Extra voltage source\nValdr aldr 7 0V\n\n');

%R7
fprintf(fileID, '* R7\nR7 7 8 %.11fk\n\n',data(7));

%Ib dependent current source
fprintf(fileID, '* Ib dependent current source\nGb 6 3 (2,5) %.11fm\n\n',data(10));

%Vc dependent supply voltage
fprintf(fileID, '* Vd dependent supply voltage\nHc 5 8 Valdr %.11fk\n\n',data(11));

% Vs supply voltage
fprintf(fileID, '* Vs supply voltage\nVs 1 0 0.0 ac 1.0 sin(0 1 1k)\n\n');

% Capacitor
fprintf(fileID, '* Capacitor\nCb 6 8 %.11fuF\n\n',data(9));

% Transient analysis settings
fprintf(fileID, '.end\n\n.op\n\n.ic v(6)= %d v(8)= %d\n\n.end\n\n.control\n\n',double(n_fronteira6),double(n_fronteira8));
%Changing graph colours
fprintf(fileID,'*makes plots in color\nset hcopypscolor=0\nset color0=white\nset color1=black\nset color2=red\nset color3=blue\nset color4=violet\nset color5=rgb:3/8/0\nset color6=rgb:4/0/0\n\n');
%Transient analysis
fprintf(fileID,'echo "********************************************"\n\necho  "Transient analysis"\n\necho "********************************************"\n\ntran 1e-5 20e-3\n\nhardcopy trans4.ps v(6) v(1)\n\necho trans4_FIG\n\n');

% Frequency analysis settings
fprintf(fileID,'echo "********************************************"\necho  "Frequency analysis"\necho "********************************************"\n\nac dec 10 0.1 1MEG\n\nhardcopy acm.ps db(v(6)) db(v(1))\necho acm_FIG\nhardcopy acp.ps 180/PI*phase(v(6)) 180/PI*phase(v(1))\necho acp_FIG\n\n');

%close file
fprintf(fileID,'quit\n\n.endc\n');

fclose(fileID); % close circ 4

