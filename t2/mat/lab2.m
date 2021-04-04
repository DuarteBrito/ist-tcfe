close all
clear all

pkg load symbolic

%% getting initial data
data = importdata('../data.txt',"=",8);
%data = importdata('data.txt',"=",8);
data = data.data;
%disp(data)

%% MATLAB   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% chamar variáveis

%pkg load symbolic

syms G1
syms G2
syms G3
syms G4
syms G5
syms G6
syms G7
syms R1
syms R2
syms R3
syms R4
syms R5
syms R6
syms R7
syms Vs
syms Kb
syms Kd
syms Id

Z = vpa(0.0);
U = vpa(1.0);

data = [1.04944227714 , 2.06296295698, 3.07855037163, 4.04814283444, 3.03583837907, 2.01824745844, 1.04357678508,  5.07638677695, 1.0053213836, 7.26693007101, 8.10223845988];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% circuit 1
    % V1 V2 V3 V5 V6 V7 V8
A1 = [U ,Z,Z,Z,Z,Z,Z;
    G1,-G1-G2-G3, G2,G3, Z,Z,Z; 
    Z,+Kb+G2,-G2, -Kb, Z, Z,Z; %
    Z, Z,Z,U,Z,Kd*G6, -U;
    Z,-Kb,Z,G5+Kb,-G5, Z, Z;
    Z,Z,Z,Z,Z,-G6-G7,G7;
    G1,-G1,Z,-G4,Z,-G6,Z];
    %Z,G3,-G3,-G4-G5,G5,G7,-G7];
 
B1 = [Vs;Z;Z;Z;Z;Z;Z];

f1 = A1\B1;
res1 = subs(f1,{G1, G2, G3, G4,G5, G6, G7,Vs, Kb, Kd},{1/data(1), 1/data(2), 1/data(3), 1/data(4),1/data(5), 1/data(6), 1/data(7), data(8), data(10), data(11)});
Vxfrom1 = res1(5)-res1(7);

%% circuit 2
syms Vx
    %V2 V3 V5 V6 V7 V8
A2 = [-G1-G2-G3, G2,G3, Z,Z,Z; 
    +Kb+G2,-G2, -Kb, Z, Z,Z; %
    Z,Z,U,Z,Kd*G6, -U;
    Z,Z,Z,U,Z,-U;
    Z,Z,Z,Z, -G6-G7, G7;
    G3-Kb,Z,-G4-G3-Kb,Z,G7,-G7];
    %-G1,Z,-G4,Z,-G6,Z];

B2 = [Z;Z;Z;Vx;Z;Z]; %VX
f2 = A2\B2;
res2 = subs(f2,{G1, G2, G3, G4,G5, G6, G7,Vs, Kb, Kd, Vx},{1/data(1), 1/data(2), 1/data(3), 1/data(4),1/data(5), 1/data(6), 1/data(7), 0, data(10), data(11), Vxfrom1}); 

res2 = double(res2);

Ix = (res2(4)-res2(3))/data(5) + data(10)*(res2(1)-res2(3));
R = (res2(4)-res2(6))/Ix;

t = 0:0.1e-5:20e-3;
V_6n = (res2(4)-res2(6))*exp((-1/(R * data(9)*1e-3))*t);

n_fronteira6 = res2(4);
n_fronteira8 = res2(6);

figure
plot(t,V_6n)

title('Natural Solution')
xlabel('t [s]')
ylabel('V_6_n [V]')
print ("circ2.png", "-dpng");

%% circuit 3
syms w
syms C
Im = vpa(1i);
freq = 1e3;
W = 2*pi*freq;
% circuito 4
    % V~1 V~2 V~3 V~5 V~6 V~7 V~8
A3 = [U ,Z,Z,Z,Z,Z,Z;
    G1,-G1-G2-G3, G2,G3, Z,Z,Z; 
    Z,+Kb+G2,-G2, -Kb, Z, Z,Z; %
    Z, Z,Z,U,Z,Kd*G6, -U;
    Z,-Kb,Z,G5+Kb,-G5-Im*w*C, Z, Im*w*C; %
    Z,Z,Z,Z,Z,-G6-G7,-G7;
    Z,G3,-G3,-G4-G5,G5+Im*w*C,G7,-G7-Im*w*C];
 
B3 = [Vs;Z;Z;Z;Z;Z;Z];

f3 = A3\B3;
res3 = subs(f3,{G1, G2, G3, G4,G5, G6, G7,Vs, C, Kb, Kd, w},{1/data(1), 1/data(2), 1/data(3), 1/data(4),1/data(5), 1/data(6), 1/data(7), exp(-1i*(pi/2)),data(9)*1e-6, data(10), data(11), W});

res3 = double(res3);

V_circ3 = abs(res3);
angle_circ3 = angle(res3);

t = 0:0.1e-5:20e-3;
% W*t 
% angle_circ3(5)
% exp(1i*(W*t-  angle_circ3(5)));
%exp(1i*(W*t - angle_circ3(5)))
V_6f = real(V_circ3(5) * exp(1i*(W*t - angle_circ3(5)))); 


figure
plot(t,real(V_6f))
title('Forced Solution')
xlabel('t [s]')
ylabel('V_6_f [V]')
print ("circ3.png", "-dpng");


%% ponto 5
freq = 1e3;

t_ = -5e-3:0.1e-5:0-0.1e-5;
vs_menor = zeros(1,length(t_));
v_6_menor = zeros(1,length(t_));

for n=1:length(t_)
    vs_menor(n) = data(8);
    v_6_menor(n) = double(res1(5));
end

t = 0:0.1e-5:20e-3;
vs_maior = sin(2*pi*freq*t);
v_6_maior = V_6n + V_6f;

t = [t_ t];
vs = [vs_menor vs_maior];
v_6 = [v_6_menor v_6_maior];
    
figure
plot(t,vs, t,v_6)
title('Final Solution')
xlabel('t [s]')
ylabel('V [V]')
print ("tot_ponto5.png", "-dpng");



%% ponto 6
V_in = exp(-1i *pi/2);
trans_w = subs(f3,{G1, G2, G3, G4,G5, G6, G7,Vs, C, Kb, Kd},{1/data(1), 1/data(2), 1/data(3), 1/data(4),1/data(5), 1/data(6), 1/data(7), exp(-1i*(pi/2)),data(9)*1e-6, data(10), data(11)});
trans_w_6 = trans_w(5);
trans_w_8 = trans_w(7);
w_s=0.1:49999.995:1e6;
T1 = zeros(1, length(w_s));
T2 = T1;
for n=1:length(w_s)
    T1(n) = double(subs(trans_w_6,{w},{w_s(n)}))/V_in;
    T2(n) = double(subs(trans_w_6,{w},{w_s(n)})-subs(trans_w_8,{w},{w_s(n)}))/V_in;
end
%T = j*w*L ./(R + j*w*L)

figure
plot (log10(w_s/2/pi), 20*log10(abs(T1)), 20*log10(abs(T2)), 20*log10(abs(V_in)));
% plot (log10(w_s/2/pi), 20*log10(abs(T2)));
% plot (log10(w_s/2/pi), 20*log10(abs(V_in)));
xlabel ("log10(w) [rad/s]");
ylabel ("|T| dB");
print ("T_abs.png", "-dpng");

figure
hold on
plot (log10(w_s/2/pi), 20*log10(angle(T1)*180/pi));
plot (log10(w_s/2/pi), 20*log10(angle(T2)*180/pi));
plot (log10(w_s/2/pi), 20*log10(angle(V_in)*180/pi));
xlabel ("log10(w) [rad/s]");
ylabel ("Phase (degrees)");
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
fprintf(fileID, '.end\n\n.op\n\n.ic v(6)= %d v(8)= %d\n\n.end\n\n.control\n\necho "********************************************"\n\necho  "Transient analysis"\n\necho "********************************************"\n\ntran 1e-5 20e-3\n\nhardcopy trans3.ps v(6)\n\necho trans3_FIG\n\nquit\n\n.endc\n',double(n_fronteira6),double(n_fronteira8));
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
fprintf(fileID, '.end\n\n.op\n\n.ic v(6)= %d v(8)= %d\n\n.end\n\n.control\n\necho "********************************************"\n\necho  "Transient analysis"\n\necho "********************************************"\n\ntran 1e-5 20e-3\n\nhardcopy trans4.ps v(6) v(1)\n\necho trans4_FIG\n\n',double(n_fronteira6),double(n_fronteira8));

% Frequency analysis settings
fprintf(fileID,'echo "********************************************"\necho  "Frequency analysis"\necho "********************************************"\n\nac dec 10 0.1 1MEG\n\nhardcopy acm.ps db(v(6)) v(1)\necho acm_FIG\nhardcopy acp.ps v(1)\necho acp_FIG\n\n');

%close file
fprintf(fileID,'quit\n\n.endc\n');

fclose(fileID); % close circ 4

