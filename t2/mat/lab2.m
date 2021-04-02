close all
clear all

pkg load symbolic

%% getting initial data
data = importdata('../data.txt',"=",8);
data = data.data;
disp(data)

%% Circuit 1
fileID = fopen('../sim/circ1.net','w');
   
%R1
fprintf(fileID, '* R1\nR1 1 2 %.11fk\n',data(1));

%R2
fprintf(fileID, '* R2\nR1 2 3 %.11fk\n',data(2));

%R3
fprintf(fileID, '* R3\nR1 2 5 %.11fk\n',data(3));

%R4
fprintf(fileID, '* R4\nR1 0 5 %.11fk\n',data(4));

%R5
fprintf(fileID, '* R5\nR1 5 6 %.11fk\n',data(5));

%R6
fprintf(fileID, '* R6\nR1 0 aldr %.11fk\n',data(6));

% extra voltage source to measure current
fprintf(fileID, '* Extra voltage source\nValdr aldr 7 0V\n');

%R7
fprintf(fileID, '* R7\nR1 7 8 %.11fk\n',data(7));

%Ib dependent current source
fprintf(fileID, '* Ib dependent current source\nGb 6 3 (2,4) %.11fm\n',data(10));

%Vc dependent supply voltage
fprintf(fileID, '* Vd dependent supply voltage\nHc 5 8 Valdr %.11fk\n',data(11));

% Vs supply voltage
fprintf(fileID, '* Vs supply voltage\nHc 1 0 %.11fV\n',data(8));

% Analysis settings
fprintf(fileID, '.control\n\nop\n\necho "********************************************"\necho  "Operating point"\necho "********************************************"\n\n\necho  "op_TAB"\nprint all\necho  "op_END"\n\nquit\n.endc\n\n.end\n');
fclose(fileID); % close circ 1

%% Circ2
fileID = fopen('../sim/circ2.net','w');
   
%R1
fprintf(fileID, '* R1\nR1 1 2 %.11fk\n',data(1));

%R2
fprintf(fileID, '* R2\nR1 2 3 %.11fk\n',data(2));

%R3
fprintf(fileID, '* R3\nR1 2 5 %.11fk\n',data(3));

%R4
fprintf(fileID, '* R4\nR1 0 5 %.11fk\n',data(4));

%R5
fprintf(fileID, '* R5\nR1 5 6 %.11fk\n',data(5));

%R6
fprintf(fileID, '* R6\nR1 0 aldr %.11fk\n',data(6));

% extra voltage source to measure current
fprintf(fileID, '* Extra voltage source\nValdr aldr 7 0V\n');

%R7
fprintf(fileID, '* R7\nR1 7 8 %.11fk\n',data(7));

%Ib dependent current source
fprintf(fileID, '* Ib dependent current source\nGb 6 3 (2,4) %.11fm\n',data(10));

%Vc dependent supply voltage
fprintf(fileID, '* Vd dependent supply voltage\nHc 5 8 Valdr %.11fk\n',data(11));

% Vs supply voltage
fprintf(fileID, '* Vs supply voltage\nHc 1 0 %.11fV\n',0);

% Replacing capacitor with a supply voltage
fprintf(fileID, '* Va supply voltage\nHc 6 8 %.11fV\n',HenriqueMagicNumber);

% Analysis settings
fprintf(fileID, '.control\n\nop\n\necho "********************************************"\necho  "Operating point"\necho "********************************************"\n\n\necho  "op_TAB"\nprint all\necho  "op_END"\n\nquit\n.endc\n\n.end\n');
fclose(fileID); % close circ 2

%% Circ3
fileID = fopen('../sim/circ3.net','w');
   
%R1
fprintf(fileID, '* R1\nR1 1 2 %.11fk\n',data(1));

%R2
fprintf(fileID, '* R2\nR1 2 3 %.11fk\n',data(2));

%R3
fprintf(fileID, '* R3\nR1 2 5 %.11fk\n',data(3));

%R4
fprintf(fileID, '* R4\nR1 0 5 %.11fk\n',data(4));

%R5
fprintf(fileID, '* R5\nR1 5 6 %.11fk\n',data(5));

%R6
fprintf(fileID, '* R6\nR1 0 aldr %.11fk\n',data(6));

% extra voltage source to measure current
fprintf(fileID, '* Extra voltage source\nValdr aldr 7 0V\n');

%R7
fprintf(fileID, '* R7\nR1 7 8 %.11fk\n',data(7));

%Ib dependent current source
fprintf(fileID, '* Ib dependent current source\nGb 6 3 (2,4) %.11fm\n',data(10));

%Vc dependent supply voltage
fprintf(fileID, '* Vd dependent supply voltage\nHc 5 8 Valdr %.11fk\n',data(11));

% Vs supply voltage
fprintf(fileID, '* Vs supply voltage\nHc 1 0 %.11fV\n',0);

% Capacitor
fprintf(fileID, '* Capacitor\nCb 6 8 %.11%fuF\n',data(9));

% Analysis settings
fprintf(fileID, 'echo "********************************************"\necho  "Transient analysis"\necho "********************************************"\n.ic v(6)= %d v(8)= %d\ntran 1e-5 20e-3\n\nhardcopy trans.ps v(V6) v(base)\necho trans_FIG\n\nquit\n.endc\n\n.end\n',numero2,numero3);
fclose(fileID); % close circ 3


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
syms Va
syms Kc
syms Kb
syms Id

Z = vpa(0.0);
U = vpa(1.0);
% An = [U, Z,Z,Z,Z,Z,Z;
%     Z, Z,Z,1,Z,Kc*G6, -U;
%     -G1,G1+G2+G3, -G2,-G3,Z,Z,Z;
%     Z,Kb+G2,-G2, -Kb, Z,Z,Z; 
%     Z,-Kb,Z,G5+Kb,-G5, Z, Z;
%     Z,Z,Z,Z,Z,G6-G7,G7;
%     -G1,G1,Z,G4, Z, G6, Z];
An = [U, Z,Z,Z,Z,Z,Z;              %1
%     Z, -G3,Z,G4+G5+G3,-G5,-G7,G7;  %4-7
    Z,Z,Z,U,Z, Kc*G6,-U;
    G1,-G1-G2-G3, G2,G3,Z,Z,Z;     %2
    Z,-Kb-G2,G2, Kb, Z,Z,Z;        %3
    Z,-Kb,Z,G5+Kb,-G5, Z, Z;       %5
    Z,Z,Z,Z,Z,-G6-G7,G7;           %6
    G1,-G1,Z,-G4, Z, -G6, Z];      %0-1
Bn=[Va;Z;Z ;Z; -Id; Z; Z];

Am=[R1+R3+R4 , R3, +R4, Z;
  Kb*R3 , Kb*R3-U, Z, Z ;
  R4, 0 , R6-Kc+R7+R4 , Z;
  Z,Z,Z,U];
Bm = [Va; Z; Z; Id];
fn = An\Bn;
fm = Am\Bm;

res_malhas = subs(fm,{R1, R2, R3, R4,R5, R6, R7,Va, Kc, Kb, Id},{1.04944227714,2.06296295698, 3.07855037163, 4.04814283444, 3.03583837907, 2.01824745844, 1.04357678508, 5.07638677695, 8.10223845988, 7.26693007101, 1.0053213836});

res_nos = subs(fn,{G1,G2, G3, G4, G5, G6, G7, Va, Kc, Kb, Id},{1/1.04944227714,1/2.06296295698, 1/3.07855037163, 1/4.04814283444, 1/3.03583837907, 1/2.01824745844, 1/1.04357678508, 5.07638677695, 8.10223845988, 7.26693007101, 1.0053213836});

res_nos = double(res_nos);
res_malhas = double(res_malhas);

i1 = res_malhas(1);
i2 = res_malhas(2);
i3 = res_malhas(1)+res_malhas(2);
i4 = res_malhas(1)+res_malhas(3);
i5 = res_malhas(2)-res_malhas(4);
i6 = res_malhas(3);
i7 = res_malhas(3);

disp([i1; i2; i3; i4; i5; i6; i7])

fidMalhas = fopen("malhas.txt","w");
fprintf(fidMalhas," ,I(mA)\n");
fprintf(fidMalhas,"Ia,%f\n",res_malhas(1));
fprintf(fidMalhas,"Ib,%f\n",res_malhas(2));
fprintf(fidMalhas,"Ic,%f\n",res_malhas(3));
fprintf(fidMalhas,"Id,%f\n",res_malhas(4));
fclose(fidMalhas);

fidNos = fopen("nos.txt","w");
fprintf(fidNos," ,V(V)\n");
fprintf(fidNos,"V1,%f\n",res_nos(1));
fprintf(fidNos,"V2,%f\n",res_nos(2));
fprintf(fidNos,"V3,%f\n",res_nos(3));
fprintf(fidNos,"V4,%f\n",res_nos(4));
fprintf(fidNos,"V5,%f\n",res_nos(5));
fprintf(fidNos,"V6,%f\n",res_nos(6));
fprintf(fidNos,"V7,%f\n",res_nos(7));
fclose(fidNos);

fidCur = fopen("cur.txt","w");
fprintf(fidCur," ,I(mA)\n");
fprintf(fidCur,"I1,%f\n",i1);
fprintf(fidCur,"I2,%f\n",i2);
fprintf(fidCur,"I3,%f\n",i3);
fprintf(fidCur,"I4,%f\n",i4);
fprintf(fidCur,"I5,%f\n",i5);
fprintf(fidCur,"I6,%f\n",i6);
fprintf(fidCur,"I7,%f\n",i7);
fclose(fidCur);

is = [i1; i2; i3; i4; i5; i6; i7];
rs=[1.04944227714;2.06296295698;3.07855037163;4.04814283444;3.03583837907;2.01824745844;1.04357678508];
vs = rs.*is;
vs_nos = [res_nos(1)-res_nos(2);
    res_nos(3)-res_nos(2);
    res_nos(2)-res_nos(4);
    res_nos(4)-0;
    res_nos(4)-res_nos(5);
    0-res_nos(6);
    res_nos(6)-res_nos(7)];


fidComp = fopen("comp.txt","w");
fprintf(fidComp," n R , I(mA), RI(V), V nodal(V), Difference(V)\n");
for k=1:length(rs)
fprintf(fidComp,"%f,%f, %f, %f,%f\n",fix(k),is(k), vs(k), vs_nos(k), abs(vs(k)- vs_nos(k)) );
end
fclose(fidComp);
