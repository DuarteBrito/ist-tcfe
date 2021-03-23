close all
clear all

pkg load symbolic

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

fidMalhas = fopen("malhas.csv","w");
fprintf(fidMalhas," ,I(mA);");
fprintf(fidMalhas,"Ia,%f\n;",res_malhas(1));
fprintf(fidMalhas,"Ib,%f\n;",res_malhas(2));
fprintf(fidMalhas,"Ic,%f\n;",res_malhas(3));
fprintf(fidMalhas,"Id,%f\n;",res_malhas(4));
fclose(fidMalhas);

fidNos = fopen("nos.csv","w");
fprintf(fidNos," ,V(V)");
fprintf(fidNos,"V1,%f\n",res_nos(1));
fprintf(fidNos,"V2,%f\n",res_nos(2));
fprintf(fidNos,"V3,%f\n",res_nos(3));
fprintf(fidNos,"V4,%f\n",res_nos(4));
fprintf(fidNos,"V5,%f\n",res_nos(5));
fprintf(fidNos,"V6,%f\n",res_nos(6));
fprintf(fidNos,"V7,%f\n",res_nos(7));
fclose(fidNos);

fidCur = fopen("cur.csv","w");
fprintf(fidCur," ,I(mA)");
fprintf(fidCur,"I1,%f\n",i1);
fprintf(fidCur,"I2,%f\n",i2);
fprintf(fidCur,"I3,%f\n",i3);
fprintf(fidCur,"I4,%f\n",i4);
fprintf(fidCur,"I5,%f\n",i5);
fprintf(fidCur,"I6,%f\n",i6);
fprintf(fidCur,"I7,%f\n",i7);
fclose(fidCur);
