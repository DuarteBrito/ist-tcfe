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
A = [U, Z,Z,Z,Z,Z,Z;
    Z, Z,Z,1,Z,Kc*G6, -U;
    -G1,G1+G2+G3, -G2,-G3,Z,Z,Z;
    Z,Kb+G2,-G2, -Kb, Z,Z,Z; 
    Z,-Kb,Z,G5+Kb,-G5, Z, Z;
    Z,Z,Z,Z,Z,G6-G7,G7;
    -G1,G1,Z,G4, Z, G6, Z];
B=[Va;Z;Kb ;Z; -Id; Z; Z];
% A=[R1+R3-R4 , R3, -R4, Z;
%   Kb*R3 , Kb*R3-1, Z, Z ;
%   -R4, 0 , R6-Kc+R7-R4 , Z;
%   Z,Z,Z,U];
% B = [Va; Z; Z; Id];
f = A\B;
%res = subs(f,{R1, R2, R3, R4,R5, R6, R7,Va, Kc, Kb, Id},{1.04944227714,2.06296295698, 3.07855037163, 4.04814283444, 3.03583837907, 2.01824745844, 1.04357678508, 5.07638677695, 8.10223845988, 7.26693007101, 1.0053213836});

res = subs(f,{G1,G2, G3, G4, G5, G6, G7, Va, Kc, Kb, Id},{1/1.04944227714,1/2.06296295698, 1/3.07855037163, 1/4.04814283444, 1/3.03583837907, 1/2.01824745844, 1/1.04357678508, 5.07638677695, 8.10223845988, 7.26693007101, 1.0053213836})

res = double(res)

