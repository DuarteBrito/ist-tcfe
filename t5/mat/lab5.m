% values

Rb1 = 1e3;
Rb2 = 1e3;


Rb3 = 110e3;
Rb4 = 1e3;

R_cost = (Rb1 + Rb2 + Rb3 + Rb4)/1000;

Cb1 = 0.22e-6;

Cb2 = 0.11e-6;

C_cost = (Cb1 + Cb2*4) * 1e6;



% point 1
% calculation for central frequency

w_L = 1/(Rb1*Cb1);
f_L = w_L/2/pi;
w_H = 1/(Rb2*Cb2);
f_H = w_H/2/pi;

w_0 = sqrt(w_L*w_H);
f_0 = w_0/2/pi;

deviation = f_0-1000;

gain_1 = (j*w_0*Rb1*Cb1)/(1+j*w_0*Rb1*Cb1);
gain_2 = 1+Rb3/Rb4;
gain_3 = 1/(1+1i*w_0*Rb2*Cb2);

gain = 20*log10(abs(gain_1 * gain_2 * gain_3));
Z_i = (1/(j*w_0*Cb1) + Rb1);
Z_o = (1/(1/Rb2+j*w_0*Cb2)); 
 

% point 2
% frequency response

f = logspace(1,8,70);
w = 2*pi*f;

f_res = ((Rb1*Cb1.*w*j)./(1+Rb1*Cb1.*w*j)).*gain_2.*(1./(1+Rb2*Cb2.*w*j));

maxi = max(abs(f_res));

figure
semilogx(f,20*log10(abs(f_res)))

xlabel("Frequency [Hz]");
ylabel("Gain [dB]");
title("Gain");
print("gain.png", "-dpng");


figure
semilogx(f,angle(f_res))
xlabel("Frequency [Hz]");
ylabel("Phase [Rad]");
title("Phase");
print("phase.png", "-dpng");

% merit
cost = 1.3323e4 + R_cost + C_cost;

merit = 1/(cost*(abs(maxi-100)+abs(deviation)+1e-6));

% resultados
fidRes = fopen("resultados.txt","w");
fprintf(fidRes,"Theoretical results ,value\n");
fprintf(fidRes,"gain ,%f (dB)\n",gain);
fprintf(fidRes,"input impedance,%f Ohm\n",Z_i);
fprintf(fidRes,"output impedance,%f Ohm\n",Z_o);
fclose(fidRes);



% resultados
fidRes2 = fopen("resultados2.txt","w");
fprintf(fidRes2,"Theoretical results ,value (Hz)\n");
fprintf(fidRes2,"$f_L$,%f\n",f_L);
fprintf(fidRes2,"$f_H$,%f\n",f_H);
fprintf(fidRes2,"$f_0$,%f\n",f_0);
fclose(fidRes2);


% merito
fidRes3 = fopen("merito.txt","w");
fprintf(fidRes3,"M & %d \\\\ \n \\hline\n",merit);


