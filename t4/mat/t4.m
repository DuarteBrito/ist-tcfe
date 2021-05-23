%% gain stage

VT=25e-3;
BFN=178.7;
VAFN=69.7;
RE1=300;  %-------------
RC1=1000; %-------------
RB1=80000; %-------------
RB2=20000; %-------------
C = 80e-6; %-------------
VBEON=0.7;
VCC=12;
RS=100;

RB=1/(1/RB1+1/RB2);  % thevenin resistance
VEQ=RB2/(RB1+RB2)*VCC; % thevenin voltage
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1);  % current through the base of the first transistor
IC1=BFN*IB1; % current through the collector of the first transistor
IE1=(1+BFN)*IB1;  % current through the emmitor of the first transistor
VE1=RE1*IE1;  % voltage drop in the resistence near the emmitor of the first transistor
VO1=VCC-RC1*IC1;  %output voltage of the gain stage
VCE=VO1-VE1;  % voltage drop of the transistor

% incremental analyses
gm1=IC1/VT;
rpi1=BFN/gm1;
ro1=VAFN/IC1;

RSB=RB*RS/(RB+RS);

% AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2);
% AVI_DB = 20*log10(abs(AV1));
% AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1);
% AVIsimple_DB = 20*log10(abs(AV1simple));

% RE1=0;  % current goes through the capacitor
% AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2);  % voltage gain
% AVI_DB = 20*log10(abs(AV1));
% AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1);
% AVIsimple_DB = 20*log10(abs(AV1simple));

f = logspace(2,8,70);
w_ = 2*pi*f;
AV1_ = [];
for w=w_
RE1 = RE1/(1+1i*w*RE1*C);
AV1_ = [AV1_ , RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)];  % voltage gain
end
AVI_DB = 20*log10(abs(AV1_));
%loglog(f, AVI_DB)
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1);
AVIsimple_DB = 20*log10(abs(AV1simple));

RE1=100;
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)));
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)));
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) );
ZO1 = 1/(1/ZX+1/RC1);

RE1=0;
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)));
ZO1 = 1/(1/ro1+1/RC1);


%% ouput stage
BFP = 227.3;
VAFP = 37.2;
RE2 = 300;
VEBON = 0.7;
VI2 = VO1;
IE2 = (VCC-VEBON-VI2)/RE2;
IC2 = BFP/(BFP+1)*IE2;
VO2 = VCC - RE2*IE2;

gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/RE2;

AV2 = gm2/(gm2+gpi2+go2+ge2);
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);
ZO2 = 1/(gm2+gpi2+go2+ge2);


%% total
gB = 1/(1/gpi2+ZO1);
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1_;
AV_DB = 20*log10(abs(AV));
ZI=ZI1;
ZO=1/(go2+gm2/gpi2*gB+ge2+gB);

figure
loglog(log10(f), AV_DB)
title('gain')
xlabel('log_1_0(f) [Hz]')
ylabel('gain ')
legend('gain')
%print ("circ2.png", "-dpng");

%%
c1 = 1e-6;
c2 = 80e-6;
c3 = 35e-6;
Rpi2 =1/gpi2;
Load = 8;
k = 1;
vin = 0.01;
for t=1:0.1:8
    w = 2*pi*power(10,t);
    Zc1 = 1 ./(j*w*c1);
    Zc2 = 1 ./(j*w*c2);
    Zc3 = 1 ./(j*w*c3);
    ZRe_C1 = 1/(1/RE1+1/Zc2);
    Zeq = 1/(1/RE2+1/(Load+Zc3));
    Ro2 = 300;
    
    A = [RS+Zc1+RB,-RB,0,0,0,0,0;
        -RB,RB+rpi1+ZRe_C1, 0 , -ZRe_C1,0,0,0;
        0,rpi1*gm1,1,0,0,0,0;
        0, ZRe_C1, -ro1, ZRe_C1+ro1+RC1,-RC1,0,0;
        0,0,0,-RC1, Rpi2+RC1+Zeq,0,-Zeq;
        0,0,0,0,Rpi2*gm2,1,0;
        0,0,0,0,-Zeq,-Ro2,Zeq+Ro2];
    
    B = [vin;0;0;0;0;0;0];
    
    X=A\B;
    Vout = ((X(7))-X(5))*Zeq;
    
    gain(k) = Vout*Load/(Load+Zc3)/vin;
    %gain_DB(k)=20*log10(abs(gain(k)));
    k = k+1;
end

gain_DB=20*log10(abs(gain));
cut_off_val = max(gain_DB)-3;
[~,cut_off] = min(abs(gain_DB-cut_off_val));
%cut_off_f = [f(cut_off), 1+(cut_off-1)*(0.1)]
cut_off_f = f(cut_off)
figure
plot(1:0.1:8, gain_DB, 1:10, cut_off_val*ones(1,10))
title('gain')
xlabel('log_1_0(f) [Hz]')
ylabel('gain ')
legend('gain')
%print ("circ2.png", "-dpng");