close all;
clear all; % Wichtig, um Annahmen zurückzusetzen
s=LagrangeMechanicalSystem(7,3);

% Koordinaten
a=sym('a(t)');
b=sym('b(t)');
c=sym('c(t)');
s.q=[a;b;c];

% sonstige nicht-geometrische Konstanten
s.g=sym('g');
assume(s.g>0)

% Linearschiene
syms mA
s.m(1)=mA;
F1=sym('-dA*diff(a(t),t)-cA*sign(diff(a(t),t))+kA*iA(t)');
s.F{1}=[F1;0;0];
s.r{1}=[a;0;0];
s.phi{1}=[0;0];
s.Theta{1}=zeros(3);


% Neigeeinheit
% 2: linker Hebelpunkt
% 3: rechter Hebelpunkt
% 4: Masse

% u: Unwucht = Abstand SP-Drehachse
% verschiebung: konstante Verschiebung (ist egal)
% delta: gedachte Hebellänge (kürzt sich raus)

syms mB delta2 verschiebung2xEgal verschiebung2yEgal verschiebung2zEgal
ursprung2=s.r{1}+[verschiebung2xEgal;verschiebung2yEgal;verschiebung2zEgal];
assume(mB>0 & delta2>0)
M2=sym('-dB*diff(b(t),t)-cB*sign(diff(b(t),t))+kB*iB(t)');
Fhebel2=M2/delta2/2;

% "Hebelpunkte" (Angriffspunkte für Kräftepaar, masselos)
s.m(2)=0;
s.r{2}=ursprung2+rotMatrixInverted(b,'x')*[0;delta2;0];
s.F{2}=rotMatrixInverted(b,'x')*[0;0;Fhebel2];
s.Theta{2}=zeros(3);
s.phi{2}=[0;0];

s.m(3)=0;
s.r{3}=ursprung2+rotMatrixInverted(b,'x')*[0;-delta2;0];
s.F{3}=rotMatrixInverted(b,'x')*[0;0;-Fhebel2];
s.Theta{3}=zeros(3);
s.phi{3}=[0;0];

% "Massepunkt" B
s.m(4)=mB;
syms uBxEgal uBy uBz;
s.r{4}=ursprung2+rotMatrixInverted(b,'x')*[uBxEgal;uBy;uBz];
s.F{4}=[0;0;0];
syms TBx TByEgal TBzEgal TBxyEgal TBxzEgal TByzEgal
s.Theta{4}=[TBx TBxyEgal TBxzEgal; TBxyEgal TByEgal TByzEgal; TBxzEgal TByzEgal TBzEgal];
s.phi{4}=[b;0];

% Drehteller
% 5: linker Hebelpunkt
% 6: rechter Hebelpunkt
% 7: Masse

% u: Unwucht = Abstand SP-Drehachse
% uCx/y/z: Abweichung SP Dreheinheit - Drehachse Neigeeinheit
% ux3: x-Abw. ...., wenn Dreharm in x-Richtung steht
% delta: gedachte Hebellänge (kürzt sich raus)
syms mC uCx uCy uCz delta3
syms verschiebung3xEgal
assume(mC>0 & delta3>0)
M3=sym('-dC*diff(c(t),t)-cC*sign(diff(c(t),t))+kC*iC(t)');
Fhebel3=M3/delta3/2;
ursprung3=ursprung2+[verschiebung3xEgal;0;0];
% Angriffspunkte für Kräftepaar
s.m(5)=0;
s.r{5}=ursprung3+rotMatrixInverted(b,'x')*rotMatrixInverted(c,'z')*[delta3;0;0];
s.F{5}=rotMatrixInverted(b,'x')*rotMatrixInverted(sym('c(t)'),'z')*[0;Fhebel3;0];
s.phi{5}=[0;0];
s.Theta{5}=zeros(3);
s.m(6)=0;
s.r{6}=ursprung3+rotMatrixInverted(b,'x')*rotMatrixInverted(c,'z')*[-delta3;0;0];
s.F{6}=rotMatrixInverted(b,'x')*rotMatrixInverted(sym('c(t)'),'z')*[0;-Fhebel3;0];
s.phi{6}=[0;0];
s.Theta{6}=zeros(3);

% "Massepunkt"
s.m(7)=mC;
s.r{7}=ursprung3+rotMatrixInverted(b,'x')*rotMatrixInverted(c,'z')*[uCx;uCy;uCz];
s.F{7}=[0;0;0];
s.phi{7}=[b;c];
% Einträge von Theta
syms TCx TCy TCz TCxy TCxz TCyz
s.Theta{7}=[TCx TCxy TCxz; TCxy TCy TCyz; TCxz TCyz TCz];
nb=s.lagrangeEquationsMupad();

% s.eKin()
% s.ePot()
% g=s.lagrangeEquationsSimplified();
texVariables={
'mA','m_A';
'mB','m_B';
'mC','m_C';
'mGes','m_{ges}';
'kA','k_A';
'kB','k_B';
'kC','k_C';
'dA','d_A';
'dB','d_B';
'dC','d_C';
'cA','c_A';
'cB','c_B';
'cC','c_C';
'iA','i_A';
'iB','i_B';
'iC','i_C';
'uCx','u_{C,x}';
'uCy','u_{C,y}';
'uCz','u_{C,z}';
'uBy','u_{B,y}';
'uBz','u_{B,z}';
'TBx','{\Theta_{B,x}}';
'TCx','{\Theta_{C,x}}';
'TCy','{\Theta_{C,y}}';
'TCz','{\Theta_{C,z}}';
'TCxy','{\Theta_{C,xy}}';
'TCxz','{\Theta_{C,xz}}';
'TCyz','{\Theta_{C,yz}}';
'JC','{J_C}';
};
texReplacements={'\frac{\partial}{\partial t} a\!\left(t\right)','\dot a(t)';
    '\frac{\partial}{\partial t} b\!\left(t\right)','\dot b(t)';
    '\frac{\partial}{\partial t} c\!\left(t\right)','\dot c(t)';
    '\frac{\partial^2}{\partial t^2} a\!\left(t\right)','\ddot a(t)';
    '\frac{\partial^2}{\partial t^2} b\!\left(t\right)','\ddot b(t)';
    '\frac{\partial^2}{\partial t^2} c\!\left(t\right)','\ddot c(t)';
    'sign','\sign'; % dazu \DeclareMathOperator{\sign}{sign} in der Präambel des Dokuments
    };
t=s.lagrangeEquationsTex(texVariables,texReplacements);

f=fopen('tex.txt','w');
fprintf(f,'%s',evalc('t'));
fclose(f);

g=s.lagrangeEquationsSimplified()

assume(b==0);
assume(uCy==0);
syms mGes JC;
assume(mA==mGes-mB-mC);
assume(TCz==JC-uCx^2*mC);

tB=s.equationsToTex(simplify(g),texVariables,texReplacements);
f=fopen('texText_Darmstadtaehnlich.txt','w');
fprintf(f,'%s',evalc('tB'));
fclose(f);
% TODO wenn man jetzt noch andere Fälle erzeugen will, muss man mit sym mA b uCy <...> clear die Annahmen entfernen...