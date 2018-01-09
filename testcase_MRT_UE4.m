% Testfall: Vorlesung MRT, Übung 4, "nachgerechnet"

clear all; % Wichtig, um Annahmen zurückzusetzen
s=LagrangeMechanicalSystem(3,3);

% Koordinaten
a=sym('q1(t)')
b=sym('q2(t)')
c=sym('c(t)') % unbenutzt
s.q=[a;b;c]

% Trägheiten
% Massen (gedanklich: im Schwerpunkt konzentriert)
s.m=[sym('m1');sym('m2');sym('m3')];
assume(s.m>0);

% sonstige nicht-geometrische Konstanten
s.g=sym('g');
assume(s.g>0)

%s.r=positionen(a,b,c)
s.r{1}=[a;0;0]
s.r{2}=[a+sym('l')*sin(b);0;-sym('l')*cos(b)]
assume(sym('l')>0)
s.r{3}=[sym(0);0;0]
s.F{1}=[sym('Fext');0;0]
s.F{2}=[0;0;0]
s.F{3}=[0;0;0]

% Drehträgheit wird hier nicht genutzt
null=[0;0;0];
Null=zeros(3);
s.phi{1}=null;
s.phi{2}=null;
s.phi{3}=null;
s.Theta{1}=null;
s.Theta{2}=null;
s.Theta{3}=null;

s.eKin()
s.ePot()

s.laplaceEquationsSimplified()
s.laplaceEquationsMupad()