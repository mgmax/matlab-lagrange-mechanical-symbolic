% Testfall: Drehung eines starren Kˆrpers mit Unwucht (Achse nicht durch
% SP)

clear all; % Wichtig, um Annahmen zur¸ckzusetzen
s=LaplaceMechanicalSystem(4,1);

% Koordinaten
phi=sym('phi(t)')
x=sym('x(t)') % unbenutzt
y=sym('y(t)') % unbenutzt
z=sym('z(t)') % unbenutzt
%s.q=[phi;x;y;z]
s.q=[phi]

% Tr‰gheiten
syms m
s.m=[sym('m/2');sym('m/2');sym('0');sym('0')];
assume(m>0);

% sonstige nicht-geometrische Konstanten
s.g=sym('g');
assume(s.g>0)

% J_sp Drehtr‰gheit um den SP
syms Jsp
% J eines Massepunktes i = r^2 * m_i
% Jsp_ges = 2* m/2 * r^2
r=sqrt(Jsp/m);
% Exzentrizit‰t
d=sym('d')
%r0=[x;y;z]
r0=[0;0;0]
null=[0;0;0];
Null=zeros(3);
s.phi{1}=null;
s.phi{2}=null;
s.phi{3}=null;
s.phi{4}=null;
s.Theta{1}=null;
s.Theta{2}=null;
s.Theta{3}=null;
s.Theta{4}=null;
s.r{1}=r0+rotMatrixInverted(phi,'z')*[d+r;0;0];
s.r{2}=r0+rotMatrixInverted(phi,'z')*[d-r;0;0];
s.r{3}=r0+rotMatrixInverted(phi,'z')*[sym('delta');0;0];
s.r{4}=r0+rotMatrixInverted(phi,'z')*[sym('-delta');0;0];
s.F{1}=[0;0;0];
s.F{2}=[0;0;0];
syms Mext delta
M=Mext%-sym('diff(phi(t),t)*kd')
% Drehmoment als Kr‰ftepaar dargestellt gem‰ﬂ Moor Skript MRT
s.F{3}=rotMatrixInverted(phi+pi/2,'z')*[M/(2*delta);0;0];
s.F{4}=rotMatrixInverted(phi-pi/2,'z')*[M/(2*delta);0;0];
s.eKin()
s.ePot()

simplify(s.laplaceEquations())
'Soll: Mext=(Jsp+m*d≤) diff(phi(t),t,t)'