% Bestimmung der Gleichung w=(...) im Abschnitt "vektorielle
% Winkelgeschwindigkeit"
% vektorielle WInkelgeschw. abhängig von phiX=b, phiZ=c
b=sym('phi_x(t)');
c=sym('phi_z(t)');
T=sym(simplify(sym(rotMatrixInverted(c,'z')^(-1) * rotMatrixInverted(b,'x')^-1)));
A=simplify(transpose(diff(T,'t'))*T)
omega=simplify([A(3,2);A(1,3);A(2,1)])