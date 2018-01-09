classdef LagrangeMechanicalSystem < handle
    % mechanisches Massepunktsystem nach Lagrange
    %   Detailed explanation goes here
    
    properties
        m=[] % Massen
        Theta={} % Trägheitstensor
        r={} % Punkte [x;y;z] abhängig von q1...qn
        phi={} % Winkel [PhiX;PhiZ] abhängig von q1...qn
        F={}
        q=[] % symbolische Variablen [sym('q1');sym('q2');...]
        n=0
        nFrei=0
        g=9.81
    end
    methods
        function this=LagrangeMechanicalSystem(n,nFrei)
            % mechanisches System mit n Punkten und nFrei Freiheitsgraden
            this=this@handle;
            this.n=n;
            this.nFrei=nFrei;
            assert(nFrei>0);
            assert(nFrei<5*n); % ein System mit 5n F.graden wäre unbestimmt und nutzlos - 5n, weil hier der 6te (phiX) immer 0 sein muss
            this.m=sym(ones(n,1)*NaN);
            this.r=cell(n, 1);
            this.phi=cell(n, 1);
            this.F=cell(n, 1);
            this.q=sym(ones(nFrei,1)*NaN);
        end
        
        
        
        function e=eKin(this)
            % kinetische Energie = 0.5 m*v^2
            e=0;
            for i=1:this.n
                e=e+0.5*this.m(i)*absSquare(diff(sym(this.r{i}),'t'));
            end
            e=simplify(e);
        end
        
        function omega=omega(this,b,c)
            % vektorielle WInkelgeschw. abhängig von phiX=b, phiZ=c
            % transformiert ins körperfeste Koordinatensystem
            % also eigentlich T*omega
            T=sym(simplify(sym(rotMatrixInverted(c,'z')^(-1) * rotMatrixInverted(b,'x')^-1)));
            A=transpose(diff(T,'t'))*T;
            omega=simplify([A(3,2);A(1,3);A(2,1)]);
            omega=T*omega
        end
        
        function e=eRot(this)
             % Rotationsenergie = 0.5 omega transponiert * Theta * omega
            e=0;
            for i=1:this.n
                omega_i=this.omega(this.phi{i}(1),this.phi{i}(2));
                e=e+0.5*transpose(omega_i)*this.Theta{i}*omega_i;
            end
            e=simplify(e);
        end
        
        function e=ePot(this)
            % potentielle Energie = m*g*h (h=z-Koordinate)
            e=0;
            for i=1:this.n
                e=e+this.m(i)*this.g*this.r{i}(3);
            end
            e=simplify(e);
        end
        
        function leftSide=lagrangeEquationLeftside(this,j)
            assert(0<j)
            assert(j<=this.nFrei)
            % TODO qj string
            % L nach (dqj/dt) ableiten
           L=this.eKin()+this.eRot()-this.ePot();
           temp=partialDiffSpecial(L,this.q(j),1);
           leftSide=diff(sym(temp),'t') - partialDiffSpecial(L,this.q(j),0);
        end
        
        function Qj=lagrangeEquationRightside(this,j)
             assert(0<j)
            assert(j<=this.nFrei)
            % verallgemeinerte Kraft Qj
            Qj=0;
            for i=1:this.n
                % Qj = summe( Fi transponiert * dri/dqj )
                Qj=Qj+transpose(this.F{i})*partialDiffSpecial(this.r{i},this.q(j),0);
            end
        end
        
        function y=lagrangeEquations(this)
            assert(isequal(size(this.q),[this.nFrei,1]));
            assert(isequal(size(this.r),[this.n,1]));
            assert(isequal(size(this.F),[this.n,1]));
            assert(isequal(size(this.m),[this.n,1]));
            y=sym(zeros(this.nFrei,1));
            for j=1:this.nFrei
                y(j)=sym(this.lagrangeEquationLeftside(j))==sym(this.lagrangeEquationRightside(j));
            end
        end
        
        function y=lagrangeEquationsSimplified(this)
            
            y=simplify(this.lagrangeEquations(),100);
            y=simple(y);
            y=simple(y);
            y=simple(y);
            for i=1:this.nFrei
                % umstellen auf Form "... = 0"
                % und ausklammern von ersten/zweiten ABleitungen
                try
                    y(i)=y(i)-evalin(symengine,['rhs(',char(y(i)),')']);
                    y(i)=evalin(symengine,['collect(',char(y(i)),',[diff(b(t),t),diff(b(t),t,t),diff(a(t),t),diff(a(t),t,t),diff(c(t),t),diff(c(t),t,t)])']);
                    y(i)=y(i)-evalin(symengine,['rhs(',char(y(i)),')']);
                    y(i)=evalin(symengine,['Simplify(',char(y(i)),',Steps=400,Valuation=length)']);
                    y(i)=y(i)-evalin(symengine,['rhs(',char(y(i)),')']);
                catch
                    warning 'Could not get RHS of expression or collect terms - missing assume(x>0) ? Expression:'
                    char(y(i))
                end
            end
        end
        function y=lagrangeEquationsTex(this,variablesTex,replacementsTex)
            % erzeuge Tex Output zu LagrangeGleichungen
            y=this.equationsToTex(this.lagrangeEquationsSimplified(),variablesTex,replacementsTex);
        end     
        function t=equationsToTex(this,eqn,variablesTex,replacementsTex)
           % erzeuge Tex Output zu gegebenen Gleichungen
           % tex-Variablennamen stehen dabei in variablesTex, z.B.:
           % variablesTex= {'a1','\alpha_1';'b','\beta'}
           % weitere Ersetzungen z.B. für Zeitableitungen in
           % replacementsTex im gleichen FOrmat. Wenn unbenutzt:
           % replacementsTex={'',''}
           % ACHTUNG! replacementsTex achtet nicht auf Befehlsgrenzen o.ä.,
           % also dort bitte die Ersetzungen als '{ersatzFoo}' passend
           % klammern
           gln=this.lagrangeEquationsSimplified();
           lhs=sym(zeros(this.nFrei,1));
           rhs=sym(zeros(this.nFrei,1));
           for i=1:this.nFrei
                % umstellen auf Form "... = 0"
                try
                    lhs(i)=evalin(symengine,['lhs(',char(gln(i)),')']);
                    rhs(i)=evalin(symengine,['rhs(',char(gln(i)),')']);
                catch
                    warning 'could not get lhs/rhs of expression:'
                    char(gln(i))
                    rhs(i)=sym(NaN);
                    lhs(i)=sym(NaN);
                end
           end
           y=cell(this.nFrei,1);
           s=size(variablesTex);
           assert(s(2)==2)
           for i=1:this.nFrei
               y{i}=['0=' latex(lhs(i)-rhs(i))];
               
               for j=1:s(1)
                    y{i}=strrep(y{i},['{' variablesTex{j,1} '}'],['{' variablesTex{j,2} '}']);
               end
               y{i}=strrep(y{i},'\mathrm','');
%                y{i}=regexprep(y{i},'\\mathrm{([^\{]*)}',);
           end
           s=size(replacementsTex);
           assert(s(2)==2)
           for i=1:this.nFrei
               for j=1:s(1)
                    y{i}=strrep(y{i},replacementsTex{j,1},replacementsTex{j,2});
               end
%                y{i}=regexprep(y{i},'\\mathrm{([^\{]*)}',);
           end
           newline=sprintf('\n');
            start=['\begin{dmath}' newline];
            ende=[newline '\end{dmath}' newline];
            sep=[ende start];
            t=['\begin{dgroup}' newline start y{1}];
            for i=2:this.nFrei
                t=[t sep y{i}];
            end
            t=[t ende  '\end{dgroup}'];
        end
        
        
        
        function nb=lagrangeEquationsMupad(this)
           d=this.lagrangeEquationsSimplified();
           q=sym(this.q);
           r=sym(this.r);
           f=this.F;
           for i=1:this.n
              f{i}=sym(f{i});
              
           end
           f=simple(simplify(sym(f)));
           ekin=this.eKin();
           erot=this.eRot();
           epot=this.ePot();
           nb=mupad('mupad.mn');
           setVar(nb,d)
           setVar(nb,r)
           setVar(nb,f)
           setVar(nb,q)
           setVar(nb,ekin)
           setVar(nb,erot)
           setVar(nb,epot)
           % im MuPad-Fenster eintippen: d <Enter>
        end
    end
end
