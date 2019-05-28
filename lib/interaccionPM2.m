function [M,P,phi_red,et,c] = interaccionPM2(b,h,As,ys,parametros,n,tolerancia,incluirRamaTraccion,varargin)
%
% Variables de salida:
%   -M: vector fila con los momentos nominales asociados a cada carga axial (en kN*m)
%   -P: vector fila con las cargas axiales nominales utilizadas (en kN)
%   -phi_red: vector fila con los factores de reduccion de resistencia
%             calculados en funci�n de la deformaci�n axial del acero a tracci�n
%   -et: vector fila con la deformaci�n unitaria de la capa de refuerzo m�s
%              alejada de la fibra m�s comprimido (negativa para tracci�n)
%   -c: vector fila con la ubicaci�n de los ejes neutros (en mm)
%
% Variables de entrada:
%   -b: ancho de la secci�n rectangular (en mmm)
%   -h: altura de la secci�n rectangular (en mmm)
%   -As: vector fila con las �reas de refuerzo por cada capa (en mmm^2)
%   -ys: vector fila con la ubicaci�n del centroide de cada capa de rufuerzo medida desde un borde (en mmm)
%   -parametros: estructura que contiene los par�metros que definen las leyes constitutivas de los materiales
%   -n: n�mero de segmentos en que se dividir� la secci�n (debe ser un entero positivo)
%   -tolerancia: error permitido en la sumatoria de fuerzas (en kN)
%   -incluirRamaTraccion: valor l�gico (true o false) que indica si se incluye en el diagrama la rama para P < 0
%   -varargin: vector con las cargas axiales (nominales) para las cuales se desea calcular el diagrama (en kN)

% par�metros de los materiales
Es1 = 200000; % m�dulo de elasticidad inicial del acero (en MPa)
fy = parametros.fy; % tensi�n de fluencia del acero (en MPa)
fcc = parametros.fc; % f'c del hormig�n (en MPa)
Ecc = 4700*sqrt(fcc); % m�dulo de elasticidad del hormig�n (en MPa)
ey = fy/Es1; % deformaci�n de fluencia del acero (0.0021 para fy=420)
eu = 0.003; % deformaci�n m�xima para la cual se estima la resistencia de la secci�n, seg�n ACI-318

% discretizaci�n elegida para el an�lisis seccional
dy = h/n; % largo de segmento
y = 0:dy:h; % distancia a cada segmento (a sus extremos) medida desde el borde inferior de la secci�n

% definici�n del rango de muestra para las cargas axiales
if ~isempty(varargin) % cargas definidas por el usuario
    P = 10^3*varargin{1}; % conversi�n de kN a N
    np = length(P);
    cargasDefinidas = true;
else
    np = 500;
    Pmin = -fy*sum(As);    
    Pmax = 0.85*fcc*(b*h-sum(As))+fy*sum(As);
    factorPmin = 0;
    if incluirRamaTraccion
        factorPmin = 0.8; % puede no haber convergencia para Pmin = -fy*sum(As)
        np = 1000;
    end    
    P = sort([linspace(factorPmin*Pmin,1.15*Pmax,np-1) 0.8*Pmax]);
    cargasDefinidas = false;
end

% inicializar variables
M = zeros(1,np); % momentos resultantes para cada curvatura
phi = zeros(1,np); % curvatura resultante para cada carga (posible variable de salida)
c = zeros(1,np); %
et = zeros(1,np);% deformaci�n en el acero a tracci�n (posible variable de salida)
phi_red = zeros(1,np); % factor de reducci�n de resiatencia (posible variable de salida)

tolerancia = 1000*tolerancia; % conversi�n de kN a N

% if parametros.modeloHormigon{1} ~= 0
Fc = zeros(1,n); % fuerzas en cada segmento de hormig�n
yc = zeros(1,n); % brazos de palanca de cada fuerza Fc  

% inicializar los par�metros de la iteraci�n
dp = 0; % error inicial considerado
eo = 0; % deformaci�n unitaria inicial considerada (en h/2)
J = Ecc*b*h+2*Es1*sum((1-ys/h).*As);    

for i = 1:np
    error = tolerancia+1; % para entrar en un nuevo ciclo de iteraciones para la carga axial sig.
    numIteraciones = 0; % contador de iteraciones
    % iteraciones (se considerar� esfuerzos/deformaciones/fuerzas positivas para la compresi�n)
    while error > tolerancia % test de convergencia
        % la no convergencia puede darse para cargas P muy altas que superan...
        % la compresi�n m�xima que es capaz de resistir la secci�n
        if J == 0 || numIteraciones > 50
            if np > 1
                % redimensionar vectores
                M = M(1:i-1);
                P = P(1:i-1)/1000; % conversi�n de N a kN
                c = c(1:i-1);
                phi_red = phi_red(1:i-1);
                et = et(1:i-1);

                % asegurar que el diagrama llegue a M=0 interpolando (o extrapolando) para P
                if ~cargasDefinidas
                    [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion,parametros,Pmin/1000,Pmax/1000);
                end
            end

            % entregar los resultados
            return
        end

        deo = J^-1*dp;
        eo = eo+deo;
        e = eo+2*(eu-eo)/h*(y-h/2); % deformaciones unitarias para los segementos de hormig�n
        es = eo+2*(eu-eo)/h*(ys-h/2); % deformaciones en el acero

        % tensiones en el acero para la carga P(i) y la deformaci�n en h/2, eo, recalculada
        [fs,Es] = curvaAcero(es,parametros);         
        Fs = fs.*As;

        % c�lculo de tensiones, rigideces tangentes, fuerzas y brazos de palanca en cada segmento de hormig�n
        [fc,Ec] = curvaHormigon(e,parametros);
        for k = 1:n
            Fc(k) = 0.5*(fc(k)+fc(k+1))*b*dy; % se usa un punto de integraci�n --> �rea del trapecio en el perfil de tensiones
            if Fc(k) ~= 0, yc(k) = y(k)+dy/3*(fc(k)+2*fc(k+1))/(fc(k)+fc(k+1)); end % yc(k) = 0.5*(y(k)+y(k+1));
        end

        % c�lculo de tensiones del hormig�n en las zonas donde est� distribuida la armadura
        % (para cuant�as de acero peque�as este c�lculo es poco relevante)
        fc_As = curvaHormigon(es,parametros);
        Fc_ficticia = fc_As.*As; % estas fuerzas son descontadas para obtener la contribuci�n real del hormig�n

        % error
        dp = P(i)-(sum(Fs)+sum(Fc)-sum(Fc_ficticia));
        error = abs(dp);

        % antesala para nueva iteraci�n
        J = 2*sum(Ec.*(1-y/h))*b*dy+2*sum(Es.*(1-ys/h).*As); % rigidez actualizada
        numIteraciones = numIteraciones+1;
    end

    % resultados para la carga P(i)
    M(i) = (sum(Fs.*ys)+sum(Fc.*yc)-sum(Fc_ficticia.*ys)-P(i)*h/2)*10^-6; % momento resultante en kN*m
    phi(i) = 2*(eu-eo)/h;
    c(i) = eu/phi(i);
    et(i) = eo+phi(i)*(min(ys)-h/2); % deformaci�n unitaria en el acero a tracci�n para cada curvatura
    phi_red(i) = phiFactor(abs(et(i)),ey);
end
% else % modelo de Whitman (este m�todo puede tener problemas de convergencia...
%      % para las cargas de tracci�n y para las compresiones bajas)
%     beta1 = max(0.65,min(0.85,0.85-0.008*(fcc-30)));
%     a = 0.1*h; % bloque de compresi�n uniforme (=beta1*c)
%     da = 0;
%     for i = 1:np
%         error = tolerancia+1;
%         numIteraciones = 0;
%         
%         while error>tolerancia % test de convergencia
%             if numIteraciones > 50
%                 if np > 1
%                     % redimensionar vectores
%                     M = M(1:i-1);
%                     P = P(1:i-1)/1000; % conversi�n de N a kN
%                     c = c(1:i-1);
%                     phi_red = phi_red(1:i-1);
%                     et = et(1:i-1);
% 
%                     % asegurar que el diagrama llegue a M=0
%                     if ~cargasDefinidas
%                         [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion,parametros,Pmin/1000,Pmax/1000);
%                     end
%                 end
% 
%                 % entregar los resultados
%                 return
%             end            
%             
%             a = min(a+da,h); % el bloque de compresi�n uniforme no puede tener una altura mayor que h
%             
%             es = eu+2*eu/(a/beta1)*(ys-h); % deformaciones en el acero
% 
%             % tensiones y fuerzas en el acero para la carga P(i)
%             fs = curvaAcero(es,parametros);
%             Fs = fs.*As;
% 
%             % c�lculo de la contribuci�n total del hormig�n a compresi�n
%             Fc = 0.85*fcc*b*a;
% 
%             % c�lculo de tensiones del hormig�n en las zonas donde est� distribuida la armadura
%             % (para cuant�as de acero peque�as este c�lculo es poco relevante)
%             fc_As = 0.85*fcc*(es>eu*(1-beta1));
%             Fc_ficticia = fc_As.*As; % estas fuerzas son descontadas para obtener la contribuci�n real del hormig�n
%             
%             % error
%             dp = P(i)-(sum(Fs)+Fc-sum(Fc_ficticia));
%             error = abs(dp);
% 
%             % antesala para nueva iteraci�n
%             da = dp/(0.85*fcc*b);
%             numIteraciones = numIteraciones+1;
%         end
% 
%         % Resultados para la carga P(i)
%         M(i) = (sum(Fs.*ys)+Fc*(h-a/2)-sum(Fc_ficticia.*ys)-P(i)*h/2)*10^-6; % momento resultante en kN*m        
%         c(i) = a/beta1;
%         phi(i) = eu/c(i);
%         et(i) = phi(i)*(h-min(ys)-c(i)); % deformaci�n unitaria en el acero a tracci�n para cada curvatura
%         phi_red(i) = phiFactor(abs(et(i)),ey);
%     end    
% end

P = P/1000; % conversi�n de N a kN

% asegurar que el diagrama llegue a M=0
if np > 1 && ~cargasDefinidas
    [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion,parametros,Pmin/1000,Pmax/1000);
end

end

% funci�n que asegura que el diagrama llegue a M=0 al interpolar (o extrapolar) para P
function [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion,parametros,Pmin,Pmax)
% l�mite M=0, P=Pmax (compresi�n m�xima)
if parametros.modeloHormigon{1} ~= 0, Pmax=interp1(M,P,0,'pchip'); end
ind = M>=0;
P = [P(ind) Pmax];
M = [M(ind) 0];
c = [c(ind) Inf];
phi_red = [phi_red(ind) 0.65];
et = [et(ind) 0.003];

% l�mite M=0, P=Pmin (tracci�n m�xima)
if length(P) > 1
    if incluirRamaTraccion
        if parametros.modeloHormigon{1} ~= 0
            try
                Pmin = interp1(M(1:2),P(1:2),0,'pchip');
            catch
            end
        end        
        et = [interp1(P,et,Pmin,'pchip') et];        
        P = [Pmin P];
        M = [0 M];
        c = [-Inf c];
        phi_red = [0.9 phi_red];
    end
end
end

% funci�n que calcula el factor de reduccion de resistencia acorde a la deformaci�n axial del acero a tracci�n
function phi = phiFactor(et,ey)
if et < ey % compresi�n controla
    phi = 0.65;
elseif et < 0.005 % zona de transici�n
    phi = 0.65+0.25*(et-ey)/(0.005-ey);
else%if et >= 0.005 % tracci�n controla
    phi = 0.9;
end
end

% funci�n que calcula el esfuerzo y rigidez tangente en el acero
% dependiendo de la deformaci�n unitaria y la ley constitutiva del material
function [fs,Es] = curvaAcero(es,parametros)

% par�metros de la curva tensi�n-deformaci�n del acero de refuerzo
opcionAcero = parametros.modeloAcero; % tipo de curva
Es1 = 200000; % m�dulo de elasticidad inicial del acero (en MPa)
fy = parametros.fy; % tensi�n de fluencia del acero (en MPa)
ey = fy/Es1; % deformaci�n unitaria de fluencia del acero

% inicializar variables
m = length(es);
fs = zeros(1,m); % tensiones en el acero
Es = zeros(1,m); % modulo de elasticidad tangente del acero (derivada de la funcion fs)

% calcular para cada capa de refuerzo
if opcionAcero{1} == 1 % modelo elastopl�stico   
    ef = parametros.ef; % m�xima deformaci�n unitaria permitida
    Es2 = Es1*parametros.Es2; % m�dulo de elasticidad despu�s tramo lineal-el�stico del acero (en MPa)    
    
    for i = 1:m
        if abs(es(i)) <= ey
            fs(i) = Es1*es(i);
            Es(i) = Es1;
        elseif abs(es(i)) <= ef
            fs(i) = sign(es(i))*(fy+Es2*(abs(es(i))-ey));
            Es(i) = Es2;
        else
            fs(i) = 0;
            Es(i) = 0;
        end
    end
else%if opcionAcero{1} == 2 % modelo de Mander
    esh = parametros.esh; % deformaci�n para la cual inicia el endurecimiento del acero
    esu = parametros.esu; % deformaci�n unitaria para la cual se genera la m�xima tensi�n
    ef = parametros.ef; % m�xima deformaci�n unitaria permitida en el acero
    Esh = Es1*parametros.Esh;  % pendiente inicial post-fluencia del acero (en MPa) 
    fsu = fy*parametros.fsu; % m�ximo esfuerzo que puede alcanzar el acero (en MPa)
    p = Esh*(esu-esh)/(fsu-fy); % par�metro que define la curva post-fluencia

    for i = 1:m
        if abs(es(i)) <= ey
            fs(i) = Es1*es(i);
            Es(i) = Es1;
        elseif abs(es(i)) <= esh
            fs(i) = sign(es(i))*fy;
            Es(i) = 0;
        elseif abs(es(i)) <= ef
            fs(i) = sign(es(i))*(fsu+(fy-fsu)*abs((esu-abs(es(i)))/(esu-esh))^p);
            Es(i) = p*(fsu-fy)*abs(esu-abs(es(i)))^(p-1)/(esu-esh)^p;
        else
            fs(i) = 0;
            Es(i) = 0;
        end
    end
end
end

% funci�n que calcula el esfuerzo y rigidez tangente en el hormig�n
% dependiendo de la deformaci�n unitaria y la ley constitutiva del material
function [fc,Ec] = curvaHormigon(e,parametros)

% par�metros de la curva esfuerzo-deformaci�n del hormig�n
opcionHormigon = parametros.modeloHormigon; % tipo de curva para la compresi�n
opcionHormigonTrac = parametros.modeloHormigonTrac; % tipo de curva para la tracci�n
fcc = parametros.fc; % f'c

e0 = parametros.e0; % deformaci�n unitaria para la cual se genera el m�ximo esfuerzo
ef = parametros.ef; % m�xima deformaci�n unitaria permitida en el hormig�n

% inicializar variables
n = length(e); % n�mero de segmentos en que se divide la secci�n
fc = zeros(1,n); % tensiones en el hormig�n
Ec = zeros(1,n); % modulo de elasticidad tangente del hormig�n (derivada de la funcion fc)

% calcular para cada segmento de hormig�n
for i = 1:n
    if e(i) >= 0
        if opcionHormigon{1} == 1 % modelo de Saenz
            if e(i) <= 2*e0
                fc(i) = fcc*(2*e(i)/e0-(e(i)/e0)^2);
                Ec(i) = 2*(fcc/e0)*(1-e(i)/e0);
            else
                fc(i) = 0;
                Ec(i) = 0;
            end
        elseif opcionHormigon{1} == 2 % modelo de Hognestad
            if e(i) <= e0 % tramo parab�lico ascendente
                fc(i) = fcc*(2*e(i)/e0-(e(i)/e0)^2);
                Ec(i) = 2*(fcc/e0)*(1-e(i)/e0);
            elseif e(i) <= ef % tramo lineal descendente
                fc(i) = fcc*(1-0.15*(e(i)-e0)/(ef-e0));
                Ec(i) = -0.15*fcc/(ef-e0);
            else
                fc(i) = 0;
                Ec(i) = 0;
            end
        else
            k = 1;
            if opcionHormigon{1} == 3 % modelo de Thorenfeldt calibrado seg�n Collins y Porasz
                r = 0.8+fcc/17;
                if e(i) > e0, k = 0.67+fcc/62; end
            else%opcionHormigon{1} == 4 % modelo de Thorenfeldt calibrado seg�n Carreira y Kuang-Han
                r = 1.55+(fcc/32.4)^3;
            end
            fc(i) = fcc*r*(e(i)/e0)/(r-1+(e(i)/e0)^(r*k));
            Ec(i) = fcc*r/e0*(r-1+(1-r*k)*(e(i)/e0)^(r*k))/(r-1+(e(i)/e0)^(r*k))^2;
        end
    else
        if opcionHormigonTrac{1} == 1 % sin resistencia a tracci�n
            fc(i) = 0;
            Ec(i) = 0;
        else%if opcionHormigonTrac{1} == 2 % con resistencia lineal-el�stica hasta la rotura 
            if e(i) <= -0.62/4700
                fc(i) = 0;
                Ec(i) = 0;
            else
                Ecc = 0.62*sqrt(4700);
                fc(i) = Ecc*e(i);
                Ec(i) = Ecc;
            end
        end
    end
end
end
