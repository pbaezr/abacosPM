function [M,P,phi_red,et,c] = interaccionPM2(b,h,As,ys,parametros,n,tolerancia,incluirRamaTraccion,varargin)
%
% Variables de salida:
%   -M: vector fila con los momentos nominales asociados a cada carga axial (en kN*m)
%   -P: vector fila con las cargas axiales nominales utilizadas (en kN)
%   -phi_red: vector fila con los factores de reduccion de resistencia
%             calculados en funcion de la deformacion axial del acero a traccion
%   -et: vector fila con la deformacion unitaria de la capa de refuerzo mas
%              alejada de la fibra mas comprimido (negativa para traccion)
%   -c: vector fila con la ubicacion de los ejes neutros (en mm)
%
% Variables de entrada:
%   -b: ancho de la seccion rectangular (en mmm)
%   -h: altura de la seccion rectangular (en mmm)
%   -As: vector fila con las areas de refuerzo por cada capa (en mmm^2)
%   -ys: vector fila con la ubicacion del centroide de cada capa de rufuerzo medida desde un borde (en mmm)
%   -parametros: estructura que contiene los parametros que definen las leyes constitutivas de los materiales
%   -n: numero de segmentos en que se dividira la seccion (debe ser un entero positivo)
%   -tolerancia: error permitido en la sumatoria de fuerzas (en kN)
%   -incluirRamaTraccion: valor logico (true o false) que indica si se incluye en el diagrama la rama para P < 0
%   -varargin: vector con las cargas axiales (nominales) para las cuales se desea calcular el diagrama (en kN)

% parametros de los materiales
Es1 = 200000; % modulo de elasticidad inicial del acero (en MPa)
fy = parametros.fy; % tension de fluencia del acero (en MPa)
fcc = parametros.fc; % f'c del hormigon (en MPa)
Ecc = 4700*sqrt(fcc); % modulo de elasticidad del hormigon (en MPa)
ey = fy/Es1; % deformacion de fluencia del acero (0.0021 para fy=420)
eu = 0.003; % deformacion maxima para la cual se estima la resistencia de la seccion, segun ACI-318

% discretizacion elegida para el analisis seccional
dy = h/n; % altura de segmento
y = dy/2:dy:h-dy/2; % distancia a cada segmento (en sus centros) medida desde el borde inferior de la seccion

% definicion del rango de muestra para las cargas axiales
if ~isempty(varargin) % cargas definidas por el usuario
    P = 10^3*varargin{1}; % conversion de kN a N
    np = length(P);
    cargasDefinidas = true;
else
    np = 500;
    Pmin = -fy*sum(As);    
    Pmax = 0.85*fcc*(b*h-sum(As))+fy*sum(As);
    factorPmin = 0;
    if incluirRamaTraccion
        % exageradamente se supone que la tracción pura ocurre cuando la
        % tensión es máxima
        if parametros.modeloAcero{1} == 1
            factorPmin = 1+parametros.Es2/ey*(parametros.ef-ey);
        else
            factorPmin = parametros.fsu;
        end
        np = 1000;
    end    
    P = sort([linspace(factorPmin*Pmin,1.15*Pmax,np-1) 0.8*Pmax]);
    cargasDefinidas = false;
end

% inicializar variables
M = zeros(1,np); % momentos resultantes para cada curvatura
phi = zeros(1,np); % curvatura resultante para cada carga
c = zeros(1,np); %
et = zeros(1,np);% deformacion en el acero a traccion
phi_red = zeros(1,np); % factor de reduccion de resiatencia
indSaltoP = zeros(1,np); % indices de cargas excluidas

tolerancia = 1000*tolerancia; % conversion de kN a N

% if parametros.modeloHormigon{1} ~= 0  

% inicializar los parametros de la iteracion
dp = 0; % error inicial considerado
eo = 0; % deformacion unitaria inicial considerada (en h/2)
J0 = Ecc*b*h+2*Es1*sum((1-ys/h).*As);    

% iteraciones
i  = 1;
J = J0;
suma = 0;
while i <= np
    error = tolerancia+1; % para entrar en un nuevo ciclo de iteraciones para la carga axial sig.
    numIteraciones = 0; % contador de iteraciones
    % iteraciones (se considerara esfuerzos/deformaciones/fuerzas positivas para la compresion)
    while error > tolerancia % test de convergencia
        numIteraciones = numIteraciones+1;

        deo = J^-1*dp;
        eo = eo+deo;
        e = eo+2*(eu-eo)/h*(y-h/2); % deformaciones unitarias para los segementos de hormigon
        es = eo+2*(eu-eo)/h*(ys-h/2); % deformaciones en el acero

        % tensiones en el acero para la carga P(i) y la deformacion en h/2, eo, recalculada
        [fs,Es] = curvaAcero(es,parametros);         
        Fs = fs.*As;

        % calculo de tensiones, rigideces tangentes, fuerzas y brazos de palanca en cada segmento de hormigon
        [fc,Ec] = curvaHormigon(e,parametros);
        Fc = fc*b*dy;

        % calculo de tensiones del hormigon en las zonas donde esta distribuida la armadura
        % (para cuantias de acero pequenas este calculo es poco relevante)
        fc_As = curvaHormigon(es,parametros);
        Fc_ficticia = fc_As.*As; % estas fuerzas son descontadas para obtener la contribucion real del hormigon

        % error
        dp = P(i)-(sum(Fs)+sum(Fc)-sum(Fc_ficticia));
        error = abs(dp);

        % actualizar rigidez J para nueva iteracion
        if numIteraciones < 11            
            Janterior = J; % almacenar J de la iteracion previa para comparar
            J = 2*sum(Ec.*(1-y/h))*b*dy+2*sum(Es.*(1-ys/h).*As);
            if J == Janterior % obviar carga y reiniciar los parametros de iteracion
                J = J0;
                indSaltoP(k) = i;
                k = k+1;
                eo = 0;
                dp = 0;
                break
            end            
        else % obviar carga si no hay convergencia y continuar con la siguiente usando el J inicial
            J = J0;
            indSaltoP(k) = i;
            k = k+1;
            break
        end
    end

    % resultados para la carga P(i)
    if k == 1 || i ~= indSaltoP(k-1)
        M(i) = (sum(Fs.*ys)+sum(Fc.*y)-sum(Fc_ficticia.*ys)-P(i)*h/2)*10^-6; % momento resultante en kN*m
        phi(i) = 2*(eu-eo)/h;
        c(i) = eu/phi(i);
        et(i) = eo+phi(i)*(min(ys)-h/2); % deformacion unitaria en el acero a traccion para cada curvatura
        phi_red(i) = phiFactor(abs(et(i)),ey);
    end
    
    i = i+1;
end
% else % modelo de Whitman (este metodo puede tener problemas de convergencia...
%      % para las cargas de traccion y para las compresiones bajas)
%     beta1 = max(0.65,min(0.85,0.85-0.008*(fcc-30)));
%     a = 0.1*h; % bloque de compresion uniforme (=beta1*c)
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
%                     P = P(1:i-1)/1000; % conversion de N a kN
%                     c = c(1:i-1);
%                     phi_red = phi_red(1:i-1);
%                     et = et(1:i-1);
% 
%                     % asegurar que el diagrama llegue a M=0
%                     if ~cargasDefinidas
%                         [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion);
%                     end
%                 end
% 
%                 % entregar los resultados
%                 return
%             end            
%             
%             a = min(a+da,h); % el bloque de compresion uniforme no puede tener una altura mayor que h
%             
%             es = eu+2*eu/(a/beta1)*(ys-h); % deformaciones en el acero
% 
%             % tensiones y fuerzas en el acero para la carga P(i)
%             fs = curvaAcero(es,parametros);
%             Fs = fs.*As;
% 
%             % calculo de la contribucion total del hormigon a compresion
%             Fc = 0.85*fcc*b*a;
% 
%             % calculo de tensiones del hormigon en las zonas donde esta distribuida la armadura
%             % (para cuantias de acero pequenas este calculo es poco relevante)
%             fc_As = 0.85*fcc*(es>eu*(1-beta1));
%             Fc_ficticia = fc_As.*As; % estas fuerzas son descontadas para obtener la contribucion real del hormigon
%             
%             % error
%             dp = P(i)-(sum(Fs)+Fc-sum(Fc_ficticia));
%             error = abs(dp);
% 
%             % antesala para nueva iteracion
%             da = dp/(0.85*fcc*b);
%             numIteraciones = numIteraciones+1;
%         end
% 
%         % Resultados para la carga P(i)
%         M(i) = (sum(Fs.*ys)+Fc*(h-a/2)-sum(Fc_ficticia.*ys)-P(i)*h/2)*10^-6; % momento resultante en kN*m        
%         c(i) = a/beta1;
%         phi(i) = eu/c(i);
%         et(i) = phi(i)*(h-min(ys)-c(i)); % deformacion unitaria en el acero a traccion para cada curvatura
%         phi_red(i) = phiFactor(abs(et(i)),ey);
%     end    
% end

% redimensionar vectores excluyendo los datos no validos
ind = ~ismember(1:np,indSaltoP) & ~isnan(M) & M~=0 & M>0;
M = M(ind);
P = P(ind)/1000; % conversión de N a kN
c = c(ind);
phi_red = phi_red(ind);
et = et(ind);

% asegurar que el diagrama llegue a M=0
if np > 1 && ~cargasDefinidas
    [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion);
end

end

% funcion que asegura que el diagrama llegue a M=0 al interpolar (o extrapolar) para P
function [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion,parametros,Pmin,Pmax)
% limite compresion maxima
Pmax = interp1(M(P>0),P(P>0),0,'linear','extrap');

% limite traccion maxima (estas interpolaciones son solo para fines graficos
% dado que en teoria ya fue incluido el ultimo par P-M en el que la seccion se encuentra en equilibrio
if incluirRamaTraccion
    Pmin = interp1(M(P<0),P(P<0),0,'linear','extrap');
    ett = interp1(P(P<0),et(P<0),Pmin,'linear','extrap');
    phit = 0.9;
    Mmin = 0;
else
    Pmin = [];
    ett = [];
    phit = [];
    Mmin = [];
end

% redimensionar vectores incluyendo los límites
M = [Mmin M 0];
c = [-Inf c Inf];
et = [ett et interp1(P(P>0),et(P>0),Pmax,'linear','extrap')];
phi_red = [phit phi_red 0.65];
P = [Pmin P Pmax];
end

% funcion que calcula el factor de reduccion de resistencia acorde a la deformacion axial del acero a traccion
function phi = phiFactor(et,ey)
if et < ey % compresion controla
    phi = 0.65;
elseif et < 0.005 % zona de transicion
    phi = 0.65+0.25*(et-ey)/(0.005-ey);
else%if et >= 0.005 % traccion controla
    phi = 0.9;
end
end

% funcion que calcula el esfuerzo y rigidez tangente en el acero
% dependiendo de la deformacion unitaria y la ley constitutiva del material
function [fs,Es] = curvaAcero(es,parametros)

% parametros de la curva tension-deformacion del acero de refuerzo
opcionAcero = parametros.modeloAcero; % tipo de curva
Es1 = 200000; % modulo de elasticidad inicial del acero (en MPa)
fy = parametros.fy; % tension de fluencia del acero (en MPa)
ey = fy/Es1; % deformacion unitaria de fluencia del acero

% inicializar variables
m = length(es);
fs = zeros(1,m); % tensiones en el acero
Es = zeros(1,m); % modulo de elasticidad tangente del acero (derivada de la funcion fs)

% calcular para cada capa de refuerzo
if opcionAcero{1} == 1 % modelo elastoplastico   
    ef = parametros.ef; % maxima deformacion unitaria permitida
    Es2 = Es1*parametros.Es2; % modulo de elasticidad despues tramo lineal-elastico del acero (en MPa)    
    
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
    esh = parametros.esh; % deformacion para la cual inicia el endurecimiento del acero
    esu = parametros.esu; % deformacion unitaria para la cual se genera la maxima tension
    ef = parametros.ef; % maxima deformacion unitaria permitida en el acero
    Esh = Es1*parametros.Esh;  % pendiente inicial post-fluencia del acero (en MPa) 
    fsu = fy*parametros.fsu; % maximo esfuerzo que puede alcanzar el acero (en MPa)
    p = Esh*(esu-esh)/(fsu-fy); % parametro que define la curva post-fluencia

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

% funcion que calcula el esfuerzo y rigidez tangente en el hormigon
% dependiendo de la deformacion unitaria y la ley constitutiva del material
function [fc,Ec] = curvaHormigon(e,parametros)

% parametros de la curva esfuerzo-deformacion del hormigon
opcionHormigon = parametros.modeloHormigon; % tipo de curva para la compresion
opcionHormigonTrac = parametros.modeloHormigonTrac; % tipo de curva para la traccion
fcc = parametros.fc; % f'c

e0 = parametros.e0; % deformacion unitaria para la cual se genera el maximo esfuerzo
ef = parametros.ef; % maxima deformacion unitaria permitida en el hormigon

% inicializar variables
n = length(e); % numero de segmentos en que se divide la seccion
fc = zeros(1,n); % tensiones en el hormigon
Ec = zeros(1,n); % modulo de elasticidad tangente del hormigon (derivada de la funcion fc)

% calcular para cada segmento de hormigon
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
            if e(i) <= e0 % tramo parabolico ascendente
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
            if opcionHormigon{1} == 3 % modelo de Thorenfeldt calibrado segun Collins y Porasz
                r = 0.8+fcc/17;
                if e(i) > e0, k = 0.67+fcc/62; end
            else%if opcionHormigon{1} == 4 % modelo de Thorenfeldt calibrado segun Carreira y Kuang-Han
                r = 1.55+(fcc/32.4)^3;
            end
            fc(i) = fcc*r*(e(i)/e0)/(r-1+(e(i)/e0)^(r*k));
            Ec(i) = fcc*r/e0*(r-1+(1-r*k)*(e(i)/e0)^(r*k))/(r-1+(e(i)/e0)^(r*k))^2;
        end
    else
        if opcionHormigonTrac{1} == 1 % sin resistencia a traccion
            fc(i) = 0;
            Ec(i) = 0;
        else%if opcionHormigonTrac{1} == 2 % con resistencia lineal-elastica hasta la rotura 
            if e(i) <= -0.62/4700
                fc(i) = 0;
                Ec(i) = 0;
            else
                Ecc = 4700*sqrt(fcc);
                fc(i) = Ecc*e(i);
                Ec(i) = Ecc;
            end
        end
    end
end
end
