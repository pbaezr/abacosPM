function [M,P,phi_red,et,c] = interaccionPM2(b,h,As,ys,parametros,n,tolerancia,incluirRamaTraccion,varargin)
%
% Variables de salida:
%   -M: vector fila con los momentos nominales asociados a cada carga axial (en kN*m)
%   -P: vector fila con las cargas axiales nominales utilizadas (en kN)
%   -phi_red: vector fila con los factores de reduccion de resistencia
%             calculados en función de la deformación axial del acero a tracción
%   -et: vector fila con la deformación unitaria de la capa de refuerzo más
%              alejada de la fibra más comprimido (negativa para tracción)
%   -c: vector fila con la ubicación de los ejes neutros (en mm)
%
% Variables de entrada:
%   -b: ancho de la sección rectangular (en mmm)
%   -h: altura de la sección rectangular (en mmm)
%   -As: vector fila con las áreas de refuerzo por cada capa (en mmm^2)
%   -ys: vector fila con la ubicación del centroide de cada capa de rufuerzo medida desde un borde (en mmm)
%   -parametros: estructura que contiene los parámetros que definen las leyes constitutivas de los materiales
%   -n: número de segmentos en que se dividirá la sección (debe ser un entero positivo)
%   -tolerancia: error permitido en la sumatoria de fuerzas (en kN)
%   -incluirRamaTraccion: valor lógico (true o false) que indica si se incluye en el diagrama la rama para P < 0
%   -varargin: vector con las cargas axiales (nominales) para las cuales se desea calcular el diagrama (en kN)

% parámetros de los materiales
Es1 = 200000; % módulo de elasticidad inicial del acero (en MPa)
fy = parametros.fy; % tensión de fluencia del acero (en MPa)
fcc = parametros.fc; % f'c del hormigón (en MPa)
Ecc = 4700*sqrt(fcc); % módulo de elasticidad del hormigón (en MPa)
ey = fy/Es1; % deformación de fluencia del acero (0.0021 para fy=420)
eu = 0.003; % deformación máxima para la cual se estima la resistencia de la sección, según ACI-318

% discretización elegida para el análisis seccional
dy = h/n; % largo de segmento
y = 0:dy:h; % distancia a cada segmento (a sus extremos) medida desde el borde inferior de la sección

% definición del rango de muestra para las cargas axiales
if ~isempty(varargin) % cargas definidas por el usuario
    P = 10^3*varargin{1}; % conversión de kN a N
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
et = zeros(1,np);% deformación en el acero a tracción (posible variable de salida)
phi_red = zeros(1,np); % factor de reducción de resiatencia (posible variable de salida)

tolerancia = 1000*tolerancia; % conversión de kN a N

% if parametros.modeloHormigon{1} ~= 0
Fc = zeros(1,n); % fuerzas en cada segmento de hormigón
yc = zeros(1,n); % brazos de palanca de cada fuerza Fc  

% inicializar los parámetros de la iteración
dp = 0; % error inicial considerado
eo = 0; % deformación unitaria inicial considerada (en h/2)
J = Ecc*b*h+2*Es1*sum((1-ys/h).*As);    

for i = 1:np
    error = tolerancia+1; % para entrar en un nuevo ciclo de iteraciones para la carga axial sig.
    numIteraciones = 0; % contador de iteraciones
    % iteraciones (se considerará esfuerzos/deformaciones/fuerzas positivas para la compresión)
    while error > tolerancia % test de convergencia
        % la no convergencia puede darse para cargas P muy altas que superan...
        % la compresión máxima que es capaz de resistir la sección
        if J == 0 || numIteraciones > 50
            if np > 1
                % redimensionar vectores
                M = M(1:i-1);
                P = P(1:i-1)/1000; % conversión de N a kN
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
        e = eo+2*(eu-eo)/h*(y-h/2); % deformaciones unitarias para los segementos de hormigón
        es = eo+2*(eu-eo)/h*(ys-h/2); % deformaciones en el acero

        % tensiones en el acero para la carga P(i) y la deformación en h/2, eo, recalculada
        [fs,Es] = curvaAcero(es,parametros);         
        Fs = fs.*As;

        % cálculo de tensiones, rigideces tangentes, fuerzas y brazos de palanca en cada segmento de hormigón
        [fc,Ec] = curvaHormigon(e,parametros);
        for k = 1:n
            Fc(k) = 0.5*(fc(k)+fc(k+1))*b*dy; % se usa un punto de integración --> área del trapecio en el perfil de tensiones
            if Fc(k) ~= 0, yc(k) = y(k)+dy/3*(fc(k)+2*fc(k+1))/(fc(k)+fc(k+1)); end % yc(k) = 0.5*(y(k)+y(k+1));
        end

        % cálculo de tensiones del hormigón en las zonas donde está distribuida la armadura
        % (para cuantías de acero pequeñas este cálculo es poco relevante)
        fc_As = curvaHormigon(es,parametros);
        Fc_ficticia = fc_As.*As; % estas fuerzas son descontadas para obtener la contribución real del hormigón

        % error
        dp = P(i)-(sum(Fs)+sum(Fc)-sum(Fc_ficticia));
        error = abs(dp);

        % antesala para nueva iteración
        J = 2*sum(Ec.*(1-y/h))*b*dy+2*sum(Es.*(1-ys/h).*As); % rigidez actualizada
        numIteraciones = numIteraciones+1;
    end

    % resultados para la carga P(i)
    M(i) = (sum(Fs.*ys)+sum(Fc.*yc)-sum(Fc_ficticia.*ys)-P(i)*h/2)*10^-6; % momento resultante en kN*m
    phi(i) = 2*(eu-eo)/h;
    c(i) = eu/phi(i);
    et(i) = eo+phi(i)*(min(ys)-h/2); % deformación unitaria en el acero a tracción para cada curvatura
    phi_red(i) = phiFactor(abs(et(i)),ey);
end
% else % modelo de Whitman (este método puede tener problemas de convergencia...
%      % para las cargas de tracción y para las compresiones bajas)
%     beta1 = max(0.65,min(0.85,0.85-0.008*(fcc-30)));
%     a = 0.1*h; % bloque de compresión uniforme (=beta1*c)
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
%                     P = P(1:i-1)/1000; % conversión de N a kN
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
%             a = min(a+da,h); % el bloque de compresión uniforme no puede tener una altura mayor que h
%             
%             es = eu+2*eu/(a/beta1)*(ys-h); % deformaciones en el acero
% 
%             % tensiones y fuerzas en el acero para la carga P(i)
%             fs = curvaAcero(es,parametros);
%             Fs = fs.*As;
% 
%             % cálculo de la contribución total del hormigón a compresión
%             Fc = 0.85*fcc*b*a;
% 
%             % cálculo de tensiones del hormigón en las zonas donde está distribuida la armadura
%             % (para cuantías de acero pequeñas este cálculo es poco relevante)
%             fc_As = 0.85*fcc*(es>eu*(1-beta1));
%             Fc_ficticia = fc_As.*As; % estas fuerzas son descontadas para obtener la contribución real del hormigón
%             
%             % error
%             dp = P(i)-(sum(Fs)+Fc-sum(Fc_ficticia));
%             error = abs(dp);
% 
%             % antesala para nueva iteración
%             da = dp/(0.85*fcc*b);
%             numIteraciones = numIteraciones+1;
%         end
% 
%         % Resultados para la carga P(i)
%         M(i) = (sum(Fs.*ys)+Fc*(h-a/2)-sum(Fc_ficticia.*ys)-P(i)*h/2)*10^-6; % momento resultante en kN*m        
%         c(i) = a/beta1;
%         phi(i) = eu/c(i);
%         et(i) = phi(i)*(h-min(ys)-c(i)); % deformación unitaria en el acero a tracción para cada curvatura
%         phi_red(i) = phiFactor(abs(et(i)),ey);
%     end    
% end

P = P/1000; % conversión de N a kN

% asegurar que el diagrama llegue a M=0
if np > 1 && ~cargasDefinidas
    [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion,parametros,Pmin/1000,Pmax/1000);
end

end

% función que asegura que el diagrama llegue a M=0 al interpolar (o extrapolar) para P
function [M,P,phi_red,et,c] = cerrarDiagrama(M,P,phi_red,et,c,incluirRamaTraccion,parametros,Pmin,Pmax)
% límite M=0, P=Pmax (compresión máxima)
if parametros.modeloHormigon{1} ~= 0, Pmax=interp1(M,P,0,'pchip'); end
ind = M>=0;
P = [P(ind) Pmax];
M = [M(ind) 0];
c = [c(ind) Inf];
phi_red = [phi_red(ind) 0.65];
et = [et(ind) 0.003];

% límite M=0, P=Pmin (tracción máxima)
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

% función que calcula el factor de reduccion de resistencia acorde a la deformación axial del acero a tracción
function phi = phiFactor(et,ey)
if et < ey % compresión controla
    phi = 0.65;
elseif et < 0.005 % zona de transición
    phi = 0.65+0.25*(et-ey)/(0.005-ey);
else%if et >= 0.005 % tracción controla
    phi = 0.9;
end
end

% función que calcula el esfuerzo y rigidez tangente en el acero
% dependiendo de la deformación unitaria y la ley constitutiva del material
function [fs,Es] = curvaAcero(es,parametros)

% parámetros de la curva tensión-deformación del acero de refuerzo
opcionAcero = parametros.modeloAcero; % tipo de curva
Es1 = 200000; % módulo de elasticidad inicial del acero (en MPa)
fy = parametros.fy; % tensión de fluencia del acero (en MPa)
ey = fy/Es1; % deformación unitaria de fluencia del acero

% inicializar variables
m = length(es);
fs = zeros(1,m); % tensiones en el acero
Es = zeros(1,m); % modulo de elasticidad tangente del acero (derivada de la funcion fs)

% calcular para cada capa de refuerzo
if opcionAcero{1} == 1 % modelo elastoplástico   
    ef = parametros.ef; % máxima deformación unitaria permitida
    Es2 = Es1*parametros.Es2; % módulo de elasticidad después tramo lineal-elástico del acero (en MPa)    
    
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
    esh = parametros.esh; % deformación para la cual inicia el endurecimiento del acero
    esu = parametros.esu; % deformación unitaria para la cual se genera la máxima tensión
    ef = parametros.ef; % máxima deformación unitaria permitida en el acero
    Esh = Es1*parametros.Esh;  % pendiente inicial post-fluencia del acero (en MPa) 
    fsu = fy*parametros.fsu; % máximo esfuerzo que puede alcanzar el acero (en MPa)
    p = Esh*(esu-esh)/(fsu-fy); % parámetro que define la curva post-fluencia

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

% función que calcula el esfuerzo y rigidez tangente en el hormigón
% dependiendo de la deformación unitaria y la ley constitutiva del material
function [fc,Ec] = curvaHormigon(e,parametros)

% parámetros de la curva esfuerzo-deformación del hormigón
opcionHormigon = parametros.modeloHormigon; % tipo de curva para la compresión
opcionHormigonTrac = parametros.modeloHormigonTrac; % tipo de curva para la tracción
fcc = parametros.fc; % f'c

e0 = parametros.e0; % deformación unitaria para la cual se genera el máximo esfuerzo
ef = parametros.ef; % máxima deformación unitaria permitida en el hormigón

% inicializar variables
n = length(e); % número de segmentos en que se divide la sección
fc = zeros(1,n); % tensiones en el hormigón
Ec = zeros(1,n); % modulo de elasticidad tangente del hormigón (derivada de la funcion fc)

% calcular para cada segmento de hormigón
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
            if e(i) <= e0 % tramo parabólico ascendente
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
            if opcionHormigon{1} == 3 % modelo de Thorenfeldt calibrado según Collins y Porasz
                r = 0.8+fcc/17;
                if e(i) > e0, k = 0.67+fcc/62; end
            else%opcionHormigon{1} == 4 % modelo de Thorenfeldt calibrado según Carreira y Kuang-Han
                r = 1.55+(fcc/32.4)^3;
            end
            fc(i) = fcc*r*(e(i)/e0)/(r-1+(e(i)/e0)^(r*k));
            Ec(i) = fcc*r/e0*(r-1+(1-r*k)*(e(i)/e0)^(r*k))/(r-1+(e(i)/e0)^(r*k))^2;
        end
    else
        if opcionHormigonTrac{1} == 1 % sin resistencia a tracción
            fc(i) = 0;
            Ec(i) = 0;
        else%if opcionHormigonTrac{1} == 2 % con resistencia lineal-elástica hasta la rotura 
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
