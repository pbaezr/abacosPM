function [As,ys] = distribuirAs(b,h,caso,gamma,rho,rho_w)
% DISTRIBUIRAS distribuye la armadura en funci�n de la cuant�a y las
% propiedades geom�tricas de la secci�n.
%
% Variables de salida:
%   -As: vector fila con las �reas de acero para cada capa
%   -ys: vector fila con las ubicaciones de cada capa de armadura medida
%       desde la base hasta su centroide
%
% Variables de entrada:
%   -b: base de la secci�n
%   -h: altura de la secci�n
%   -caso: n�mero indicando el tipo de elemento (1, 2, 3, 4 � 5)
%   -gamma: distancia entre centroides de las capas de refuerzo extremas
%       cuantificada coomo fracci�n de h (valor entre 0 y 1)
%   -rho: cuant�a de acero respecto del �rea bruta (casos 1, 2, 3 y 5) o
%       del �rea del elemento de borde (caso 4)
%   -rho_w: cuant�a de la armadura de reparticion aplicable para el caso 4
%       (usualmente 0.0025 � 0.005)

if caso == 1 % columna con armadura en los bordes extremos
    As = rho*b*h*[0.5 0.5];
    ys = 0.5*h*[1-gamma 1+gamma];
elseif caso == 2 % columna con armadura lateral
    nc = 100; % n�mero de capas (se usa un n�mero elevado solo para lograr una buena distribuci�n)
    As = rho*b*h*repmat(1/nc,1,nc);
    ys = h*(0.5*(1-gamma)+gamma*(0:nc-1)/(nc-1));
elseif caso == 3 % columna con armadura perimetral
    nc = 100; % n�mero de capas (e igual al n�mero de barras por arista)
    As_int = rho*b*h*0.5/(nc-1); % armadura para una capa interna (que se replica para nc-2 capas)
    As_borde = rho*b*h*0.25*nc/(nc-1); % se considera aprox. un 25% de la armadura total para cada borde
    As = [As_borde repmat(As_int,1,nc-2) As_borde];
    ys = h*(0.5*(1-gamma)+gamma*(0:nc-1)/(nc-1));
elseif caso == 4 % muro con elementos de borde
    Lb = h*(1-gamma); % longitud del elemento de borde
    Lw = h-2*Lb; % longitud de la zona central con la armadura de repartici�n
    nci = 100; % n�mero de capas para la armadura de repartici�n
    s = Lw/nci; % espaciamiento de la armadura de repartici�n
    As_int = rho_w*b*s; % armadura para una capa interna (que se replica para nci capas)
    As_borde = rho*b*Lb;
    As = [As_borde repmat(As_int,1,nci) As_borde];    
    s1 = 0.5*s; % distancia entre el l�mite del elemento de borde y la primera (o �ltima) capa de armadura interna
    ys = [0.5*Lb Lb+s1+s*(0:nci-1) h-0.5*Lb];
else%if caso == 5 % muro con armadura distribuida uniformemente
    nc = 100; % n�mero de capas
    s = h/nc; % espaciamiento entre capas (medido entre centroides)
    As = rho*b*h*repmat(1/nc,1,nc);
    rec = 0.5*s;
    ys = rec+s*(0:nc-1);
end
