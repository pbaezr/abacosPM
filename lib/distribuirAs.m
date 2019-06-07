function [As,ys] = distribuirAs(b,h,caso,gamma,rho,rho_w)
% DISTRIBUIRAS distribuye la armadura en funcion de la cuantia y las
% propiedades geometricas de la seccion.
%
% Variables de salida:
%   -As: vector fila con las areas de acero para cada capa
%   -ys: vector fila con las ubicaciones de cada capa de armadura medida
%       desde la base hasta su centroide
%
% Variables de entrada:
%   -b: base de la seccion rectangular (o diametro si es circular)
%   -h: altura de la seccion
%   -caso: numero indicando el tipo de elemento (1, 2, 3, 3.5, 4 o 5)
%   -gamma: distancia entre centroides de las capas de refuerzo extremas
%       cuantificada coomo fraccion de h (valor entre 0 y 1)
%   -rho: cuantia de acero respecto del area bruta (casos 1, 2, 3, 3.5 y 5) o
%       del area del elemento de borde (caso 4)
%   -rho_w: cuantia de la armadura de reparticion aplicable para el caso 4
%       (usualmente 0.0025 o 0.005)

if caso == 1 % columna con armadura en los bordes extremos
    As = rho*b*h*[0.5 0.5];
    ys = 0.5*h*[1-gamma 1+gamma];
elseif caso == 2 % columna con armadura lateral
    nc = 100; % numero de capas (se usa un numero elevado solo para lograr una buena distribucion)
    As = rho*b*h*repmat(1/nc,1,nc);
    ys = h*(0.5*(1-gamma)+gamma*(0:nc-1)/(nc-1));
elseif caso == 3 % columna con armadura perimetral
    nc = 100; % numero de capas (e igual al numero de barras por arista)
    As_int = rho*b*h*0.5/(nc-1); % armadura para una capa interna (que se replica para nc-2 capas)
    As_borde = rho*b*h*0.25*nc/(nc-1); % se considera aprox. un 25% de la armadura total para cada borde
    As = [As_borde repmat(As_int,1,nc-2) As_borde];
    ys = h*(0.5*(1-gamma)+gamma*(0:nc-1)/(nc-1));
elseif caso == 3.5 % columna circular con armadura perimetral
    nb = 100; % numero de barras
    Ab = rho*(0.25*pi*h^2)/nb; % area de un solo fierro
    As = Ab*[1 2*ones(1,nb/2-1) 1];
    ys = 0.5*h*(1-gamma*cos(0:2*pi/nb:pi));
elseif caso == 4 % muro con elementos de borde
    Lb = h*(1-gamma); % longitud del elemento de borde
    Lw = h-2*Lb; % longitud de la zona central con la armadura de reparticion
    nci = 100; % numero de capas para la armadura de reparticion
    s = Lw/nci; % espaciamiento de la armadura de reparticion
    As_int = rho_w*b*s; % armadura para una capa interna (que se replica para nci capas)
    As_borde = rho*b*Lb;
    As = [As_borde repmat(As_int,1,nci) As_borde];    
    s1 = 0.5*s; % distancia entre el limite del elemento de borde y la primera (o ultima) capa de armadura interna
    ys = [0.5*Lb Lb+s1+s*(0:nci-1) h-0.5*Lb];
else%if caso == 5 % muro con armadura distribuida uniformemente
    nc = 100; % numero de capas
    s = h/nc; % espaciamiento entre capas (medido entre centroides)
    As = rho*b*h*repmat(1/nc,1,nc);
    rec = 0.5*s;
    ys = rec+s*(0:nc-1);
end
