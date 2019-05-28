function [As,ys] = distribuirAs(b,h,caso,gamma,rho,rho_w)
% DISTRIBUIRAS distribuye la armadura en función de la cuantía y las
% propiedades geométricas de la sección.
%
% Variables de salida:
%   -As: vector fila con las áreas de acero para cada capa
%   -ys: vector fila con las ubicaciones de cada capa de armadura medida
%       desde la base hasta su centroide
%
% Variables de entrada:
%   -b: base de la sección
%   -h: altura de la sección
%   -caso: número indicando el tipo de elemento (1, 2, 3, 4 ó 5)
%   -gamma: distancia entre centroides de las capas de refuerzo extremas
%       cuantificada coomo fracción de h (valor entre 0 y 1)
%   -rho: cuantía de acero respecto del área bruta (casos 1, 2, 3 y 5) o
%       del área del elemento de borde (caso 4)
%   -rho_w: cuantía de la armadura de reparticion aplicable para el caso 4
%       (usualmente 0.0025 ó 0.005)

if caso == 1 % columna con armadura en los bordes extremos
    As = rho*b*h*[0.5 0.5];
    ys = 0.5*h*[1-gamma 1+gamma];
elseif caso == 2 % columna con armadura lateral
    nc = 100; % número de capas (se usa un número elevado solo para lograr una buena distribución)
    As = rho*b*h*repmat(1/nc,1,nc);
    ys = h*(0.5*(1-gamma)+gamma*(0:nc-1)/(nc-1));
elseif caso == 3 % columna con armadura perimetral
    nc = 100; % número de capas (e igual al número de barras por arista)
    As_int = rho*b*h*0.5/(nc-1); % armadura para una capa interna (que se replica para nc-2 capas)
    As_borde = rho*b*h*0.25*nc/(nc-1); % se considera aprox. un 25% de la armadura total para cada borde
    As = [As_borde repmat(As_int,1,nc-2) As_borde];
    ys = h*(0.5*(1-gamma)+gamma*(0:nc-1)/(nc-1));
elseif caso == 4 % muro con elementos de borde
    Lb = h*(1-gamma); % longitud del elemento de borde
    Lw = h-2*Lb; % longitud de la zona central con la armadura de repartición
    nci = 100; % número de capas para la armadura de repartición
    s = Lw/nci; % espaciamiento de la armadura de repartición
    As_int = rho_w*b*s; % armadura para una capa interna (que se replica para nci capas)
    As_borde = rho*b*Lb;
    As = [As_borde repmat(As_int,1,nci) As_borde];    
    s1 = 0.5*s; % distancia entre el límite del elemento de borde y la primera (o última) capa de armadura interna
    ys = [0.5*Lb Lb+s1+s*(0:nci-1) h-0.5*Lb];
else%if caso == 5 % muro con armadura distribuida uniformemente
    nc = 100; % número de capas
    s = h/nc; % espaciamiento entre capas (medido entre centroides)
    As = rho*b*h*repmat(1/nc,1,nc);
    rec = 0.5*s;
    ys = rec+s*(0:nc-1);
end
