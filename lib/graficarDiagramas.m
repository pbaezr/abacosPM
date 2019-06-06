function graficarDiagramas(handles,varargin)
% Esta funcion grafica los diagramas para las cuantias de refuerzo de 1 a 8%

% Licenciado bajos los terminos del MIT.
% Copyright (c) 2019 Pablo Baez R.

handles.figure1.Pointer = 'watch';

% borrar resultados antiguos
axes(handles.axes1), cla
handles.output_rho.setText('');
handles.output_As.setText('');

% definir valores arbitrarios para la seccion
% (no influyen en el calculo pues M y P son normalizados por el area de la seccion)
b = 100; % ancho de la seccion
h = 100; % altura de la seccion

% obtener valores de los inputs
gamma = handles.input_gamma.UserData; % ubicacion de la armadura As2 (en mm)
rho_w = handles.input_rho_w.UserData;
n = handles.input_n.UserData;
tol = handles.input_error.UserData;
parametros = getappdata(handles.figure1,'parametrosMateriales');
fc = parametros.fc;
fy = parametros.fy;

if ~isempty(n) && ~isempty(tol)
    opcion = handles.uibuttongroup1.UserData;    
    
    grid on, hold on
    set(handles.axes1,'xlimmode','auto','ylimmode','auto')
    
    % verificar la congruencia de los inputs asociados a la distribucion de la armadura
    if (opcion == 4 && (isempty(gamma) || isempty(rho_w))) || (opcion ~= 5 && isempty(gamma))
        handles.figure1.Pointer = 'arrow';
        handles.datosValidos = false;
        return
    else
        handles.datosValidos = true;
    end
    
    % verificar si se require incluir la rama con P < 0
    incluirTraccion = false;
    if ~isempty(varargin), incluirTraccion = varargin{1}; end
    
    ecuacion = text(0.96,0.95,handles.ecuacion_cuantia{opcion},'units','normalized','interpreter','latex',...
        'horizontalalignment','right','verticalalignment','top','edgecolor','k','backgroundcolor','w');
    
    % inicializar variables que que se usaran para identificar zonas de interes en el diagrama
    ymax = zeros(1,8);
    xlim08Pmax = zeros(1,8);
    xyLim_et0 = zeros(8,2);
    xyLim_et_ey = zeros(8,2);
    xyLim_et0005 = zeros(8,2);
    
    % determinar si la seccion es circular
    Ag = b*h; % area gruesa de la seccion
    factor = 0.8; % factor de reducion de carga maxima por excentricidad accidental para columnas con estribos
    if opcion == 3 && handles.uibuttongroup4.UserData == 2 % columna circular
        opcion = 3.5;
        Ag = 0.25*pi*h^2;
        factor = 0.85; % se supone que se usan espirales para confinar la armadura longitudinal
        ecuacion.String = '$$\rho_g = \frac{A_s}{A_g} = \frac{A_s}{\pi(h/2)^2}$$';
    end
    
    % graficar los diagramas para las cuantias de 1 a 8%
    for i=1:8
        % determinar la distribucion del refuerzo por capas
        [As,ys] = distribuirAs(b,h,opcion,gamma,i*0.01,rho_w);
        
        % determinacion de resultados para los inputs considerados y la cuantia correspondiente
        [M,P,phi_red,et] = interaccionPM2(b,h,As,ys,parametros,n,tol,incluirTraccion,opcion == 3.5);
        
        % transformar M y P a variables adimensionales
        x = phi_red.*M*10^6/(fc*Ag*h);
        y = phi_red.*P*10^3/(fc*Ag);
        
        ymax(i) = y(end);
        
        % incorporar el limite de compresion maxima si el elemento se trata de una columna
        if handles.radiobutton1.Value == 1
            % identificar el limite superior de compresion por concepto de la excentricidad accidental
            Pmax = (0.85*fc*(Ag-sum(As))+fy*sum(As))*10^-3;
            [~,ind] = min(abs(P-factor*Pmax));
            ymax(i) = y(ind);
            xlim08Pmax(i) = x(ind);

            % graficar curva reducida por el factor phi
            plot([x(1:ind) 0],[y(1:ind) y(ind)],'color','b','userdata',i)
%             plot(x(ind:end),y(ind:end),'--','color','b')
        else
            plot(x,y,'color','b','userdata',i)
        end
        
        % determinar zonas limites en el diagrama
        [~,ind_fs0] = min(abs(et)); % tension nula en la capa mas cercana a la fibra mas traccionada
        [~,ind_inicio_fluencia] = min(abs(-et-fy/200000)); % et es negativa para traccion
        [~,ind_et0005] = min(abs(-et-0.005));        
        
        xyLim_et0(i,:) = [x(ind_fs0) y(ind_fs0)];
        xyLim_et_ey(i,:) = [x(ind_inicio_fluencia) y(ind_inicio_fluencia)];
        xyLim_et0005(i,:) = [x(ind_et0005) y(ind_et0005)];
    end
    
    % mostrar las zonas limites en que fs = 0 y fs = fy para el acero en traccion 
    % (se deja oculto el limite en et > 0.005 i.e phi = 0.9)
    plot(xyLim_et0(:,1),xyLim_et0(:,2),'--','color','k')
    plot(xyLim_et_ey(:,1),xyLim_et_ey(:,2),'--','color','k','tag','lim_phi_065')
    plot(xyLim_et0005(:,1),xyLim_et0005(:,2),'--','color','k','tag','lim_phi_09','visible','off')
    
    tipoTexto = {'verticalalignment','bottom','horizontalalignment','left','fontsize',7.5};
    text(xyLim_et0(end,1),xyLim_et0(end,2),'f_s = 0',tipoTexto{:});
    text(xyLim_et_ey(end,1),xyLim_et_ey(end,2),'f_s = f_y',tipoTexto{:});
    
    % mostrar textos indicando las cuantias de acero extremas del diagrama
    % (1 y 8% del area gruesa de la seccion completa o del elemento de borde segun corresponda)
    if handles.radiobutton1.Value == 1        
        handles.ylim_08Pmax = ymax;
        handles.xlim_08Pmax = xlim08Pmax;
%         ylimite = ylim;
%         if ymax(8)/ylimite(2) > 0.9, ylim([ylimite(1) 1.1*ymax(8)]); end
        text(0,ymax(1),' \rho_g = 1%','verticalalignment','top',tipoTexto{3:end});
        text(0,ymax(8),' \rho_g = 8%',tipoTexto{:});        
    else        
        ratioXYZ = handles.axes1.DataAspectRatio;
        alfa = -180/pi*atan(ratioXYZ(1)/ratioXYZ(2)*(y(end)-y(ind_fs0))/x(ind_fs0));
        dy = 0.45;
        textoInicial = ' \rho_g = ';
        if handles.radiobutton5.Value == 1, textoInicial = ' \rho = '; dy = 0.4; end
        text(xyLim_et0(1,1)/2,ymax(1)-0.45*(ymax(1)-xyLim_et0(1,2)),[textoInicial,'1%'],...
            'verticalalignment','top','horizontalalignment','center','fontsize',7.5,'rotation',alfa);
        text(xyLim_et0(8,1)/2,ymax(8)-dy*(ymax(8)-xyLim_et0(8,2)),[textoInicial,'8%'],...
            'verticalalignment','bottom','horizontalalignment','center','fontsize',7.5,'rotation',alfa);
    end
    
    handles.ymin = y(1);
else
    handles.datosValidos = false;
end

handles.figure1.Pointer = 'arrow';
guidata(handles.figure1,handles)
