function calcularAs(handles)
% Esta funci�n determina la armadura dada la geometr�a y esfuerzos en el elemento

% Licenciado bajos los t�rminos del MIT.
% Copyright (c) 2019 Pablo Baez R.

handles.figure1.Pointer = 'watch';

% borrar resultados antiguos
borrarResultadosAs(handles)

% verificar que los inputs sean v�lidos
if handles.datosValidos    
    % obtener valores almacenados de las variables
    opcion = handles.uibuttongroup1.UserData; % 1, 2, 3, 4 � 5
    fc = handles.input_fc.UserData;
    fy = handles.input_fy.UserData;
    gamma = handles.input_gamma.UserData;
    rho_w = handles.input_rho_w.UserData;
    n = handles.input_n.UserData;
    tol = handles.input_error.UserData;    
    
    b = handles.input_b.UserData;
    h = handles.input_h.UserData;
    Mu = handles.input_Mu.UserData;
    Pu = handles.input_Pu.UserData;
    
    if ~isempty(b) && ~isempty(h) && ~isempty(Mu) && ~isempty(Pu)        
        axes(handles.axes1)
        
        % si la carga axial es negativa (tracci�n), trazar nuevamente los diagramas para incluir la rama con P<0
        % (por defecto solo se incluye la flexo-compresi�n)
        incluirTraccion = false;
        ylims = ylim;
        if Pu<0            
            incluirTraccion = true;
            if handles.ymin >= 0
                graficarDiagramas(handles,incluirTraccion);
                handles.ymin = min(get(findobj(handles.axes1,'userdata',8),'ydata'));
            end
            ylim([handles.ymin ylims(2)])
        elseif handles.ymin < 0
            ylim([0 ylims(2)])
        end 
        
        % transformar Mu y Pu a variables adimensionales e incorporar al gr�fico
        xu = 10^6*Mu/(fc*b*h^2);
        yu = 10^3*Pu/(fc*b*h);        
        plot(xu,yu,'o','markerfacecolor','r','tag','punto')
        
        % determinar cual es la curva m�s cercana al par (xu,yu)
        x = cell(8,1);
        y = cell(8,1);
        dist = zeros(1,8);
        for i = 1:8 % cuant�as de referencia de 1 a 8%
            curvas = findobj(handles.axes1,'userdata',i);
            x{i} = curvas.XData;
            y{i} = curvas.YData;
            dist(i) = min(sqrt((x{i}-xu).^2+(y{i}-yu).^2));
        end
        [~,ind1] = min(dist);

        l = h; % largo de referencia para el c�lculo de la cuant�a de acero
        if opcion == 4, l = h*(1-gamma); end % para muros con elementos de borde
        
        % determinar la cuant�a y �rea requerida de acero para las dimensiones y esfuerzos especificados
        % valor que indica si el par (xu,yu) est� bajo la curva definida por ind1
        flag = inpolygon(xu,yu,[0 x{ind1} 0],[0 y{ind1} 0]);
        if ind1 == 1 && flag % no calcular si rho es menor a 1%
            handles.output_rho.setText('< 1.00');
            handles.output_As.setText(['< ',num2str(0.01*b*l*0.01,'%.2f')]);
        elseif ind1 == 8 && ~flag % no calcular ni graficar si rho es mayor a 8%
            delete(findobj(handles.axes1,'tag','punto'));
            handles.output_rho.setText('> 8.00');
            handles.output_As.setText(['> ',num2str(0.08*b*l*0.01,'%.2f')]);
        else
            % determinar en que zona se encuentra el par (xu,yu)
            xlims = xlim;            
            xmax = xlims(2);
            ymin = min(y{8});            
            lim_phi_09 = findobj(handles.axes1,'tag','lim_phi_09');
            lim_phi_065 = findobj(handles.axes1,'tag','lim_phi_065');
            if inpolygon(xu,yu,[0 lim_phi_09.XData xmax xmax 0],...
                    [lim_phi_09.YData(1) lim_phi_09.YData lim_phi_09.YData(end) ymin ymin])
                caso = 1; % control por tracci�n
                phi = 0.9;
            elseif inpolygon(xu,yu,[0 lim_phi_065.XData xmax xmax 0],...
                    [lim_phi_065.YData(1) lim_phi_065.YData lim_phi_065.YData(end) ymin ymin])
                caso = 2; % zona de transici�n (0.65 < phi < 0.9)
                phi = 0.9:-0.005:0.65;
            else
                caso = 3; % control por compresi�n
                phi = 0.65;
            end            
            
            % iterar para calcular rho y su diagrama correspondiente
            parametros = getappdata(handles.figure1,'parametrosMateriales');
            if caso ~=3 || handles.radiobutton1.Value == 0 || ~inpolygon(xu,yu,[0 handles.xlim_08Pmax 0],...
                    [handles.ylim_08Pmax(1) handles.ylim_08Pmax handles.ylim_08Pmax(end)])
                
                % definir cargas axiales nominales para calcular los diagramas PM                
                Pitera = Pu./phi;
                
                % para phi = 0.65 la funci�n 'interaccionPM2' puede requerir 
                % un mayor n�mero de cargas de muestra para arrojar resultados
                if caso == 3, Pitera = linspace(0,Pu/phi,10); end
                
                ind = length(Pitera);
                error = tol+1;
                rho = 0.01*ind1; % se itera a partir de la cuant�a asociada a la curva m�s cercana al par (xu,yu)
                dAs = 0; % diferencial de �rea de acero para el inicio de la iteraci�n
                
                % se usa el mismo error definido para el an�lisis seccional
                % (aunque aquel est� asociado a la sumatoria de cargas axiales)
                while abs(error) > tol
                    rho = rho+dAs/(b*l);
                    [As,ys] = distribuirAs(b,h,opcion,gamma,rho,rho_w);
                    [M,P,phi_red] = interaccionPM2(b,h,As,ys,parametros,n,tol,false,Pitera);
                    
                    % determinar el �ndice para el cual Pu = phi*P
                    if caso == 2, [~,ind] = min(abs(Pu-phi_red.*P)); end
                    
                    % si el vector de momentos M fue redimensionado (al no haber convergencia para ciertas cargas con la
                    % funci�n 'interaccionPM2') probar con un vector de cargas m�s "tupido" y reiniciar las iteraciones
                    if ind > length(M)
                        if caso == 3, Pitera = linspace(0,Pu/phi,length(Pitera)+10); end
                        rho = 0.01*ind1;
                        dAs = 0;
                        ind = length(Pitera);
                        continue
                    end
                    
                    % calcular el error asociado al momento y el diferencial de �rea para la nueva iteraci�n
                    error = Mu-phi_red(ind)*M(ind); % error en kN*m
                    dAs = 10^6*error/(fy*0.9*h); % diferencial de �rea en mm^2
                end
                
                % calcular y trazar el diagrama para la cuant�a y As previamente determinadas
                [M,P,phi_red] = interaccionPM2(b,h,As,ys,parametros,n,tol,incluirTraccion);
                xobj = 10^6*phi_red.*M/(fc*b*h^2);
                yobj = 10^3*phi_red.*P/(fc*b*h);
                
                % si el elemento es una columna, identificar el l�mite superior de
                % compresi�n por concepto de la excentricidad accidental
                if handles.radiobutton1.Value == 1
                    Pnmax = (0.85*fc*(b*h-sum(As))+fy*sum(As))*10^-3;
                    [~,ind08Pmax] = min(abs(P-0.8*Pnmax));                        
                    plot([xobj(1:ind08Pmax) 0],[yobj(1:ind08Pmax) yobj(ind08Pmax)],'r','tag','solucion')
                else
                    plot(xobj,yobj,'r','tag','solucion')
                end
            else % para columnas con Pu = 0.8*phi*Pnmax no es necesario iterar
                rho = ((10^3*Pu/0.8/0.65 -0.85*fc*b*h)/(fy-0.85*fc))/(b*h);
                [As,ys] = distribuirAs(b,h,opcion,gamma,rho,rho_w);
                [M,P,phi_red] = interaccionPM2(b,h,As,ys,parametros,n,tol,incluirTraccion);                
                
                xobj = 10^6*phi_red.*M/(fc*b*h^2);
                yobj = 10^3*phi_red.*P/(fc*b*h);
                
                Pnmax = (0.85*fc*(b*h-sum(As))+fy*sum(As))*10^-3;
                [~,ind08Pmax] = min(abs(P-0.8*Pnmax));
                plot([xobj(1:ind08Pmax) 0],[yobj(1:ind08Pmax) yobj(ind08Pmax)],'r','tag','solucion')
            end

            % actualizaci�n de outputs (objetos Java de tipo JLabels)
            handles.output_rho.setText(num2str(rho*100,'%.2f'));
            handles.output_As.setText(num2str(rho*b*l*0.01,'%.2f'));
        end
    end
end
handles.figure1.Pointer = 'arrow';
