function varargout = abacosPM(varargin)
% ABACOSPM es una interfaz de usuario que permite calcular los diagramas
% de interaccion P-M de una seccion rectangular o circular de hormigon armado y
% determinar el area requerida de acero en funcion de la geometria y los
% esfuerzos ultimos.

% Licenciado bajos los terminos del MIT.
% Copyright (c) 2019 Pablo Baez R.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @abacosPM_OpeningFcn, ...
                   'gui_OutputFcn',  @abacosPM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function abacosPM_OpeningFcn(hObject, ~, handles, varargin)

movegui(handles.figure1,'center')

path = pwd;
addpath(path)
if isdeployed, path = fileparts(mfilename('fullpath')); end
addpath(fullfile(path,'lib'))

axes(handles.axes1)
xlabel('$$\mathbf{\frac{\phi M_n}{f''_c A_g h}}$$','interpreter','latex')
ylabel('$$\mathbf{\frac{\phi P_n}{f''_c A_g}}$$','interpreter','latex','rotation',360,...
    'horizontalalignment','right','verticalalignment','baseline')

handles.ecuacion_cuantia = {'$$\rho_g = \frac{A_s}{A_g} = \frac{A_s}{b\times h}$$'
    '$$\rho_g = \frac{A_s}{A_g} = \frac{A_s}{b\times h}$$'
    '$$\rho_g = \frac{A_s}{A_g} = \frac{A_s}{b\times h}$$'
    '$$\rho = \frac{A_s}{b\times L_b} = \frac{A_s}{b\times L\times(1-\gamma)}$$'
    '$$\rho_g = \frac{A_s}{A_g} = \frac{A_s}{b\times L}$$'};

% definicion de variable que indica si los datos ingresados son validos
handles.datosValidos = true;

% creacion de labels (los "static texts" de Matlab no permiten el uso de html...)
textos(1) = javacomponent(javax.swing.JLabel('<html>f''<sub>c'),[10 45 23 15],handles.uipanel1);
textos(2) = javacomponent(javax.swing.JLabel('<html>f<sub>y'),[10 15 23 15],handles.uipanel1);

[aux1,handles.texto_gamma] = javacomponent(javax.swing.JLabel('<html>&gamma;'),[10 96 23 15],handles.uibuttongroup3);
[aux2,handles.texto_rho_w] = javacomponent(javax.swing.JLabel('<html>&rho;<sub>w'),[10 85 23 15],handles.uibuttongroup3);
handles.texto_rho_w.Visible = 'off';
textos(3) = aux1;
textos(4)= aux2;

textos(5) = javacomponent(javax.swing.JLabel('<html>&#1013;'),[10 17 20 15],handles.uipanel2);

textos(6) = javacomponent(javax.swing.JLabel('<html>M<sub>u'),[10 107 23 15],handles.uipanel3);
textos(7) = javacomponent(javax.swing.JLabel('<html>P<sub>u'),[10 77 23 15],handles.uipanel3);
textos(8) = javacomponent(javax.swing.JLabel('<html><i>&rho;</i><sub>g'),[10 46 23 20],handles.uipanel3);
textos(9) = javacomponent(javax.swing.JLabel('<html>A<sub>s'),[10 17 23 15],handles.uipanel3);
handles.texto_rho = textos(8);

% estandarizacion de fuente (segoe ui de 11px)
for i = 1:length(textos), textos(i).setFont(java.awt.Font('segoe ui',0,11)); end

% creacion y personalizacion de los outputs
outputs = {'output_rho','output_As'};
tooltips = {[] '<html>area total de acero requerido<br>(dividir por 2 para obtener lo de cada borde)'};
pos = [33 45 49 21];
for i = 1:2
    handles.(outputs{i}) = javacomponent(javax.swing.JLabel(''),pos,handles.uipanel3);
    handles.(outputs{i}).setBackground(java.awt.Color(0.94,0.94,0.94));
    handles.(outputs{i}).setBorder(javax.swing.BorderFactory.createLineBorder(java.awt.Color(0.67,0.68,0.7),1));
    handles.(outputs{i}).setFont(java.awt.Font('Segoe ui',0,11));
    handles.(outputs{i}).setHorizontalAlignment(0)
    handles.(outputs{i}).setToolTipText(tooltips{i});
    pos = pos-[0 30 0 0];
end

% definicion de parametros iniciales para las curvas tension-deformacion de los materiales
parametros.modeloAcero = {1 'Elastoplastico'};
parametros.fy = 420;
parametros.ef = 0.15;
parametros.Es2 = 0;
parametros.modeloHormigon = {1 'Saenz'};
parametros.fc = 30;
parametros.e0 = 0.002;
parametros.eu = 0.004;
parametros.modeloHormigonTrac = {1 'sin traccion'};
setappdata(handles.figure1,'parametrosMateriales',parametros); %handles.parametros = parametros;

% creacion de figuras explicativas de los inputs
axes(handles.axes2) % columna
tipoMarcador = {'o','color','k','markerfacecolor','k','tag','fierros'};
grid on, hold on, axis off, axis equal
plot([0.03 2 2 0.03 0.03],[0 0 2 2 0],'color','k','linewidth',1.5)
plot(0.33+1.33/3*[0 1 2 3 0 1 2 3],[0.33*ones(1,4) 1.66*ones(1,4)],tipoMarcador{:},'markersize',4)

plot([-0.7 0],[0 0],'color',0.7*[1 1 1])
plot([-0.7 0],[2 2],'color',0.7*[1 1 1])
plot([1.7 2.7],[1.66 1.66],'color',0.7*[1 1 1])
plot([1.7 2.7],[0.33 0.33],'color',0.7*[1 1 1])
plot([0 0],[-0.5 0],'color',0.7*[1 1 1])
plot([2 2],[-0.5 0],'color',0.7*[1 1 1])

xlim([-0.6 3])
ylim([-0.7 2])

text(1,-0.7,'b','fontsize',7.5,'horizontalalignment','center')
text(-0.77,1.1,'h','fontsize',7.5)
text(2.55,1.1,'\gammah','fontsize',7.5)
text(1,1,'A_s','fontsize',6.5,'horizontalalignment','center','tag','texto_As')
text(3.3,1.1,')','fontsize',10)

tipoFlecha = {'head1width',3,'head2width',3,'head1length',3,'head2length',3,'headstyle','vback3'};
annotation(handles.uipanel5,'doublearrow',[0.195 0.195],[0.22 0.74],tipoFlecha{:})
annotation(handles.uipanel5,'doublearrow',[0.745 0.745],[0.31 0.66],tipoFlecha{:})
annotation(handles.uipanel5,'doublearrow',[0.27 0.67],[0.12 0.12],tipoFlecha{:})
annotation(handles.uipanel5,'arrow',[0.92 0.9],[0.55 0.57],'headwidth',5,'headlength',3,'headstyle','vback3')

axes(handles.axes7) % columna circular
grid on, hold on, axis off, axis equal

plot([-0.7 1],-0.03*[1 1],'color',0.7*[1 1 1])
plot([-0.7 1],[2 2],'color',0.7*[1 1 1])
plot([1 2.7],[1.66 1.66],'color',0.7*[1 1 1])
plot([1 2.7],[0.33 0.33],'color',0.7*[1 1 1])

alfa = linspace(0,2*pi);
alfa2 = linspace(0,2*pi,9);
plot(1+cos(alfa),1+sin(alfa),'color','k','linewidth',1.5)
plot(1+0.665*cos(alfa2),1+0.665*sin(alfa2),'o','color','k','markerfacecolor','k','markersize',4)

xlim([-0.6 3])
ylim([-0.7 2])

text(-0.77,1.1,'h','fontsize',7.5)
text(2.55,1.1,'\gammah','fontsize',7.5)
text(1,1,'A_s','fontsize',6.5,'horizontalalignment','center')
text(3.3,1.1,')','fontsize',10)

annotation(handles.uipanel8,'doublearrow',[0.195 0.195],[0.22 0.74],tipoFlecha{:})
annotation(handles.uipanel8,'doublearrow',[0.745 0.745],[0.31 0.66],tipoFlecha{:})
annotation(handles.uipanel8,'arrow',[0.92 0.9],[0.55 0.57],'headwidth',5,'headlength',3,'headstyle','vback3')

axes(handles.axes5) % muro
grid on, hold on, axis off, axis equal
plot([0 1 1 0 0],[0 0 0.2 0.2 0],'color','k','linewidth',1.5)
plot(0.075*[1 1],[-0.1 0],'color',0.7*[1 1 1],'tag','borde_As')
plot(0.925*[1 1],[-0.1 0],'color',0.7*[1 1 1],'tag','borde_As')

plot([0 0],[0.22 0.35],'color',0.7*[1 1 1])
plot([1 1],[0.22 0.35],'color',0.7*[1 1 1])

plot([0.15 0.15],[0 0.2],'color','k','tag','borde_As')
plot([0.85 0.85],[0 0.2],'color','k','tag','borde_As')

plot([0.05 0.1 0.1 0.05 0.9 0.95 0.95 0.9],[0.05 0.05 0.15 0.15 0.05 0.05 0.15 0.15],tipoMarcador{:},'markersize',3)
plot(0.175+0.05*(0:13),0.05*ones(1,14),tipoMarcador{:},'markersize',1)
plot(0.175+0.05*(0:13),0.15*ones(1,14),tipoMarcador{:},'markersize',1)

xlim([-0.15 1.15])
ylim([-0.2 0.45])

text(0.5,0.35,'L','fontsize',7.5,'horizontalalignment','center')
text(0.5,-0.15,'\gammaL','fontsize',7.5,'horizontalalignment','center','tag','texto_gamma')
text(2.9,1.1,')','fontsize',10)
text(-0.1,0.1,'A_s','fontsize',6.5,'horizontalalignment','center','tag','texto_As')
text(1.1,0.1,'A_s','fontsize',6.5,'horizontalalignment','center','tag','texto_As')
text(0.5,0.42,')','fontsize',10,'rotation',90,'horizontalalignment','center','verticalalignment','baseline')

handles.flecha_gamma = annotation(handles.uipanel7,'doublearrow',[0.18 0.82],[0.17 0.17],tipoFlecha{:});
annotation(handles.uipanel7,'doublearrow',[0.12 0.88],[0.55 0.55],tipoFlecha{:})
annotation(handles.uipanel7,'arrow',[0.45 0.43],[0.68 0.66],'headwidth',5,'headlength',3,'headstyle','vback3')

% crear un boton y asociarlo a la definicion de las leyes constitutivas de los materiales
[handles.boton,~]=javacomponent(javax.swing.JToggleButton('Leyes Constitutivas'),[7 127 125 25],handles.uipanel2);
handles.boton.setFont(java.awt.Font('Segoe ui',0,11));
handles.boton.setMargin(java.awt.Insets(0,0,0,0));
cargarInterfazMateriales(handles,path);

% Choose default command line output for abacosPM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% iniciar calculando los diagramas para los valores por defecto
graficarDiagramas(handles)

function varargout = abacosPM_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

function input_gamma_Callback(hObject, ~, handles) %#ok<*DEFNU>
verificarInput(hObject,handles,1,@(x)x>0 && x<=1,...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un real entre 0 y 1')

function input_rho_w_Callback(hObject, ~, handles)
verificarInput(hObject,handles,1,@(x)x>=0,...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un real no negativo.')

function input_n_Callback(hObject, ~, handles)
verificarInput(hObject,handles,1,@(x)rem(x,1) == 0 && x > 0,...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un numero entero positivo.')

function input_error_Callback(hObject, ~, handles)
verificarInput(hObject,handles,10,@(x)x>0,...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un real positivo.')

function input_b_Callback(hObject, ~, handles)
verificarInput2(hObject,handles,10,@(x)x>0,...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un real positivo.')

function input_h_Callback(hObject, ~, handles)
verificarInput2(hObject,handles,10,@(x)x>0,...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un real positivo.')

function input_Mu_Callback(hObject, ~, handles)
verificarInput2(hObject,handles,10,@(x)x>=0,...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un real positivo.')

function input_Pu_Callback(hObject, ~, handles)
verificarInput2(hObject,handles,10,@(x)isreal(x),...
    'El valor ingresado (directamente o mediante el uso de funciones) debe ser un numero real.')

function verificarInput(hObject,handles,factorConversion,condicion,mensajeError)
hObject.UserData = [];
try
    input = evalin('base',hObject.String);
    hObject.String = num2str(input);
    if feval(condicion,input)        
        hObject.UserData = factorConversion*input;
        graficarDiagramas(handles) % recalcular el diagrama
    else
        axes(handles.axes1), cla 
        errordlg(mensajeError,'Error en el ingreso de datos','modal')
    end
catch
    axes(handles.axes1), cla
    errordlg(mensajeError,'Error en el ingreso de datos','modal')
end

function verificarInput2(hObject,handles,factorConversion,condicion,mensajeError)
hObject.UserData = [];
try
    input = evalin('base',hObject.String);
    hObject.String = num2str(input);
    if feval(condicion,input)
        hObject.UserData = factorConversion*input;        
        calcularAs(handles)
    else % borrar resultados antiguos 
        borrarResultadosAs(handles)
        errordlg(mensajeError,'Error en el ingreso de datos','modal')
    end
catch
    borrarResultadosAs(handles)
    errordlg(mensajeError,'Error en el ingreso de datos','modal')
end

function input_fc_Callback(hObject, ~, handles)
fc = [16 20:5:55]; % valores de f'c en MPa mas usuales
hObject.UserData = fc(hObject.Value);
setappdata(handles.figure1,'parametrosMateriales',...
    setfield(getappdata(handles.figure1,'parametrosMateriales'),'fc',hObject.UserData))
graficarDiagramas(handles)

function input_fy_Callback(hObject, ~, handles)
fy = [280 420]; % valores de fy en MPa tipicos para el acero utilizado en hormigon armado
hObject.UserData = fy(hObject.Value);
setappdata(handles.figure1,'parametrosMateriales',...
    setfield(getappdata(handles.figure1,'parametrosMateriales'),'fy',hObject.UserData))
graficarDiagramas(handles)

function radiobutton1_KeyPressFcn(~, eventdata, handles)
if ~isempty(strfind(eventdata.Key,'arrow'))
    handles.radiobutton2.Value = 1;
    uibuttongroup1_SelectionChangedFcn([],[],handles)
end

function radiobutton2_KeyPressFcn(~, eventdata, handles)
if ~isempty(strfind(eventdata.Key,'arrow'))
    handles.radiobutton1.Value = 1;
    uibuttongroup1_SelectionChangedFcn([],[],handles)
end

function uibuttongroup1_SelectionChangedFcn(~, ~, handles)
% dejar por defecto siempre la primera sub-opcion
handles.radiobutton5.Value = 1;
% borrar inputs asociados al calculo de As
inputs = {'input_b','input_h','input_Mu','input_Pu'};
for i = 1:4
    handles.(inputs{i}).String = '';
    handles.(inputs{i}).UserData = [];
end
uibuttongroup3_SelectionChangedFcn(handles.radiobutton5, [], handles)

function uibuttongroup3_SelectionChangedFcn(hObject, ~, handles)

% determinar si se escogio la opcion 'Columna'
flag = handles.radiobutton1.Value == 1;

% mostrar y ocultar elementos segun corresponda
flags={'off' 'on'};
handles.uibuttongroup4.Visible = 'off';
handles.input_b.Enable = 'on';
handles.input_gamma.Visible = 'on';
handles.texto_rho_w.Visible = 'off';
handles.texto_gamma.Visible = 'on';
handles.radiobutton7.Visible = flags{flag+1};
handles.uipanel5.Visible = flags{flag+1};
handles.uipanel7.Visible = flags{~flag+1};
handles.uipanel8.Visible = 'off';
handles.input_rho_w.Visible = flags{~flag+1};
handles.input_As2.Visible = flags{~flag+1};
handles.texto_rho.setText('<html><i>&rho;</i><sub>g')
handles.output_rho.setToolTipText([])
handles.output_As.setToolTipText('<html>area total de acero requerido<br>(dividir por 2 para obtener lo de un borde)')

sub_opcion = hObject.UserData;
if flag % columna
    handles.uibuttongroup1.UserData = sub_opcion;
    handles.radiobutton6.String = 'Lateral';
    handles.text9.String = 'h';
    delete(findobj(handles.axes2,'type','line','tag','fierros'))
    set(findobj(handles.axes2,'type','text','tag','texto_As'),'visible','on')
    handles.texto_gamma.Position = [10 96 23 15];
    handles.input_gamma.Position = [33 93 49 21];    
    
    % trazar esquema explicativo segun la sub opcion escogida
    axes(handles.axes2)
    tipoMarcador1 = {'o','color','k','markerfacecolor','k','markersize',4,'tag','fierros'};
    if sub_opcion == 1
        plot(0.33+1.33/3*[0 1 2 3 0 1 2 3],[0.33*[1 1 1 1] 1.66*[1 1 1 1]],tipoMarcador1{:})
    elseif sub_opcion == 2
        plot([0.33*[1 1 1 1] 1.66*[1 1 1 1]],0.33+1.33/3*[0 1 2 3 0 1 2 3],tipoMarcador1{:})
    else%if sub_opcion == 3
    handles.uibuttongroup4.Visible = 'on';
        handles.uibuttongroup4.UserData = 1;
        handles.checkbox3.Value = 1;
        handles.radiobutton12.Value = 0;
        handles.output_As.setToolTipText('area total de acero requerido')
        plot(0.33+1.33/3*[0 1 2 3 0 1 2 3],[0.33*[1 1 1 1] 1.66*[1 1 1 1]],tipoMarcador1{:})
        plot([0.33*[1 1] 1.66*[1 1]],0.33+1.33/3*[1 2 1 2],tipoMarcador1{:})        
    end

else% handles.radiobutton2.Value == 1 % muro
    handles.uibuttongroup1.UserData = sub_opcion+3;
    handles.radiobutton6.String = 'Uniforme';
    handles.text9.String = 'L';    
    delete(findobj(handles.axes5,'type','line','tag','fierros'))
    handles.texto_gamma.Position = [10 118 23 15];
    handles.input_gamma.Position = [33 115 49 21];
    handles.flecha_gamma.Visible = flags{3-sub_opcion};
    set(findobj(handles.axes5,'type','text','tag','texto_gamma'),'visible',flags{3-sub_opcion});
    
    textos_As = findobj(handles.axes5,'type','text','tag','texto_As');
    for i = 1:2, textos_As(i).Visible = flags{3-sub_opcion}; end
    
    lineas = findobj(handles.axes5,'type','line','tag','borde_As');
    for i = 1:4, lineas(i).Visible = flags{3-sub_opcion}; end
    
    % trazar esquema explicativo segun la sub opcion escogida    
    axes(handles.axes5)
    tipoMarcador2 = {'o','color','k','markerfacecolor','k','markersize',1,'tag','fierros'};
    if sub_opcion == 1
        handles.texto_rho_w.Visible = 'on';
        plot([0.05 0.1 0.1 0.05 0.9 0.95 0.95 0.9],[0.05 0.05 0.15 0.15 0.05 0.05 0.15 0.15],tipoMarcador2{:},'markersize',3)
        plot(0.175+0.05*(0:13),0.05*ones(1,14),tipoMarcador2{:})
        plot(0.175+0.05*(0:13),0.15*ones(1,14),tipoMarcador2{:})
        handles.texto_rho.setText('<html><i>&rho;')
        handles.output_rho.setToolTipText('cuantia del elemento de borde')
        handles.output_As.setToolTipText('area de acero requerida para un solo borde')
    elseif sub_opcion == 2
        plot(0.175+0.05*(-3:16),0.05*ones(1,20),tipoMarcador2{:})
        plot(0.175+0.05*(-3:16),0.15*ones(1,20),tipoMarcador2{:})
        
        handles.texto_gamma.Visible = 'off';
        handles.input_gamma.Visible = 'off';
        handles.input_rho_w.Visible = 'off';
    end
end
graficarDiagramas(handles) % recalcular diagramas

% opcion columna rectangular
function checkbox3_Callback(hObject, ~, handles)
hObject.Value = 1;
handles.radiobutton12.Value = 0;
handles.uipanel5.Visible = 'on';
handles.uipanel8.Visible = 'off';
if handles.uibuttongroup4.UserData == 2
    handles.input_b.Enable = 'on';
    handles.uibuttongroup4.UserData = 1;    
    graficarDiagramas(handles)
end

% opcion columna circular
function radiobutton12_Callback(~, ~, handles)
handles.checkbox3.Value = 0;
handles.uipanel5.Visible = 'off';
handles.uipanel8.Visible = 'on';
if handles.uibuttongroup4.UserData == 1
    handles.input_b.Enable = 'off';
    handles.uibuttongroup4.UserData = 2;    
    graficarDiagramas(handles)
end
