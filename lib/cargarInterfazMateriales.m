function cargarInterfazMateriales(handles,path)
% Esta función crea la interfaz para ver y editar las leyes constitutivas
% del acero y el hormigón, con las cuales se calculan los diagramas de
% interacción PM.

% Licenciado bajos los términos del MIT.
% Copyright (c) 2019 Pablo Baez R.

% cargar el directorio donde se encuentra la clase 'JavaAxes'
folder = fullfile(path,'lib');
if ~ismember(folder,javaclasspath('-all')) && isempty(which('JavaAxes.class'))
    if ~isempty(javaclasspath('-dynamic'))
        warning(['Limpiar las clases Java anteriormente agregadas o ',...
            'crear un ''classpath'' personalizado incluyendo el directorio ',folder]);
        return
    else
        javaaddpath(folder);
    end
end

font = java.awt.Font('segoe ui',0,11);

% crear popup que servirá como una ventana anidada a la figura base
popup = javaObjectEDT('javax.swing.JPopupMenu');
popup = handle(popup,'CallbackProperties');
popup.setPreferredSize(java.awt.Dimension(365,365));

% crear un panel con pestañas para cada material (acero y hormigón)
tabgp = javax.swing.JTabbedPane;
popup.add(tabgp);

%% acero
% crear panel que será el contenedor para la pestaña asociada a la curva del acero
panelTab1 = javax.swing.JPanel;
panelTab1.setLayout(javax.swing.BoxLayout(panelTab1,javax.swing.BoxLayout.Y_AXIS));

% crear un sub-panel para las opciones de modelos e inputs
panel_top = javax.swing.JPanel;
panel_top.setLayout(java.awt.GridBagLayout);
panel_top.setPreferredSize(java.awt.Dimension(365,110));

% crear opciones para definir la ley constitutiva del acero
radiopanel = javax.swing.JPanel;
radiopanel.setLayout(java.awt.GridLayout(3,1));
grupo = javax.swing.ButtonGroup;

opcion = {'Elastoplástico' 'Mander'};
for i = 1:2
    opcionAcero = javaObjectEDT('javax.swing.JRadioButton',opcion{i});    
    opcionAcero.setFont(font);
    opcionAcero.setActionCommand(num2str(i));
    grupo.add(opcionAcero);
    radiopanel.add(opcionAcero);
    handles.modelosAcero(i) = handle(opcionAcero,'CallbackProperties');
end
handles.modelosAcero(1).setSelected(true); % dejar el modelo elastoplástico como opción por defecto
handles.opcionesAcero = grupo;

grid_top = java.awt.GridBagConstraints;
grid_top.insets = java.awt.Insets(10,0,0,10);
grid_top.anchor = java.awt.GridBagConstraints.WEST;
grid_top.weightx = 1;
grid_top.weighty = 1;
panel_top.add(radiopanel,grid_top);

% crear un sub-sub-panel para aglutinar los inputs asociados al modelo elastoplástico
panel_op1 = javax.swing.JPanel;
panel_op1.setLayout(java.awt.GridBagLayout);
panel_op1.setBorder(javax.swing.border.EmptyBorder(0,0,28,105))
grid_op1= java.awt.GridBagConstraints;
grid_op1.anchor = java.awt.GridBagConstraints.NORTHWEST;
grid_op1.insets = java.awt.Insets(8,0,0,0);

% crear un sub-sub-panel para aglutinar los inputs asociados al modelo de Mander
panel_op2 = javax.swing.JPanel;
panel_op2.setLayout(java.awt.GridBagLayout);
grid_op2 = java.awt.GridBagConstraints;
grid_op2.anchor = java.awt.GridBagConstraints.NORTHWEST;
grid_op2.insets = java.awt.Insets(8,0,0,0);

% crear los inputs y sus textos asociados y agregarlos a los sub-sub-paneles correspondientes
textos1 = {'<html>&epsilon;<sub>f','<html>E<sub>s2',...
    '<html>&epsilon;<sub>sh','<html>&epsilon;<sub>su','<html>&epsilon;<sub>f','<html>E<sub>sh','<html>f<sub>su'};
textos2 = {'','<html>&nbsp; &times; E<sub>s',...
    '','','','<html>&nbsp; &times; E<sub>s','<html>&nbsp; &times; f<sub>y',};
valores_inputs = [0.15 0.01 0.01 0.1 0.15 0.033 1.5];
for i = 1:7
    texto = javax.swing.JLabel(textos1{i});
    input = javaObjectEDT('javax.swing.JTextField',num2str(valores_inputs(i)));    
    texto2 = javax.swing.JLabel(textos2{i});
    
    texto.setPreferredSize(java.awt.Dimension(25,21));
    input.setPreferredSize(java.awt.Dimension(49,21));
    texto2.setPreferredSize(java.awt.Dimension(30,21));
    
    input.setHorizontalAlignment(0) % alinear al centro el texto de los inputs
    
    texto.setFont(font);
    input.setFont(font);
    texto2.setFont(font);
    
    if i > 2
        grid = grid_op2;
        pane = panel_op2;
        if i <= 5
            grid.gridy = i-3;
        else
            grid.gridy = i-6;
        end
        input.setFocusTraversalKeysEnabled(false); % para definir un 'tab order' personalizado
    else
        grid = grid_op1;
        pane = panel_op1;
        grid.gridy = i-1;
    end

    pane.add(texto,grid);
    pane.add(input,grid);
    pane.add(texto2,grid);
    
    handles.inputsMateriales(i) = handle(input,'CallbackProperties');
end
panel_op1.setVisible(true)
panel_op2.setVisible(false)

grid_top.anchor = java.awt.GridBagConstraints.NORTHWEST;
grid_top.insets = java.awt.Insets(10,1,-10,10);
panel_top.add(panel_op1,grid_top);
panel_top.add(panel_op2,grid_top);
panelTab1.add(panel_top);

% crear gráfico para mostrar la curva del modelo actualmente seleccionado
% (por defecto se establece el modelo elastoplástico con ey = 420/200000 = 0.0021)
graficoAcero = javaObjectEDT('JavaAxes',[0 0.0021 0.15],[0 1 1]);
graficoAcero.setAxisLabels('<html>&epsilon;','<html>f<sub>s</sub><hr>f<sub>y')
graficoAcero.setPreferredSize(java.awt.Dimension(365,200));
handles.curvaAcero = graficoAcero;
panelTab1.add(graficoAcero);

tabgp.addTab('Acero', panelTab1);

%% hormigón
% crear panel que será el contenedor para la pestaña asociada a la curva del hormigón
panelTab2 = javax.swing.JPanel;
panelTab2.setLayout(javax.swing.BoxLayout(panelTab2,javax.swing.BoxLayout.Y_AXIS));

% crear un sub-panel para las opciones de modelos e inputs
grupoTab2 = javax.swing.ButtonGroup;
panelTopTab2 = javax.swing.JPanel;
panelTopTab2.setLayout(java.awt.GridBagLayout);
panelTopTab2.setPreferredSize(java.awt.Dimension(365,110));

% crear opciones para definir la ley constitutiva del hormigón a compresión
radiopanelTab2=javax.swing.JPanel;
radiopanelTab2.setLayout(java.awt.GridLayout(4,1));

opcion = {'Saenz' 'Hognestad' 'Thorenfeldt 1' 'Thorenfeldt 2'};
tooltips = {[] [] 'calibrado por Collins y Porasz' 'calibrado por Carreira y Kuang-Han'};
for i = 1:4
    opcionHormigon = javaObjectEDT('javax.swing.JRadioButton',opcion{i});    
    opcionHormigon.setFont(font);
    opcionHormigon.setActionCommand(num2str(i));
    opcionHormigon.setToolTipText(tooltips{i});
    grupoTab2.add(opcionHormigon);
    radiopanelTab2.add(opcionHormigon);
    handles.modelosHormigon(i) = handle(opcionHormigon,'CallbackProperties');
end
handles.modelosHormigon(1).setSelected(true); % dejar el modelo elastoplástico como opción por defecto
handles.opcionesHormigon = grupoTab2;

grid2 = java.awt.GridBagConstraints;
grid2.anchor = java.awt.GridBagConstraints.NORTHWEST;
grid2.insets = java.awt.Insets(5,0,0,0);
grid2.weightx = 0.145;
panelTopTab2.add(radiopanelTab2,grid2);

% % crear opciones para definir la ley constitutiva del hormigón a tracción
% grupo2Tab2 = javax.swing.ButtonGroup;
% radiopanel2Tab2 = javax.swing.JPanel;
% radiopanel2Tab2.setLayout(java.awt.GridLayout(2,1));
% 
% opcion = {'<html>Sin resistencia <br>a tracción' '<html>Lineal hasta <br>la rotura'};
% for i = 1:2
%     opcionHormigonTrac = javaObjectEDT('javax.swing.JRadioButton',opcion{i});    
%     opcionHormigonTrac.setFont(font);
%     opcionHormigonTrac.setActionCommand(num2str(i));
%     grupo2Tab2.add(opcionHormigonTrac);
%     radiopanel2Tab2.add(opcionHormigonTrac);
%     handles.modelosHormigonTrac(i) = handle(opcionHormigonTrac,'CallbackProperties');
% end
% handles.modelosHormigonTrac(1).setSelected(true); % por defecto no considerar resistencia a tracción en el hormigón
% handles.opcionesHormigonTrac = grupo2Tab2;
% 
% dejar el modelo sin resistencia a tracción como opción por defecto
% handles.modelosHormigonTrac(1).setSelected(true);
% 
% separador = javax.swing.JLabel;
% separador.setBorder(javax.swing.BorderFactory.createLineBorder(java.awt.Color(0.67,0.68,0.7),1));
% separador.setPreferredSize(java.awt.Dimension(1,90))
% panelTopTab2.add(separador,grid2);
% % panelTopTab2.add(javax.swing.JLabel('<html><hr width=1 size=85></hr>'),grid2);
% panelTopTab2.add(radiopanel2Tab2,grid2);

% crear un sub-sub-panel para aglutinar los inputs asociados a los modelos del hormigón
panelInputsTab2 = javax.swing.JPanel;
panelInputsTab2.setLayout(java.awt.GridBagLayout);
gridInputsTab2= java.awt.GridBagConstraints;
gridInputsTab2.anchor = java.awt.GridBagConstraints.NORTHWEST;
gridInputsTab2.insets = java.awt.Insets(8,0,0,0);

% crear los inputs y sus textos asociados y agregarlos al panel 'panelInputsTab2'
textos = {'<html>&epsilon;<sub>0','<html>&epsilon;<sub>u'};
valoresInputsTab2 = [0.002 0.0038];
for i = 1:2
    gridInputsTab2.gridy = i-1;

    texto = javax.swing.JLabel(textos{i});
    input = javaObjectEDT('javax.swing.JTextField',num2str(valoresInputsTab2(i)));    
    
    texto.setPreferredSize(java.awt.Dimension(25,21));
    input.setPreferredSize(java.awt.Dimension(49,21));
    
    input.setHorizontalAlignment(0)
    
    texto.setFont(font);
    input.setFont(font);
    
    panelInputsTab2.add(texto,gridInputsTab2);
    panelInputsTab2.add(input,gridInputsTab2);

    handles.inputsMateriales(i+7) = handle(input,'CallbackProperties');
end
texto.setVisible(false);
handles.texto_eu = texto;
input.setVisible(false);

grid2.weightx = 1-grid2.weightx;
panelTopTab2.add(panelInputsTab2,grid2);
panelTab2.add(panelTopTab2);

% crear gráfico para mostrar la curva del modelo actualmente seleccionado
% (por defecto se considera el modelo de Saenz con e0 = 0.002)
e = linspace(0,0.004);
graficoHormigon = javaObjectEDT('JavaAxes',e,2*(e/0.002)-(e/0.002).^2);
graficoHormigon.setAxisLabels('<html>&epsilon;','<html>f<sub>c</sub><hr>f''<sub>c')
graficoHormigon.setPreferredSize(java.awt.Dimension(365,205));
panelTab2.add(graficoHormigon);
handles.curvaHormigon = graficoHormigon;

tabgp.addTab('Hormigón', panelTab2);

%% crear botones y asignar callbacks
% crear botones para guardar, cancelar y establecer los valores predetermiandos
panel_botones = javax.swing.JPanel;
panel_botones.setLayout(java.awt.GridLayout(1,3,20,0));

boton_guardar = javaObjectEDT('javax.swing.JButton','Guardar');
boton_cancelar = javaObjectEDT('javax.swing.JButton','Cancelar');
boton_defecto = javaObjectEDT('javax.swing.JButton','Val. Predeterm.');

boton_guardar = handle(boton_guardar,'CallbackProperties');
boton_cancelar = handle(boton_cancelar,'CallbackProperties');
boton_defecto = handle(boton_defecto,'CallbackProperties');

boton_defecto.setMargin(java.awt.Insets(0,0,0,0));
panel_botones.setBorder(javax.swing.border.EmptyBorder(0,20,3,20));

panel_botones.setPreferredSize(java.awt.Dimension(365,30));
panel_botones.add(boton_guardar);
panel_botones.add(boton_cancelar);
panel_botones.add(boton_defecto);

popup.add(panel_botones);

% asignar callbacks a opciones, inputs y botones
for i = 1:2, handles.modelosAcero(i).ActionPerformedCallback = @(~,~)graficarCurvaAcero(handles); end
for i = 1:4, handles.modelosHormigon(i).ActionPerformedCallback = @(~,~)graficarCurvaHormigon(handles); end
% for i = 1:2, handles.modelosHormigonTrac(i).ActionPerformedCallback = @(~,~)graficarCurvaHormigon(handles); end

nTab = 1;
for i = 1:9
    if i > 7, nTab = 2; end 
    handles.inputsMateriales(i).KeyPressedCallback = @(src,evt)recalcularPost_Enter_Tab(src,evt,handles,nTab,i);
    handles.inputsMateriales(i).FocusLostCallback = @(src,evt)recalcularPostEdicion(src,evt,handles,popup,nTab);
    handles.inputsMateriales(i).FocusGainedCallback = @(src,~)seleccionarEditText(src);    
end

boton_guardar.ActionPerformedCallback = @(~,~)guardarParametros(handles,popup);
boton_cancelar.ActionPerformedCallback = @(~,~)popup.setVisible(false);
boton_defecto.MousePressedCallback = @(~,~)resetearValores(handles);

% handles.boton.setComponentPopupMenu(popup);
handles.boton.ActionPerformedCallback = @(src,~)botonLeyes(src,handles,popup,tabgp);
handles.boton.PropertyChangeCallback = @(src,~)src.setSelected(false);
popup.PropertyChangeCallback = @(~,evt)handles.boton.setSelected(~(strcmp(evt.getPropertyName,'ancestor') && ...
    isa(evt.getOldValue,'javax.swing.JPanel') && isempty(evt.getNewValue) && handles.boton.isSelected));

end
%% funciones asoiadas a los callbacks

% función que se ejecuta al presionar el botón 'Leyes Constitutivas'
function botonLeyes(boton,handles,popup,panelTab)
if boton.isSelected
    % mover la ventana principal si no hay suficiente espacio para
    % desplegar la interfaz sin solapar parte de esta
    pantalla = get(0,'screensize');
    posFigura = handles.figure1.Position;
    if pantalla(3) >= 830+365 && posFigura(1)+830 > pantalla(3)-365
        handles.figure1.Position = [pantalla(3)-830-365,posFigura(2:4)];
    end
    panelTab.setSelectedIndex(0)
    popup.show(boton,boton.getWidth+30,0)
    verLeyesConstitutivas(handles) % mostrar los parámetros actualmente vigentes
else
    popup.setVisible(false)
end
end

% función que se ejecuta cuando los inputs pierden el foco
function recalcularPostEdicion(textfield,evt,handles,popup,nTab)
evaluarInput(textfield);
elemento = evt.getOppositeComponent;

% no recalcular si se cambia de modelo, pues esto se hará directamente
% a través del 'ActionPerformedCallback' de los radiobuttons
if ~isa(elemento,'javax.swing.JRadioButton') && popup.isVisible
    if nTab == 1
        graficarCurvaAcero(handles)
    else%if nTab == 2
        graficarCurvaHormigon(handles)
    end
end
end

% función que se ejecuta después de presionar las teclas enter o tab mientras se edita un input
function recalcularPost_Enter_Tab(textfield,evt,handles,nTab,id)
global altPressed; % variable que indica si la tecla Alt fue presionada

% actualizar las curvas al presionar las teclas enter o tab
if ismember(evt.getKeyCode,[9 10])
    altPressed = false;
    evaluarInput(textfield);
    if nTab == 1
        opcion = str2double(handles.opcionesAcero.getSelection.getActionCommand);
        if evt.getKeyCode == 9 && opcion == 2
            inputs = handles.inputsMateriales;
            if id == 7, id = 2; end
            inputs(id+1).requestFocus;
        end        
        graficarCurvaAcero(handles);
    else%if nTab == 2
        graficarCurvaHormigon(handles);
    end
elseif evt.getKeyCode == 18 % alt
    % por defecto los popups pierden el foco y desaparecen al presionar la tecla Alt
    % por lo que hay que volver a hacer visible la interfaz y devolver el foco al textfield
    altPressed = true;    
    popup = textfield.getParent.getParent.getParent.getParent.getParent;
    popup.setVisible(true);    
    textfield.requestFocus;
    textfield.setCaretPosition(textfield.getDocument.getLength)
end
end

% función que se ejecuta cuando un input captura el foco
function seleccionarEditText(textfield)
global altPressed;
if ~altPressed, textfield.selectAll; end
end

% función que verifica que el valor ingresado es los textfields son válidos
function evaluarInput(textfield)
mensajeError = ['El valor ingresado (directamente o mediante el uso de funciones)'...
    'debe ser un real no negativo.'];
try
    input = evalin('base',char(textfield.getText));
    if input >= 0
        textfield.setText(num2str(input));
    else
        errordlg(mensajeError,'Error en el ingreso de datos','modal')
    end
catch
    errordlg(mensajeError,'Error en el ingreso de datos','modal')
end
end

% función que grafica la curva esfuerzo-deformación del acero dependiendo 
% del modelo y los parámetros establecidos
function graficarCurvaAcero(handles)
ey = handles.input_fy.UserData/200000;
inputs = handles.inputsMateriales;

% mostrar los paneles con los inputs correspondientes
opcion = str2double(handles.opcionesAcero.getSelection.getActionCommand);
flag = (opcion == 1);
inputs(1).getParent.setVisible(flag);
inputs(3).getParent.setVisible(~flag);

% graficar la ley constitutiva del acero según la opción y los inputs elegidos
if flag % modelo elastoplástico
    ef = str2double(inputs(1).getText);
    Es2 = str2double(inputs(2).getText);
    
    if ~isnan(ef) && ~isnan(Es2)
        handles.curvaAcero.plot([0 ey ef],[0 1 1+Es2*(ef-ey)/ey]);
    end
else%if opcion == 2 % modelo de Mander
    esh = str2double(inputs(3).getText);
    esu = str2double(inputs(4).getText);
    ef = str2double(inputs(5).getText);
    Esh = str2double(inputs(6).getText);
    fsu = str2double(inputs(7).getText);
    
    if ~isnan(esh) && ~isnan(esu) && ~isnan(ef) && ~isnan(Esh) && ~isnan(fsu)
        p = Esh*(esu-esh)/(fsu-1)/ey;
        e = linspace(esh,ef);
        f = fsu+(1-fsu)*abs((esu-e)/(esu-esh)).^p;
        handles.curvaAcero.plot([0 ey e],[0 1 f]);
    end
end
end

% función que grafica la curva esfuerzo-deformación del hormigón
% dependiendo del modelo y los parámetros establecidos
function graficarCurvaHormigon(handles)
opcion = str2double(handles.opcionesHormigon.getSelection.getActionCommand);
flag = opcion ~= 1;
handles.texto_eu.setVisible(flag);
handles.inputsMateriales(9).setVisible(flag);

% graficar la ley constitutiva del hormigón según la opción y los inputs elegidos
e0 = str2double(handles.inputsMateriales(8).getText);
if ~isnan(e0)
    fc = handles.input_fc.UserData;
    if flag
        eu = str2double(handles.inputsMateriales(9).getText);
        if ~isnan(eu)
            if opcion == 2 % modelo de Hognestad
                e = linspace(0,e0);
                f = 2*e/e0-(e/e0).^2;
                e = [e e0 eu];
                f = [f 1 0.85];                
            else % modelos de Thorenfeldt                
                if opcion == 3 % calibrado según Collins y Porasz
                    r = 0.8+fc/17;
                    k = [ones(1,50),(0.67+fc/62)*ones(1,50)];
                else%if opcionHormigon == 4 % calibrado según Carreira y Kuang-Han
                    r = 1.55+(fc/32.4)^3;
                    k = 1;
                end
                e = [linspace(0,e0,50) linspace(e0,eu,50)];
                f = r*(e/e0)./(r-1+(e/e0).^(r*k));                
            end
        end
    else%if opcion == 1 % modelo de Saenz
        eu = 2*e0;
        e = linspace(0,eu);
        f = 2*(e/e0)-(e/e0).^2;        
    end
    handles.curvaHormigon.plot(e,f)

%   % tracción en el hotmigón? (incide poco y nada en el cálculo y puede provocar problemas de convergencia)
%     opcionHormigonTrac = str2double(handles.opcionesHormigonTrac.getSelection.getActionCommand);
%     if opcionHormigonTrac == 2
%         fcr = 0.62/sqrt(fc); % normalizada por f'c
%         handles.curvaHormigon.plot([-0.62/4700 e],[-fcr f])
%     end
end
end

% función que permite visualizar los últimos modelos establecidos
% (los que fueron guardados a través del botón 'Guardar')
function verLeyesConstitutivas(handles)
parametros = getappdata(handles.figure1,'parametrosMateriales');
inputs = handles.inputsMateriales;

% mostrar los parámetros actuales de la curva esfuerzo-deformación del acero
modeloAcero = parametros.modeloAcero;
handles.modelosAcero(modeloAcero{1}).setSelected(true);
if modeloAcero{1} == 1 % modelo elastoplástico
    inputs(1).setText(num2str(parametros.ef));
    inputs(2).setText(num2str(parametros.Es2));
else%if modeloAcero{1} == 2 % modelo de Mander
    inputs(3).setText(num2str(parametros.esh));
    inputs(4).setText(num2str(parametros.esu));
    inputs(5).setText(num2str(parametros.ef));
    inputs(6).setText(num2str(parametros.Esh));
    inputs(7).setText(num2str(parametros.fsu));
end

% mostrar los parámetros actuales de la curva esfuerzo-deformación del hormigón
modeloHormigon = parametros.modeloHormigon;
handles.modelosHormigon(modeloHormigon{1}).setSelected(true);
inputs(8).setText(num2str(parametros.e0));
if modeloHormigon{1} ~= 1, inputs(9).setText(num2str(parametros.eu)); end

% modeloHormigonTrac = parametros.modeloHormigonTrac;
% handles.modelosHormigonTrac(modeloHormigonTrac{1}).setSelected(true);

graficarCurvaAcero(handles)
graficarCurvaHormigon(handles)
end

% función que se ejecuta al presionar el botón 'Val. Predeterm.'
% y que muestra los valores por defecto de los inputs de los modelos
function resetearValores(handles)
inputs = handles.inputsMateriales;
valores_inputs = {'0.15' '0' '0.01' '0.1' '0.15' '0.033' '1.5' '0.002' '0.0038'};
for i = 1:9, inputs(i).setText(valores_inputs{i}); end

graficarCurvaAcero(handles)
graficarCurvaHormigon(handles)
end

% función que se ejecuta cuando el botón 'Guardar' es presionado
function guardarParametros(handles,popup)
popup.setVisible(false);

opcionAcero = str2double(handles.opcionesAcero.getSelection.getActionCommand);
opcionHormigon = str2double(handles.opcionesHormigon.getSelection.getActionCommand);
% opcionHormigonTrac = str2double(handles.opcionesHormigonTrac.getSelection.getActionCommand);

inputs = handles.inputsMateriales;

% definir los parámetros del acero según el modelo utilizado
modelosAcero = {'Elastoplástico' 'Mander'};
parametros.modeloAcero = {opcionAcero modelosAcero{opcionAcero}};
if opcionAcero == 1
    parametros.ef = str2double(inputs(1).getText);
    parametros.Es2 = str2double(inputs(2).getText);
else%if opcionAcero == 2    
    parametros.esh = str2double(inputs(3).getText);
    parametros.esu = str2double(inputs(4).getText);
    parametros.ef = str2double(inputs(5).getText);
    parametros.Esh = str2double(inputs(6).getText);
    parametros.fsu = str2double(inputs(7).getText);
end
parametros.fy = handles.input_fy.UserData;

% definir los parámetros del hormigón
modelosHormigon = {'Saenz' 'Hognestad' 'Thorenfeldt 1' 'Thorenfeldt 2'};
parametros.modeloHormigon = {opcionHormigon modelosHormigon{opcionHormigon}};
parametros.e0 = str2double(inputs(8).getText);
if opcionHormigon ~= 1, parametros.eu = str2double(inputs(9).getText); end
parametros.fc = handles.input_fc.UserData;

parametros.modeloHormigonTrac = {1 'sin tracción'};
% modelosHormTrac = {'sin tracción' 'con tracción'};
% parametros.modeloHormigonTrac = {opcionHormigonTrac modelosHormTrac{opcionHormigonTrac}};

% almacenar los nuevos parámetros definidos
setappdata(handles.figure1,'parametrosMateriales',parametros);
handles.text48.String = ['Acero: ',parametros.modeloAcero{2}];
handles.text49.String = ['Hormigón: ',parametros.modeloHormigon{2},', ',parametros.modeloHormigonTrac{2}];

% recalcular el diagrama PM con los nuevos parámetros para las leyes constitutivas de los materiales
graficarDiagramas(handles)
end
