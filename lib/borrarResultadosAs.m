function borrarResultadosAs(handles)
% Esta funci�n borra los resultados asociados al c�lculo de As

handles.output_rho.setText('')
handles.output_As.setText('')
delete(findobj(handles.axes1,'tag','punto'))
delete(findobj(handles.axes1,'tag','solucion'))
