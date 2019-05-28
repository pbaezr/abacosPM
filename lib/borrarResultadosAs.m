function borrarResultadosAs(handles)
% Esta función borra los resultados asociados al cálculo de As

handles.output_rho.setText('')
handles.output_As.setText('')
delete(findobj(handles.axes1,'tag','punto'))
delete(findobj(handles.axes1,'tag','solucion'))
