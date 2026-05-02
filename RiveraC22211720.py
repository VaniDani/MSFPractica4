"""
Práctica 4: Sistema endocrino

Departamento de Ingeniería Eléctrica y Electrónica, Ingeniería Biomédica
Tecnológico Nacional de México [TecNM - Tijuana]
Blvd. Alberto Limón Padilla s/n, C.P. 22454, Tijuana, B.C., México

Nombre del alumno: Vania Daniela Rivera Durán
Número de control: C22211720
Correo institucional: l22211720@tectijuana.edu.mx

Asignatura: Modelado de Sistemas Fisiológicos
Docente: Dr. Paul Antonio Valle Trujillo; paul.valle@tectijuana.edu.mx
"""


# Librerías para cálculo numérico y generación de gráficas
import numpy as np
import matplotlib.pyplot as plt
import control as ctrl


x0,t0,tend,dt,w,h= 0,0,15,1e-3,10,5
N = round((tend-t0)/dt) + 1
t = np.linspace(t0,tend,N)
u = np.ones(N) #Step

def endo(R1,R2,C,L):
    num=[R2]
    den=[R1*C*L,L+R1*R2*C,R1+R2]
    sys=ctrl.tf(num,den)
    return sys

#%%
#Control
R1,R2,L,C= 1e3,100e3,100e-3,1e-6
num = [R2]
den = [R1*C*L,L+R1*R2*C,R1+R2]
sys = ctrl.tf(num,den)
print(f"Función de transferencia del sistema: {sys}\n")
print(f"Lambda 1 Control: {np.roots(den)[0]}")
print(f"Lambda 2 Control: {np.roots(den)[1]}\n")

#Caso
R1,R2,L,C= 1e3,1e3,100e-3,1000e-6
num = [R2]
den = [R1*C*L,L+R1*R2*C,R1+R2]
sys = ctrl.tf(num,den)
print(f"Función de transferencia del sistema: {sys}\n")
print(f"Lambda 1 Caso: {np.roots(den)[0]}")
print(f"Lambda 2 Caso: {np.roots(den)[1]}\n")
#%%
#Función de transferencia: Normotenso
R1,R2,C,L= 1e3,100e3,1e-6,100e-3
syscontrol= endo(R1,R2,C,L)
print(f'Función de transferencia de normotenso (control): {syscontrol}')

#Función de transferencia: Hipotenso
R1,R2,C,L= 1e3,1e3,1000e-6,100e-3
syscaso= endo(R1,R2,C,L)
print(f'Función de transferencia de hipotenso (caso 1): {syscaso}')
#%%

#Respuestas en lazo abierto
_,Pp0=ctrl.forced_response(syscontrol,t,u,x0)
_,Pp1=ctrl.forced_response(syscaso,t,u,x0)

#Controlador PID
def controlador(KP,KI,KD,sys):
    Cr= 1e-6
    Re=1/(KI*Cr)
    Rr=KP*Re
    Ce=KD/Rr
    numPID = [Re*Rr*Cr*Ce,(Re*Ce+Rr*Cr),1]
    denPID = [Re*Cr,0]
    PID=ctrl.tf(numPID,denPID)
    X=ctrl.series(PID,sys)
    sysPID = ctrl.feedback(X,1,sign=-1)
    print(f"El velor de Capacitancia Cr es de {Cr} Faradios \n")
    print(f"El valor de resistencia de Re es de {Re} Ohms \n")
    print(f"El valor de resistencia de Rr es de {Rr} Ohms \n")
    print(f"El velor de Capacitancia Ce es de {Ce} Faradios \n")
    return sysPID    
casoPID= controlador(1315.6501,314901.7711,0.31035,syscaso)
print(f'Función de transferencia de hipotenso en lazo cerrado: {casoPID}')


#%%
#Respuestas del sistema de control en lazo cerrado
_,PID1=ctrl.forced_response(casoPID,t,Pp0,x0)


fg1=plt.figure()
plt.plot(t,Pp0,'-',linewidth=1,color=[0.54, 0.52, 0.20],label='Vs(t): Control')
plt.plot(t,Pp1,'--',linewidth=1,color=[0.66, 0.16, 0.11],label='Vs(t): Caso')
plt.plot(t,PID1,':',linewidth=1, color=[0.46, 0.10, 0.10],label='PID(t): Caso')
plt.grid(False)
plt.xlim(0,10); plt.xticks(np.arange(0,11,1))
plt.ylim(0,1.1); plt.yticks(np.arange(0,1.1,0.2))
plt.xlabel('t[s]')
plt.ylabel('Vs(t) [V]')
plt.legend(bbox_to_anchor=(0.5,-0.2),loc='center',ncol=3)
plt.title('Control vs Caso', fontsize=10)
plt.show()
fg1.set_size_inches(w,h)
fg1.tight_layout()
fg1.savefig('Endocrino lazo cerrado python.pdf')
