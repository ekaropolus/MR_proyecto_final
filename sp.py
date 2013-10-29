#proyecto final
from pendulumlib import *
#Varialbe para el pendulo amortiguado
DPB = array([.5])
#Variable para el pendulto con torca
TPA = array([.5])
#Variable para el pendulo no-lineal
NPA = array([.25])

def chaos_e(t,x):
    x0 = x[0]
    x1 = x[1]
    x[0] = x1
    x[1] = .22*x1-sin(x0)+2.7*cos(t)
    return x
#Jacobiano para el pendulo amortiguado    
def chaos_de(t,x):
    x0 = x[0]
    x1 = x[1]
    a = matrix ([[0,1.],[-cos(x0),.22]])
    return a
#Sistema para el pendulo simple
def g(t,x):
    x0 = x[0]
    x1 = x[1]
    x[0] = x1
    x[1] = -sin(x0)
    return x
#Jacobiano para el pendulo simple
def dg(t,x):
    x0 = x[0]
    x1 = x[1]
    a = matrix ([[0,1.],[-cos(x0),0.]])
    return a
#Sistema para el pendulo amortiguado
def e(t,x):
    x0 = x[0]
    x1 = x[1]
    x[0] = x1
    x[1] = -DPB*x1-sin(x0)
    return x
#Jacobiano para el pendulo amortiguado    
def de(t,x):
    x0 = x[0]
    x1 = x[1]
    a = matrix ([[0,1.],[-cos(x0),-DPB]])
    return a
#Sistema para el pendulo con torca    
def f(t,x):
    x0 = x[0]
    x1 = x[1]
    x[0] = x1
    x[1] =TPA-sin(x0)
    return x
#Jacobiano para el pendulo con torca    
def df(t,x):
    x0 = x[0]
    x1 = x[1]
    a = matrix ([[0,1.],[-cos(x0),0]])
    return a
#Sistema para el pendulo no lineal
def h(t,x):
    x0 = x[0]
    x1 = x[1]
    x[0] = x1
    x[1] =-x1*(1+NPA*cos(x0))-sin(x0)
    return x
#Jacobiano para el pendulo no lineal    
def dh(t,x):
    x0 = x[0]
    x1 = x[1]
    a = matrix ([[0,1.],[NPA*x1*sin(x0)-cos(x0),-1-NPA*cos(x0)]])
    return a


#Sistema para el pendulo no lineal
def i(t,x):
    x0 = x[0]
    x1 = x[1]
    x[0] = x1
    x[1] =-pow(x0,3)
    return x
#Jacobiano para el pendulo no lineal    
def di(t,x):
    x0 = x[0]
    x1 = x[1]
    a = matrix ([[0,1],[0,0]])
    return a

def graph_u_pendulum(net_size,net_center,max_x_net,max_y_net,function,derivate,t0,tf,h):
	l = list()	
	creared (l,net_size,net_center,max_x_net,max_y_net)
	fpl = list()
	fpl = fixed_points (l,function,derivate,newton_raphson)
	fpcl = list()
	messages = list()
	fpcl = fixed_points_clasification(fpl,g,dg,messages)
	usml = list()
	messages2 = list()
	usml = unstable_manifolds(fpcl,messages2)
	
	    
	print 'comienza plano fase'
	plano_fase(t0,tf,h,usml,function)
	print 'termina plano fase'
	return fpl, fpcl, messages, usml

def graph_s_pendulum(net_size,net_center,max_x_net,max_y_net,function,derivate,t0,tf,h):
	l = list()
	creared (l,net_size,net_center,max_x_net,max_y_net)
	fpl = list()
	fpl = fixed_points (l,function,derivate,newton_raphson)
	fpcl = list()
	messages = list()
	fpcl = fixed_points_clasification(fpl,function,derivate,messages)
	sml = list()
	messages2 = list()
	sml = stable_manifolds(fpcl,messages2)
	
	    
	print 'comienza plano fase'
	plano_fase(t0,tf,h,sml,function)
	print 'termina plano fase'
	return fpl, fpcl, messages, sml


def fixp_pendulum(net_size,net_center,max_x_net,max_y_net,function,derivate):
	l = list()
	creared (l,net_size,net_center,max_x_net,max_y_net)
	fpl = list()
	fpl = fixed_points (l,function,derivate,newton_raphson)
	fpcl = list()
	messages = list()

	fpcl = fixed_points_clasification(fpl,function,derivate,messages)
	i = 1
	
	return messages
	#print 'Se imprimen los puntos fijos'
	
	#for x in messages:
	#	print i, x
	#	i = i +1
	
def graph_pendulum(net_size,net_center,max_x_net,max_y_net,function,derivate,t0,tf,h):
   l = list()
   creared (l,net_size,net_center,max_x_net,max_y_net)
   print 'comienza plano fase'
   plano_fase(t0,tf,h,l,function)
   print 'finaliza plano fase'
