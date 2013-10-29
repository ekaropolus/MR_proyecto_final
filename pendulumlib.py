#(i) y (ii) y (iv)
from math import *
from numpy import *
from time import *
from scipy.linalg import *
from pylab import *
NAME = 'dummy'
FNAME = 'fdummy'
SY = ''
##str(localtime()[3])+str(localtime()[4])+str(localtime()[5])

def newton_raphson(x,f,df,toleration,iterations):
    k = 1    
    while k<=iterations:
        xi=x.copy()
        a = df(0,xi)
        b = -f(0,xi)
        y = solve(a,b)
        x = x + y
        if sqrt(vdot(y,y)) <= toleration:
            k = iterations
            return x
            
        k = k +1
def kutta(t,x,h,f):
    global NAME
    NAME = 'kutta'
    x0=x.copy()
    k1 =h*f(t,x0)
    x0=x.copy()
    k2 = h*f(t,x0+.5*k1)
    x0=x.copy()
    k3 = h*f(t+.5*h,x0+.5*k2)
    x0=x.copy()
    k4 = h*f(t+h,x0+k3)
    x  = x + (k1 + 2*k2 + 2*k3 + k4)/6
    return x

def backkutta(t,x,h,f):
    global NAME
    NAME = 'kutta'
    x0=x.copy()
    k1 =h*f(t,x0)
    x0=x.copy()
    k2 = h*f(t,x0-.5*k1)
    x0=x.copy()
    k3 = h*f(t-.5*h,x0-.5*k2)
    x0=x.copy()
    k4 = h*f(t-h,x0-k3)
    x  = x - (k1 + 2*k2 + 2*k3 + k4)/6
    return x

def integrar(t0,tf,h,x0,f,imeth,lx,ly):
    global NAME
    global FNAME
    global SY
    tp_list = arange(t0,tf,h)
    tp_lx = arange(t0,tf,h)
    tp_ly = arange(t0,tf,h)
    inx = 0
    for x in tp_list:
        x0 = imeth(x,x0,h,f)
        lx.append(x0[0])
        ly.append(x0[1])
##    t = t0
##    x0 = imeth(t,x0,h,f)
##    salida = 'c:\e'+ NAME + FNAME + SY +'.txt'
##    outto = open(salida,"a")
##    outto.write(str(x0[0])+' '+str(x0[1])+'\n')
##    while t<tf:
##        t = t + h
##        x0 = imeth(t,x0,h,f)
##        outto.write(str(x0[0])+' '+str(x0[1])+'\n')
##    outto.close()
##    tf = t - h
def dintegrar(t0,tf,h,x0,f,imeth):
    global NAME
    global FNAME
    global SY
    t = t0
    x0 = imeth(t,x0,h,f)
    salida = 'c:\e'+ NAME + FNAME + SY +'.txt'
    outto = open(salida,"a")           
    outto.write(str(cos(x0[0]))+' '+str(sin(x0[0]))+' '+str(x0[1])+'\n')
    while t<tf:
        t = t + h
        x0 = imeth(t,x0,h,f)
        outto.write(str(cos(x0[0]))+' '+str(sin(x0[0]))+' '+str(x0[1])+'\n')
    outto.close()
    tf = t - h
def dplano_fase(t0,tf,h,l,f):
    for x in l:
        x0 = array([.0,.0])
        x0[0] = x[0]
        x0[1] = x[1]
        dintegrar(t0,tf,h,x0,f,kutta)
def plano_fase(t0,tf,h,l,f):
    lx = list()
    ly = list()
    for x in l:
        x0 = array([.0,.0])
        x0[0] = x[0]
        x0[1] = x[1]
        integrar(t0,tf,h,x0,f,backkutta,lx,ly)
        integrar(t0,tf,h,x0,f,kutta,lx,ly)
    a = array(lx)
    b = array(ly)
    plot(a,b,'.',markersize=.1,)
	
def creared(red,size,x0,x,y):
    dummy = x0.copy()
    xi = dummy[0]
    x = xi + x
    y = dummy[1] + y
    while xi<=x:
        dummy = x0.copy()
        yi = dummy[1]
        while yi<=y:
            elemento = (xi,yi)
            red.append(elemento)
            if -xi <> xi:
                elemento = (-xi,yi)
                red.append(elemento)
            if -yi<>yi:
                elemento = (xi,-yi)
                red.append(elemento)
            if -xi<>xi and -yi <> yi:
                elemento = (-xi,-yi)
                red.append(elemento)
            yi = yi + size
        xi = xi + size
def crealinea (linea,t0,tf,h,x,y):
    while t0 <= tf:
        z = matrix([x])
        lt = z + t0*y.copy()
        elemento = (lt[0,0],lt[0,1])
        linea.append(elemento)
        t0 = t0 + h
def fixed_points (l,f,df,nr):
    z = array([0.,0.])
    s = list()
    for x in l:
        z[0]=x[0]
        z[1]=x[1]
        zi = nr(z,f,df,.0001,1000)
        if zi <> 'None':
            elemento = (zi[0],zi[1])
            if elemento not in s:
               s.append(elemento)        
    return s
def fixed_points_clasification(l,f,df,messages):
    z0 = array([0.,0.])
    s = list()
    for x in l:
        z0[0] = x[0]
        z0[1] = x[1]
        z = z0.copy()
        w = df(0.,z)
        b = eig(w)
        alpha1 =  b[0][0]
        alpha2 =  b[0][1]
        delta = alpha1*alpha2
        tau = real(alpha1) + real(alpha2)
        nau = pow(tau,2)-4*delta
        if delta < 0:
            elemento = (z,b,'s')
        if delta > 0:
            if tau == 0.:
                elemento = (z,b,'c')
            else:
                if nau < 0.:
                    elemento = (z,b,'e')
                if nau > 0.:
                    elemento = (z,b,'n')
                if nau == 0.:
                    elemento = (z,b,'d')
        s.append(elemento)
    for x in s:
        tipo = x[2]
        ml = x[0].copy()
        if tipo == 's':
            component = ((ml[0],ml[1]),'Saddles Point')
        if tipo == 'c':
            component = ((ml[0],ml[1]),'Center Point')
        if tipo == 'e':
            component = ((ml[0],ml[1]),'Spiral Point')
        if tipo == 'n':
            component = ((ml[0],ml[1]),'Node Point')
        if tipo == 'd':
            component = ((ml[0],ml[1]),'Degenerate Point')
        messages.append(component)
    return s
def stable_manifolds(fpcl,messages):
    #print 'paso por aca'
    l = list()
    for x in fpcl:        
        option = x[2]
        vi = matrix([[1.,0.]])
        vi = vi.T
        vj = matrix([[0.,1.]])
        vj = vj.T
        #print option
        if option == 's':
            a = x[0].copy()
            evalue1 = x[1][0][0]
            evector1 =x[1][1].copy()
            evalue2 = x[1][0][1]
            evector2 =x[1][1].copy()
            #print evector2,'vi',vi,'vj',vj
            evector2 = evector2*vj
            evector1 = evector1*vi
            evector1 = evector1.T
            evector2 = evector2.T
            #print a
            if  evalue1 < 0:
                crealinea(l,-.1,.1,.1,a,evector1)
                component = ((a[0],a[1]),'Variedad estable para tipo Saddles',(evector1[0,0],evector1[0,1]))
                messages.append(component)
            else:
                crealinea(l,-.1,.1,.1,a,evector2)
                component = ((a[0],a[1]),'Variedad estable para tipo Saddles',(evector2[0,0],evector2[0,1]))
                messages.append(component)
        if option == 'n':
            a = x[0].copy()
            evalue1 = x[1][0][0]
            evector1 =x[1][1].copy()
            evalue2 = x[1][0][1]
            evector2 =x[1][1].copy()
            evector2 = evector2*vj
            evector1 = evector1*vi
            evector1 = evector1.T
            evector2 = evector2.T
            tau = real(evalue1) + real(evalue2)
            #print a
            if tau < 0:
                crealinea(l,-.1,.1,.1,a,evector1)
                crealinea(l,-.1,.1,.1,a,evector2)
                component = ((a[0],a[1]),'Variedad estable para tipo Node',(evector1[0,0],evector1[0,1]))
                messages.append(component)
                component = ((a[0],a[1]),'Variedad estable para tipo Node',(evector2[0,0],evector2[0,1]))
                messages.append(component)
        if option == 'e':
            a = x[0].copy()
            evalue1 = x[1][0][0]
            evector1 =x[1][1].copy()
            evector1 = evector1*vi
            evector1 = evector1.T
            evalue2 = x[1][0][1]
            tau = real(evalue1)*real(evalue2)
            #print a
            if tau < 0.:
                crealinea(l,0,.1,.1,a,evector1)
                component = ((a[0],a[1]),'Variedad estable para tipo Spiral',(evector1[0,0],evector1[0,1]))
                messages.append(component)
        if option == 'c':
            a = x[0].copy()
            evector1 =x[1][1].copy()
            evector1 = evector1*vi
            evector1 = evector1.T
            crealinea(l,0,.5,.5,a,evector1)
            component = ((a[0],a[1]),'Variedad estable para tipo central',(evector1[0,0],evector1[0,1]))
            messages.append(component)
    return l
def unstable_manifolds(fpcl,messages):
    l = list()
    for x in fpcl:        
        option = x[2]
        vi = matrix([[1.,0.]])
        vi = vi.T
        vj = matrix([[0.,1.]])
        vj = vj.T
        if option == 's':
            a = x[0].copy()
            evalue1 = x[1][0][0]
            evector1 =x[1][1].copy()
            evalue2 = x[1][0][1]
            evector2 =x[1][1].copy()
            evector2 = evector2*vi
            evector1 = evector1*vi
            evector1 = evector1.T
            evector2 = evector2.T
            #print a
            if evalue1 > 0:
                crealinea(l,-.1,.1,.1,a,evector1)
                component = ((a[0],a[1]),'Variedad no estable para tipo Saddles',(evector1[0,0],evector1[0,1]))
                messages.append(component)
            else:
                crealinea(l,-.1,.1,.1,a,evector2)
                component = ((a[0],a[1]),'Variedad no estable para tipo Saddles',(evector2[0,0],evector2[0,1]))
                messages.append(component)
        if option == 'n':
            a = x[0].copy()
            evalue1 = x[1][0][0]
            evector1 =x[1][1][0].copy()
            evalue2 = x[1][0][1]
            evector2 =x[1][1][1].copy()
            evector2 = evector2*vj
            evector1 = evector1*vi
            evector1 = evector1.T
            evector2 = evector2.T
            tau = real(evalue1) + real(evalue2)
            if tau > 0:
                crealinea(l,-.1,.1,.1,a,evector1)
                crealinea(l,-.1,.1,.1,a,evector2)
                component = ((a[0],a[1]),'Variedad no estable para tipo Node',(evector1[0,0],evector1[0,1]))
                messages.append(component)
                component = ((a[0],a[1]),'Variedad no estable para tipo Node',(evector2[0,0],evector2[0,1]))
                messages.append(component)
        if option == 'e':
            a = x[0].copy()
            evalue1 = x[1][0][0]
            evector1 =x[1][1][0].copy()
            evalue2 = x[1][0][1]
            evector1 = evector1*vi
            evector1 = evector1.T
            tau = real(evalue1)*real(evalue2)
            if tau > 0:
                crealinea(l,0,.5,.5,a,evector1)
                component = ((a[0],a[1]),'Variedad no estable para tipo Spiral',(evector1[0,0],evector1[0,1]))
                messages.append(component)
    return l
