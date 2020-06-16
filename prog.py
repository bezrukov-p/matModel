import plotly.express as px
import plotly.graph_objects as go
import numpy as np
#from scipy import integrate


a = 0
A = 0
b = 0
h = 0
o = 0
x0 = 0
y0 = 0
z0 = 0
points = []
pointsOnB = []
pointsOnB_Y = []
t = 0
eps = 0.000001
dt = 0.01
stepsMax = 10000
stop = False

def WhereZ(z): # показывает в какой области находится z. если выше плоскости b, то 1, если между b и -b, то 2, если ниже -b то 3
    if z > b:
        return 1
    elif (z > -b) and (z < b):
        return 2
    else:
        return 3

def ChangeSystem(whereZOld, upOrDown): # смена системы
    whereZNew = 0
    if upOrDown == "up":
        whereZNew = 1 if whereZOld == 2 else 2
    else:
        whereZNew = 3 if whereZOld == 2 else 2
    return whereZNew

def Division(t1, t2, whereZ, b): # метод деления пополам
    while np.abs(t2 - t1) > 2*eps:
        t = (t1 + t2)/2
        f_t1 = Function2(t1, whereZ).z - b
        f_t = Function2(t, whereZ).z - b
        f_t2 = Function2(t2, whereZ).z - b
        if f_t1*f_t < 0:
            t2 = t
        elif f_t*f_t2 < 0:
            t1 = t
    t = (t1 + t2)/2

    return t

def Division2(t1, t2, C1, C2, beta, b, aA): # метод деления пополам
    global t
    if t1 > t2:
        t1 = 0
    while np.abs(t2 - t1) > 2*eps:
        t = (t1 + t2)/2
        f_t1 = XY(t1, C1, C2, beta, b).y - aA
        f_t = XY(t, C1, C2, beta, b).y - aA
        f_t2 = XY(t2, C1, C2, beta, b).y - aA
        if f_t1*f_t < 0:
            t2 = t
        elif f_t*f_t2 < 0:
            t1 = t
    t = (t1 + t2)/2

    return t

def Function2(t, whereZ): # вычисляет следующую точку в зависимости от распложения предыдущей точки
    k = 0
    if whereZ == 1:
        k = -A
    elif whereZ == 2:
        k = 0
    elif whereZ == 3:
        k = A
    
    x = 1/(o*a)*(2*(np.exp(-h*t)*((h*h*k+0.5*h*x0-0.5*k+0.5*y0)*a+h*(x0-z0)/2)*np.sin(o*t)+(((h*k+0.5*x0)*a+0.5*x0-0.5*z0)*np.exp(-h*t)*np.cos(o*t)-(h-0.5*t)*k*a-0.5*x0+0.5*z0)*o))
    y = -1/(o*a)*(np.exp(-h*t)*((h*k+h*y0+x0)*a+x0-z0)*np.sin(o*t)+(np.exp(-h*t)*(k-y0)*np.cos(o*t)-k)*a*o)
    z = 1/(o*a)*(2*np.exp(-h*t)*((h*h*k+0.5*h*x0-0.5*k+0.5*y0)*a+h*(x0-z0)/2)*np.sin(o*t)+(2*((h*k+x0/2)*a+x0/2-z0/2)*np.exp(-h*t)*np.cos(o*t)+a*a*k*t+(-2*h*k+k*t-x0+z0)*a-x0+z0)*o)
    return Point(x,y,z,t)

def XY(t,C1, C2,beta, b):
    x = np.exp(-h*t)*(C1*np.cos(beta*t)+C2*np.sin(beta*t)) + b/(1+a)
    y = -h*np.exp(-h*t)*(C1*np.cos(beta*t)+C2*np.sin(beta*t)) + np.exp(-h*t)*(-C1*beta*np.sin(beta*t)+C2*beta*np.cos(beta*t))
    return Point(x, y, b, t)



def FunctionXY(upOrDown, b):
    global stop
    global t
    x0 = points[-1].x
    y0 = points[-1].y
    beta = np.sqrt((1+a)/a - h*h)
    C1 = x0 - b/(1+a)
    C2 = (y0+h*C1)/beta
    
    t = 0
    t += dt

    #######################
    pointTMP = XY(t,C1, C2, beta, b)
    points.append(pointTMP)
    points[-1].PrintPoint()
    t+=dt
    ##########################
    if (b == 0):
        while (points[-1].y + a*A) * (points[-1].y - a*A) < 0: # пока производные противоположных знаков
            pointTMP = XY(t,C1, C2, beta, b)
            points.append(pointTMP)
            points[-1].PrintPoint()
            t+=dt
            if t > 200:
                stop = True
                return points[-1]
    elif (upOrDown == "up"):
        while points[-1].y * (points[-1].y - a*A) < 0:  # пока производные противоположных знаков
            pointTMP = XY(t,C1, C2, beta, b)
            points.append(pointTMP)
            points[-1].PrintPoint()
            t+=dt
            if t > 200:
                stop = True
                return points[-1]
    elif upOrDown == "down":
        while points[-1].y * (points[-1].y + a*A) < 0:  # пока производные противоположных знаков
            pointTMP = XY(t,C1, C2, beta, b)
            points.append(pointTMP)
            points[-1].PrintPoint()
            t+=dt
            if t > 200:
                stop = True
                return points[-1]
    
    if (points[-1].y*points[-2].y < 0):#значит переход в нуле
        t = Division2(points[-2].t, points[-1].t, C1, C2, beta, b, 0)
        y0 = 0
    elif (points[-1].y - a*A)*(points[-2].y - a*A) < 0:# переход в a*A
        t = Division2(points[-2].t, points[-1].t, C1, C2, beta, b, a*A)
        y0 = a*A
    elif (points[-1].y + a*A)*(points[-2].y + a*A) < 0:# переход в -a*A
        t = Division2(points[-2].t, points[-1].t, C1, C2, beta, b, -a*A)
        y0 = -a*A
    
    pointTMP = XY(t,C1, C2,beta, b)
    pointTMP.y = y0
    pointsOnB_Y.append(pointTMP)
    points.append(pointTMP)
    points[-1].PrintPoint()

    return points[-1]

def PositionOnB(upOrDown,pointNew):
    global t
    global x0, y0, z0
    global stop
    _b = b
    if upOrDown == "down":
        _b = -b
    whereZOld = WhereZ(points[-1].z)        # запоминаем где находился z
    t = Division(points[-1].t, pointNew.t, whereZOld, _b) # вычисляем t при котором близко подойдёт к плоскости с помощью ЧМД/2
    pointOnB = Function2(t, whereZOld) # вычисляем точку, которая макс. близко к B
    pointOnB.z = _b
    points.append(pointOnB)
    pointOnB.PrintPoint()#
    pointsOnB.append(pointOnB)

    x0 = pointOnB.x
    y0 = pointOnB.y
    z0 = pointOnB.z
    t = 0
    if _b == 0:
        whereZNew = 1 if (whereZOld == 3) else 3
    else:
        whereZNew = ChangeSystem(whereZOld, upOrDown) # меняем чтобы вычислились новые значения в другой системе
    #C1, C2, C3 = Consts(pointOnB.x, pointOnB.y, pointOnB.z, whereZNew)  # новые константы при новых нач. условиях
    t += dt
    pointNew = Function2(t, whereZNew) #вычисляем след. точку в новой системе 0+dt


    if whereZOld == WhereZ(pointNew.z): # СКОЛЬЖЕНИЕ
        FunctionXY(upOrDown, _b)
        if stop:
            return pointNew
        x0 = points[-2].x # ставим начальные той точки, которая перешла границу, чтобы узнать в какое пространство пойдёт траектория
        y0 = points[-2].y
        z0 = points[-2].z
        t = dt
        whereZNew = WhereZ(Function2(t, whereZOld).z)  ### неважно какую подставляем, всё равно в одну сторону будет движение
        x0 = points[-1].x # ставим начальные той точки, которая на границе
        y0 = points[-1].y
        z0 = points[-1].z
        pointNew = Function2(t, whereZNew)
    
    points.append(pointNew)
    pointNew.PrintPoint()#
    return pointNew

class Point(object):
    x = 0
    y = 0
    z = 0
    t = 0

    def __init__(self, x, y, z, t):
        self.x = x
        self.y = y
        self.z = z
        self.t = t

    def PrintPoint(self):
        print("X = " + str(self.x) + " Y = " + str(self.y) + " Z = " + str(self.z) + " T = " + str(self.t))

def Start(_a,_A,_b,_h,_x0,_y0,_z0):  #начало рассчета при определенных параметрах
    global a,A,b,h,o,x0,y0,z0  
    global t
    global focus, stop
    a = _a
    A = _A
    b = _b
    h = _h
    o = np.sqrt(1 - h*h)
    x0 = _x0
    y0 = _y0
    z0 = _z0

    tailPointsOnB.clear() ######
    points.clear() ######
    pointsOnB.clear()#####
    pointsOnB_Y.clear()#######
    t = 0

    points.append(Point(x0, y0, z0, t))
    points[-1].PrintPoint()

    while True:
        t += dt
        if t > 150:
            stop = True
            focus = True
            print(a)
        if len(pointsOnB) > 150  or stop:
            break

        pointNew = Function2(t, WhereZ(points[-1].z))
        if (points[-1].z - b)*(pointNew.z - b) < 0: # значит перешагнул плоскость b  
            pointNew = PositionOnB("up",pointNew)
        elif (points[-1].z + b)*(pointNew.z + b) < 0: # перешагнул плоскость -b
            pointNew = PositionOnB("down",pointNew)
        else:
            points.append(pointNew)
            pointNew.PrintPoint()#


    tail = (int)(len(pointsOnB)/5)
    for i in range(tail):
        tailPointsOnB.append(pointsOnB[-i - 1])
    tailPointsOnB_All.append(tailPointsOnB.copy())
    return

focus = False

a1 = 0.1
a2 = 1.5   ######## a = 1.986   A = 1
a_array = []
A1 = 0.1
A2 = 6
A_array = []
count = 10
step_a = (a2-a1)/count
step_A = (A2-A1)/count
for i in range(count + 1):
    a_array.append(a1 + i*step_a)
    A_array.append(A1 + i*step_A)

tailPointsOnB_All = []
tailPointsOnB = []
if True:
    for i in range(count + 1):
        if focus:
            break
        Start(a_array[i], 1, 0.5, 0.1, 3, 2, 5)
else:
    Start(3, 1, 0.5, 0.1, 3, 2, 5)

x_array_array_out = []
y_array_array_out = []
a_array_array_out = []
A_array_array_out = []
for i in range(len(tailPointsOnB_All)):
    x_tmp = []
    y_tmp = []
    a_tmp = []
    A_tmp = []
    for j in range(len(tailPointsOnB_All[i])):
        x_tmp.append(tailPointsOnB_All[i][j].x)
        y_tmp.append(tailPointsOnB_All[i][j].y)
        a_tmp.append(a_array[i])
        A_tmp.append(A_array[i])
    x_array_array_out.append(x_tmp.copy())
    a_array_array_out.append(a_tmp.copy())
    y_array_array_out.append(y_tmp.copy())
    A_array_array_out.append(A_tmp.copy())



fig = go.Figure()
for i in range(len(tailPointsOnB_All)):
    fig.add_trace(go.Scatter(x=a_array_array_out[i], y=x_array_array_out[i],
                        mode='markers',
                        name='markers'))
fig.write_html('x(a).html')

fig = go.Figure()
for i in range(len(tailPointsOnB_All)):
    fig.add_trace(go.Scatter(x=a_array_array_out[i], y=y_array_array_out[i],
                        mode='markers',
                        name='markers'))
fig.write_html('y(a).html')





x_out = []
y_out = []
z_out = []
t_out = []
t_out_tmp = -dt
b_out = []
b_out_minus = []

for i in range(len(points)):
    #if i % 10 == 0:
    t_out_tmp += dt
    x_out.append(points[i].x)
    y_out.append(points[i].y)
    z_out.append(points[i].z)
    t_out.append(t_out_tmp)
    b_out.append(b)
    b_out_minus.append(-b)

fig = go.Figure()
fig.add_trace(go.Scatter(x=t_out, 
y=x_out,
mode='lines',
name='lines'))
fig.write_html('2dLine_xt_1.html')

fig = go.Figure()
fig.add_trace(go.Scatter(x=t_out, 
y=y_out,
mode='lines',
name='lines'))
fig.write_html('2dLine_yt_1.html')

fig = go.Figure()
fig.add_trace(go.Scatter(x=t_out, 
y=z_out,
mode='lines',
name='lines'))
fig.add_trace(go.Scatter(x=t_out, 
y=b_out,
mode='lines',
name='lines'))
fig.add_trace(go.Scatter(x=t_out, 
y=b_out_minus,
mode='lines',
name='lines'))
fig.write_html('2dLine_zt_1.html')

fig = go.Figure(go.Scatter3d(
x = x_out,
y = y_out,
z = z_out
))

fig.write_html('3dLine1.html')

