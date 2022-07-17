#Throughtout this code we are converting all the quantities to respective SI units
#Convereting all the quantities to the standard units eases out calculations
from matplotlib.axis import Axis
import cmath
import warnings
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
from matplotlib.widgets import Slider, Button

# Create subplot
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.35)

TR=10#float(input("Enter the repitition pulse"))
PD=0.01#float(input("Enter the time duration of the pulse "))
Theta=90#float(input("Enter the flip angle in degrees"))
Lx=1#float(input("Enter the field of view along x-axis in centimetre"))
Ly=1#float(input("Enter the field of view along y-axis in centimetre"))
Lz=1#float(input("Enter the field of view along z-axis in centimetre"))
Dx=2#int(input("Enter the number of discretizations along x-axis"))
Dy=2#int(input("Enter the number of discretizations along y-axis"))
Dz=1#int(input("Enter the number of discretizations along z-axis"))

B0=0.3E-6 #Static magnetic field
M0=1 #initial magnetization
fs=5E1 #Sampling frequency
f=1/TR #frequency of the pulse
simdur=20#simulation time
y_g=42.58*1E6 #gyromagnetic ratio
T1=3 #T1 relaxation
T2=3 #T2 relaxation

#Here we are taking the field of view Lx,Ly,Lz in centimetre and to convert it to metres we divide by 100
if(Dx!=1):
    x=np.linspace(-Lx/200, Lx/200, Dx, endpoint=True)
else:
    x = []
if(Dy!=1):
    y=np.linspace(-Ly/200, Ly/200, Dy, endpoint=True)
else:
    y = []
if(Dz!=1):
    z=np.linspace(-Lz/200, Lz/200, Dz, endpoint=True)
else:
    z = []
t=np.arange(0,simdur,1/fs)  #time vector
B1=Theta/(360*y_g*PD)#calculation of excitation magnetic field amplitude

#pre-allocation of x,y,z of the magnetization vector
D1=len(x)
D2=len(y)
D3=len(z)
cols=len(t)

Mx = np.zeros((D1,D2,D3,cols))
My = np.zeros((D1,D2,D3,cols))
Mz = np.ones((D1,D2,D3,cols))
Comp_Sig =np.zeros((D1,D2,D3,cols))
B_shim=np.zeros((D1,D2,D3))
iota=cmath.sqrt(-1)

#generation of B1 excitation pulse
base_signal1=((signal.square(2 * np.pi * f * (t+1*TR), duty=PD/TR)+1)*0.5)
B_1=B1*base_signal1

#Functions to calculate magnetization
def func1(x1,y1,z1,t1,B):
    Mx1=(x1*np.cos(2*np.pi*y_g*B*t1)-y1*np.sin(2*np.pi*y_g*B*t1))*np.e**(-t1/T2)
    My1=(x1*np.sin(2*np.pi*y_g*B*t1)+y1*np.cos(2*np.pi*y_g*B*t1))*np.e**(-t1/T2)
    Mz1=((-z1+M0/T1)*t1) + z1
    M1=(Mx1,My1,Mz1)
    return M1;
def func2(x2,y2,z2,t2):
    Mx2=x2*np.cos(2*np.pi*y_g*B1*t2)+z2*np.sin(2*np.pi*y_g*B1*t2)
    My2=y2
    Mz2=-x2*np.sin(2*np.pi*y_g*B1*t2)+z2*np.cos(2*np.pi*y_g*B1*t2)
    M2=(Mx2,My2,Mz2)
    return M2;
#Function for calculation of Mx,My,Mz
def compute(c0,c1,c2,c3,c4,c5,c6,c7,c8):
    for i in range(cols - 1):
        for j in range(D1 - 1):
            for k in range(D2 - 1):
                for l in range(D3 - 1):
                    #if(Dx==1):
                       # B_shim[j][k][l] = c0 + (c1 * x) + (c2 * y[k]) + (c3 * z[l]) + (c4 * (3 * ((x * x) - (y[k] * y[k])))) + (c5 * ((z[l] * z[l]) - 0.5 * ((x * x) + (y[k] * y[k])))) + (c6 * 6 * 0 * y[k]) + (c7 * 3 * y[k] * z[l]) + (c8 * 3 * z[l] * x)
                    #if(Dy==1):
                      #  B_shim[j][k][l] = c0 + (c1 * x[j]) + (c2 * y) + (c3 * z[l]) + (c4 * (3 * ((x[j] * x[j]) - (y * y)))) + (c5 * ((z[l] * z[l]) - 0.5 * ((x[j] * x[j]) + (y * y)))) + (c6 * 6 * x[j] * y) + (c7 * 3 * y * z[l]) + (c8 * 3 * z[l] * x[j])
                    #if(Dz==1):
                       # B_shim[j][k][l] = c0 + (c1 * x[j]) + (c2 * y[k]) + (c3 * z) + (c4 * (3 * ((x[j] * x[j]) - (y[k] * y[k])))) + (c5 * ((z * z) - 0.5 * ((x[j] * x[j]) + (y[k] * y[k])))) + (c6 * 6 * x[j] * y[k]) + (c7 * 3 * y[k] * z) + (c8 * 3 * z * x[j])
                    #else:
                    B_shim[j][k][l] = c0 + (c1 * x[j]) + (c2 * y[k]) + (c3 * z[l]) + (c4 * (3 * ((x[j] * x[j]) - (y[k] * y[k])))) + (c5 * ((z[l] * z[l]) - 0.5 * ((x[j] * x[j]) + (y[k] * y[k])))) + (c6 * 6 * x[j] * y[k]) + (c7 * 3 * y[k] * z[l]) + (c8 * 3 * z[l] * x[j])
                    # Shim fields where we are considering zero order,first order and second order
                    B = B0 + B_shim[j][k][l]
                    dt = t[i + 1] - t[i]
                    if (B_1[i] == 0):
                        args = func1(Mx[j][k][l][i], My[j][k][l][i], Mz[j][k][l][i], dt, abs(B))
                        Mx[j][k][l][i + 1] = args[0]
                        My[j][k][l][i + 1] = args[1]
                    else:
                        args = func2(Mx[j][k][l][i], My[j][k][l][i], Mz[j][k][l][i], dt)
                        Mx[j][k][l][i + 1] = args[0]
                        My[j][k][l][i + 1] = args[1]
    for i in range(len(t) - 1):
        for j in range(D1 - 1):
            for k in range(D2 - 1):
                for l in range(D3 - 1):
                    Mx[0][0][0][i] = Mx[0][0][0][i] + Mx[j][k][l][i]
                    My[0][0][0][i] = My[0][0][0][i] + My[j][k][l][i]
    warnings.simplefilter("ignore", np.ComplexWarning)
    Comp_Sig[0][0][0] = Mx[0][0][0] + (iota * My[0][0][0])
    return Comp_Sig[0][0][0];

#Formation of graph with all Shim coefficients as zeros
fft_sig=np.abs(fftpack.fft(compute(0,0,0,0,0,0,0,0,0)))
freqx=fftpack.fftfreq(len(compute(0,0,0,0,0,0,0,0,0)))*fs
p,=plt.plot(freqx,fft_sig)

#Sliders Creation
axH0=plt.axes([0.15, 0.29, 0.65, 0.03])
axHx=plt.axes([0.15, 0.26, 0.65, 0.03])
axHy=plt.axes([0.15, 0.23, 0.65, 0.03])
axHz=plt.axes([0.15, 0.2, 0.65, 0.03])
axHx2y2=plt.axes([0.15, 0.17, 0.65, 0.03])
axHz2=plt.axes([0.15, 0.14, 0.65, 0.03])
axHxy=plt.axes([0.15,0.11, 0.65, 0.03])
axHyz=plt.axes([0.15,0.08 , 0.65, 0.03])
axHzx=plt.axes([0.15, 0.05, 0.65, 0.03])

H0val= Slider(axH0, 'H0',-10*1E-6, 10*1E-6, 0*1E-6)
Hxval= Slider(axHx, 'Hx',-10*1E-6, 10*1E-6, 0*1E-6)
Hyval= Slider(axHy, 'Hy',-10*1E-6, 10*1E-6, 0*1E-6)
Hzval= Slider(axHz, 'Hz',-10*1E-6, 10*1E-6, 0*1E-6)
#-30*1E-6, 30*1E-6, 0*1E-6
Hx2y2val= Slider(axHx2y2, 'Hx2y2',-18*1E-4,18*1E-4, 0*1E-4)
Hz2val= Slider(axHz2, 'Hz2',-18*1E-4,18*1E-4, 0*1E-4)
#-20*1E-6,20*1E-6,0*1E-6
Hxyval= Slider(axHxy, 'Hxy',-18*1E-4,18*1E-4, 0*1E-4)
Hyzval= Slider(axHyz, 'Hyz',-18*1E-4, 18*1E-4, 0*1E-4)
Hzxval= Slider(axHzx, 'Hzx',-18*1E-4,18*1E-4, 0*1E-4)
def update(val):
    H0=H0val.val
    Hx=Hxval.val
    Hy=Hyval.val
    Hz=Hzval.val
    Hx2y2=Hx2y2val.val
    Hz2=Hz2val.val
    Hxy=Hxyval.val
    Hyz=Hyzval.val
    Hzx=Hzxval.val

    fft_sig = np.abs(fftpack.fft(compute(H0,Hx,Hy,Hz,Hx2y2,Hz2,Hxy,Hyz,Hzx)))
    freqx = fftpack.fftfreq(len(compute(H0,Hx,Hy,Hz,Hx2y2,Hz2,Hxy,Hyz,Hzx))) * fs

    p.set_ydata(fft_sig)
    fig.canvas.draw_idle()
    #p.xlim([Fc - 6, Fc + 6])



#Run-time error handling
import warnings
warnings.simplefilter("ignore", np.ComplexWarning)

for r in range(len(fft_sig)):
    if (fft_sig[r] == max(fft_sig)):
        fc = r
        break
Fc = (fc / len(fft_sig)) * fs
print(Fc)
H0val.on_changed(update)
Hxval.on_changed(update)
Hyval.on_changed(update)
Hzval.on_changed(update)
Hx2y2val.on_changed(update)
Hz2val.on_changed(update)
Hxyval.on_changed(update)
Hyzval.on_changed(update)
Hzxval.on_changed(update)

# Create axes for reset button and create button
resetax=plt.axes([0.01, 0.9, 0.1, 0.04])
setax=plt.axes([0.12,0.9,0.1,0.04])
button=Button(resetax, 'Reset', color='gold',hovercolor='skyblue')
button1=Button(setax,'Set',color='gold',hovercolor='skyblue')

# Create a function resetSlider to set slider to 0
# initial values when Reset button is clicked
def setSlider(event):
    H0val.set_val(H0val.val)
    Hxval.set_val(Hxval.val)
    Hyval.set_val(Hyval.val)
    Hzval.set_val(Hzval.val)
    Hx2y2val.set_val(Hx2y2val.val)
    Hz2val.set_val(Hz2val.val)
    Hxyval.set_val(Hxyval.val)
    Hyzval.set_val(Hyzval.val)
    Hzxval.set_val(Hzxval.val)
button1.on_clicked(setSlider)

def resetSlider(event):
    H0val.reset()
    Hxval.reset()
    Hyval.reset()
    Hzval.reset()
    Hx2y2val.reset()
    Hz2val.reset()
    Hxyval.reset()
    Hyzval.reset()
    Hzxval.reset()
# Call resetSlider function when clicked on reset button
button.on_clicked(resetSlider)
plt.show()