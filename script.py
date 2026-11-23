import numpy as np 

rho= 1000
Diameter= 0.15
mu= 0.001
Delta_P= 100000
K=0.4
epsilon = 0.00015

def flow_solver(u0,P,D, percent_margin=1):
    for i in range(10) :
        f =  4/(-3.6*np.log10( (6.9/(rho * u0 * D)/mu) + (epsilon/(3.7*D))**(10/9) ))**2
        u1 = np.sqrt((3*P)/(rho*(980*f + 6*K)))
        if (abs(u1-u0)/u0) *100 > percent_margin :
            break
        u0 = u1
    return ((np.pi*u1*D**2)/4)*3600

print(f'Flow rate for original system is {flow_solver(20/3,Delta_P, Diameter):.3f} m^3/s')
print(f'Flow rate for 20% weaker pump is {flow_solver(20/3,Delta_P*0.8, Diameter):.3f} m^3/s')
print(f'Flow rate for 10% larger diameter is {flow_solver(200/33,Delta_P, Diameter*1.1):.3f} m^3/s')



import matplotlib.pyplot as plt
from scipy.optimize import brentq
import matplotlib.ticker as mticker

Re_l=np.linspace(1,4000)
Re_t = np.logspace(np.log10(4000),np.log10(1e8),num=100)
eps_r = np.array([0.15e-3,0.3e-3, 0.075e-3,0.75e-3])/Diameter


def colebrook_solver(Re, e_r):
    def colebrook(f):
        return 1/np.sqrt(f) + 2*np.log10(e_r/3.7 + 2.51/(Re*np.sqrt(f)))
    return brentq(colebrook, 0.01, 0.1)

plt.figure()

colors = plt.cm.Set3(np.linspace(0,1,5))

plt.plot(Re_l,64/Re_l, label="64/Re")

for eps in eps_r:
    plt.loglog(Re_t,[colebrook_solver(Re,eps) for Re in Re_t], 
                label= rf"$\epsilon_{np.where(eps_r==eps)[0][0]}$",
                c=colors[np.where(eps_r==eps)[0][0]])

ax = plt.gca()
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.yaxis.get_major_formatter().set_scientific(False)   

plt.grid(which='both',visible=True)

plt.xlim(1000,1e8)
plt.ylim(0.008,0.1)

plt.xlabel('Reynolds number')
plt.ylabel('Darcy friction factor')

plt.legend()
plt.show()



from scipy.optimize import minimize, LinearConstraint

f0=0.02015
f1=0.02376
f2=0.01750
f3=0.03059

def friction_func(L):
    return abs(-f0*10 + f1*abs(L[0]) + f2*abs(L[1]) + f3*abs(L[2]))

result = minimize(friction_func,
                    [3.33,3.33,3.33],
                    method='trust-constr',
                    constraints=LinearConstraint([1,1,1],10,10))

print(result)
L1, L2, L3= [ 3.417,  6.192,  0.3901]
l1,l2,l3= [ 1.229,  7.334,  1.437]

print(-f0*10+f1*L1+f2*L2+f3*L3)
print(-f0*10+f1*l1+f2*l2+f3*l3)