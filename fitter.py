import numpy as np # numerical python--handles arrays
import scipy # scientific python--scientific computing package
import scipy.optimize # optimization package within scipy -- use for curve fitting
import matplotlib.pyplot as plt # plotting package with all the same functions as matlab
from mpl_toolkits.mplot3d import Axes3D
import math

# This program fits all the IV curves (Ids vs Vds) with the same parameter to the Spice JFET model by minimizing chi2.
# To use it you just need to change the path and the name of your data file and eventually the starting parameters


filename = 'Run2_Hemt_1_31May2013_12h30'
with open(filename+'.txt', 'r') as f:
    first_line = f.readline()
    second_line = f.readline()
Nmax = second_line.split()
ke2220Nmax = int(Nmax[0])
ke2401Nmax = int(Nmax[1])
NumPerHEMT = ke2220Nmax*ke2401Nmax
data = np.loadtxt(open(filename+'.txt'),skiprows=2)

HEMTID = ['1A','1B','1C','1D','2A','2B','2C','2D','3A','3B','3C','3D','4A','4B','4C','4D','5A','5B','5C','5D']

DoF = NumPerHEMT - 3 # number of data points - number of free parameters
p=[6.2e-5,-207,1e-5] # p gives the starting parameters for scipy.optimize.curve_fit in form [beta,Vt,lambda]

# function to calculate chi2
#def chi2(p, Vds, Vgs, Ids, Idserr):
#
#        chi2 = 0
#
#        fit = np.zeros(np.size(Vds))
#        for i in range(np.size(Vds)):
#            if Vgs[i]-p[1]<=0:
#                fit[i]=0
#            elif Vds[i]>=Vgs[i]-p[1]:
#                fit[i] = p[0]*(Vgs[i]-p[1])**2*(1+p[2]*Vds[i])
#            else:
#                fit[i] = p[0]*Vds[i]*(2*(Vgs[i]-p[1])-Vds[i])*(1+p[2]*Vds[i])
#
#            delta = fit[i]-Ids[i]
#            chi2 += (delta/Idserr[i])**2
#
#        return chi2

def _f(Vgs, Vds, beta, vthr ,lamb):
    if Vgs-vthr<=0:
        return 0
    elif Vds>=Vgs-vthr:
        return beta*(Vgs-vthr)**2*(1+lamb*Vds)
    else:
        return beta*Vds*(2*(Vgs-vthr)-Vds)*(1+lamb*Vds)

f = np.vectorize(_f, excluded=["beta", "vthr", "lamb"])

def f2(args, beta, vthr, lamb):
    return f(args[0], args[1], beta, vthr, lamb)

# run once for each HEMT
for x in range(0,1):
    # read in this HEMT's data
    Vgs = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),0]*(-1)
    Vds = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),1]
    cut = Vgs > -85

    Ids = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),2]

    Vgs, Vds, Ids = Vgs[cut], Vds[cut], Ids[cut]
    Vgserr = abs(Vgs*0.0003)+10                          # calculate Vgs uncertianty using the ke 2220 datasheet info on voltage accuracy
    Vdserr = np.zeros(Vds.shape)
    Idserr = np.zeros(Ids.shape)

    for g in range(len(Ids)):
        if Vds[g]<=211:
            Vdserr[g] = abs(Vds[g]*0.00012)+0.3          # calculate Vds uncertianty using the ke 2401 datasheet info on voltage ranges and accuracy
        elif Vds[g]<=2110:
            Vdserr[g] = abs(Vds[g]*0.00012)+0.3
        elif Vds[g]<=21100:
            Vdserr[g] = abs(Vds[g]*0.00015)+1.5
        if Ids[g]<=0.00105:
            Idserr[g] = abs(Ids[g]*0.00035)+0.0000006+0.000000005    # calculate Ids uncertianty using the ke 2401 datasheet info on current ranges and accuracy and noise
        elif Ids[g]<=0.0105:
            Idserr[g] = abs(Ids[g]*0.00033)+0.0000020+0.0000050
        elif Ids[g]<=0.105:
            Idserr[g] = abs(Ids[g]*0.00031)+0.0000200+0.0000500
        elif Ids[g]<=1.05:
            Idserr[g] = abs(Ids[g]*0.00034)+0.0002000+0.0005000
        elif Ids[g]<=10.5:
            Idserr[g] = abs(Ids[g]*0.00045)+0.0020000+0.0500000
        elif Ids[g]<=105:
            Idserr[g] = abs(Ids[g]*0.00066)+0.0200000+0.0010000
        elif Ids[g]<=1050:
            Idserr[g] = abs(Ids[g]*0.00270)+0.9000000+0.1000000

    VgsCount = [s for s in Vgs if s >= -120]
    Count = math.ceil(len(VgsCount)/float(ke2401Nmax))
    color=iter(plt.cm.rainbow(np.linspace(0,1,Count)))

    # do the minimization for real tho
    xdata=np.vstack((Vgs,Vds))
    ydata=Ids
    vh = scipy.optimize.curve_fit(f2,xdata,ydata,p0=p,sigma=Idserr)

    # get best-fit parameters and errors on those parameters
    v = vh[0]
    print v
    beta_opt = float(v[0])
    Vt_opt = float(v[1])
    lamb_opt = float(v[2])

    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(Vgs, Vds, zs=Ids, zdir='z', s=20, c='r', alpha=0.5)
    X = np.linspace(0,-250, num=200)
    Y = np.linspace(0,250, num=200)
    #XY = np.array([[X],[Y]])
    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, f(X, Y, beta_opt, Vt_opt ,lamb_opt), alpha=0.3)
    plt.show()

    plt.figure(2)
    plt.xlim([0,250])

    for j in range(ke2220Nmax):
        Vgsp = np.zeros(ke2401Nmax)
        Vgsperr = np.zeros(ke2401Nmax)
        Vdsp = np.zeros(ke2401Nmax)
        Vdsperr = np.zeros(ke2401Nmax)
        Idsp = np.zeros(ke2401Nmax)
        Idsperr = np.zeros(ke2401Nmax)
        fitFuncp = np.zeros(ke2401Nmax)
        for k in range(ke2401Nmax):
            Vgsp[k] = data[(NumPerHEMT*x)+(ke2401Nmax*j)+k,0]*(-1)
            Vdsp[k] = data[(NumPerHEMT*x)+(ke2401Nmax*j)+k,1]
            Idsp[k] = data[(NumPerHEMT*x)+(ke2401Nmax*j)+k,2]
            Vgsperr[k] = Vgserr[(ke2401Nmax*j)+k]
            Vdsperr[k] = Vdserr[(ke2401Nmax*j)+k]
            Idsperr[k] = Idserr[(ke2401Nmax*j)+k]
            if Vgsp[k]-Vt_opt<=0:
                fitFuncp[k] = 0
            elif Vdsp[k]>=Vgsp[k]-Vt_opt:
                fitFuncp[k] = beta_opt*(Vgsp[k]-Vt_opt)**2*(1+lamb_opt*Vdsp[k])
            else:
                fitFuncp[k] = beta_opt*Vdsp[k]*(2*(Vgsp[k]-Vt_opt)-Vdsp[k])*(1+lamb_opt*Vdsp[k])
        if Vgsp[1]>=-120:
            c=next(color)
            plt.plot(Vdsp,fitFuncp ,'r',linewidth=3)
            plt.errorbar(Vdsp,Idsp,yerr=Idsperr,xerr=Vdsperr,fmt='b.',c=c,label=str(round(Vgsp[1],3))+'$\pm$'+str(round(Vgsperr[1],3)))

        plt.legend(title='$V_{gs}$[mV]',bbox_to_anchor=(1.15, 1))

    plt.title('Fit to $I_{ds}$ vs $V_{ds}$ for all $V_{gs}$ $\geq$ -100mV ' + HEMTID[x])
    plt.ylabel('$I_{ds}$[mA]')
    plt.xlabel('$V_{ds}$[mV]')
#    paramStringChi2 = ' $\chi^2$ = '+str(round(min_chi2,2))
    paramStringBeta = r'$\beta$ = '+str(round(beta_opt,8))+' mA/m$V^2$'
    paramStringVt = '$V_{threshold}$ = '+str(round(Vt_opt,3))+' mV'
    paramStringLambda = '$\lambda$ = '+str(round(lamb_opt,6))+' $\Omega$'
#    plt.figtext(.2,.2,paramStringChi2)
    plt.figtext(.2,.15,paramStringBeta)
    plt.figtext(.6,.2,paramStringVt)
    plt.figtext(.6,.15,paramStringLambda)
    plt.savefig(filename+'optimizeChi2_'+HEMTID[x]+'.png')
    plt.close()
    color=iter(plt.cm.rainbow(np.linspace(0,1,Count)))
