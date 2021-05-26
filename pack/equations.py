import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sympy import rotations



from pack import constant


def gradient(r,element,xlabel,ylabel):
    'calculate metallicity gradient. r: radius of measured gradient.'
    'Polynomial coefficients, p, highest power first. If y was 2-D, the coefficients for k-th data set are in p[:,k].'
    'the covariance matrix, V, for which the square root of the diagonals are the estimated standard-deviation for each of the fitted coefficients.'
    'plotting the data and the its gradient plot'
    p,V=np.polyfit(r,element,1,cov=True)
    err=np.sqrt(V[0][0])
    print('d[{}]/d{} = {:.3f} +/- {:.3f} dex/kpc'.format(ylabel,xlabel,p[0],err))
    
    grad_line=np.poly1d(p)
    #plotting
    x=np.linspace(r.min(),r.max(),100)
    plt.scatter(r,element,c='k',s=1)
    plt.plot(x,grad_line(x),'r--',lw=0.5,label='d[{}]/d{} = {:.3f} +/- {:.3f} dex/kpc'.format(ylabel,xlabel,p[0],err))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best',prop={'size':8})
    plt.savefig('fig_output/gradient.png',bbox_inches='tight')
    plt.show()
    return p,V,err

def wienTemperature(lambda_max):
    'lambda max should be in meters (1 m = 1e2 cm = 1e9 nm)'
    return print('{:.3e} m ~ {:.1f} K'.format(lambda_max,constant.wien_b/lambda_max))

def radiationCurve(T):
    'T is temperature in K (for example; 5777.0 or [5777.0,6500.0,6700])'
    'plotting radiation curves given temperatures'
    wavelength_range=np.array(np.linspace(1e-9,2000.e-9,10000))
    E=np.array(constant.h_plank*constant.c/wavelength_range,dtype=np.float128)
    for i in T:
        B=((2.*constant.h_plank*constant.c**2.)/(wavelength_range**5.))/(np.exp((E)/(constant.k_boltzman*i))-1.)
        B=np.array(B,dtype=np.float128)
        plt.plot(wavelength_range*1e9,B*(1./1e13),lw=0.75,label=str(i)+' K')#wavelenght in nm and spectral energy density in W.m-3
    cmap=plt.cm.get_cmap('nipy_spectral')
    normalize = mpl.colors.Normalize(vmin=400, vmax=700)
    for i in np.arange(400,750,50):
        plt.fill_betweenx(B,i,i+50,color=cmap(normalize(i)),alpha=0.15,edgecolor=None)
    for i,j in zip([250,500,1000],['ultra-violet','visible','infrared']):
        plt.text(i,(max(B)/1e13)/2,j,rotation=90,alpha=0.5,fontsize=15)


    plt.xlabel('wavelength (nm)')
    plt.ylabel('Spectral Energy Density (10$^{13}$ W m$^{-3}$)')
    plt.ylim(0.,max(B)/1e13+0.1)
    plt.legend(loc='best')
    plt.savefig('fig_output/radiationCurve.png',bbox_inches='tight')
    plt.show()