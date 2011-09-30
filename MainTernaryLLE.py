#!/usr/bin/env python

from scipy import *
from pylab import *
from os import mkdir, path, remove, listdir
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize
import scipy.interpolate 

import glob
import tables
import vanDerWaals
import GibbsClasses
import AnalyticSystemsClasses
#import AnalyticPhaseCompSystemsClasses


class ResultsFile1(tables.IsDescription):
    T = tables.Float64Col()
    Tieline = tables.Float64Col()
    ModelParams = tables.Float64Col(shape = (6)) 
    Actual = tables.Float64Col(shape = (4))
    PureCompParams = tables.Float64Col(shape = (3))
    Converged = tables.Float64Col()

class Mixture:
        
    def __init__(self, Compounds, MixtureDataDir, PureDataDir):
        
        self.Compounds = Compounds
        self.Data = dict((Compound, {}) for Compound in Compounds)
        for Compound in Compounds:
            h5file = tables.openFile(PureDataDir+'/'+Compound+'.h5', 'r')
            properties = h5file.root.Properties
            self.Data[Compound] = dict(((field, row[field]) for row in properties.iterrows() for field in properties.colnames))        
            h5file.close()
        self.Name = '-'.join(Compounds)
        h5file = tables.openFile(MixtureDataDir+'/'+ self.Name +'.h5', 'r')
        self.M = dict(Compounds = h5file.root.Compounds.read(), Temperatures = h5file.root.ExperimentalData.TielineData.T.read())
        h5file.close()
        self.vdWaalsInstance = vanDerWaals.Slope(self.Data, Compounds)
        [Slope, VPData] = self.vdWaalsInstance.BestFit()
                    
    def Plotter(self, DWPMParams, Plane1, Plane2, Cell_s, Model, ModelInstance, Actual, c, T, Converged):

        XRange = arange(0.01, 1.00, 0.01)
        GibbsSurface = array([ModelInstance.deltaGmix(array([x1, x2]), c) for x1 in XRange for x2 in XRange])
        GradientProjection = array([ModelInstance.FirstDerivatives(array([x1, x2]), c) for x1 in XRange for x2 in XRange])
        PlotGmix = transpose(reshape(GibbsSurface, (size(XRange), -1)))
        PlotGradient = transpose(reshape(GradientProjection, (size(XRange), -1)))
        x_1, x_2 = meshgrid(XRange, XRange)
        
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        ax = Axes3D(fig)

        Gx, Gy = np.gradient(PlotGmix) # gradientS with respect to x and y
        G = (Gx**2+Gy**2)**.5  # gradient magnitude
        N = G/G.max()  # normalie 0..1
        #for x2 in range(len(XRange)):
        #    for x1 in range(len(XRange)):
        #        if ModelInstance.FirstDerivatives(array([x1, x2]), c)==1:
        #            SurfaceColours[x1, x2] = 1
        #        else:
        #            SurfaceColours[x1, x2] = 0
        
        ax.plot_surface(x_1, x_2, PlotGmix, rstride=1, cstride=1, cmap = matplotlib.cm.jet, linewidth = 1.0, antialiased =False)        
        
        Tieline1x = [Actual[0, 0], Actual[0, 2]]
        Tieline1y = [Actual[0, 1], Actual[0, 3]]
        Tieline2x = [Actual[1, 0], Actual[1, 2]]
        Tieline2y = [Actual[1, 1], Actual[1, 3]]

        ax.plot(Tieline1x, Tieline1y, zs = -1, zdir ='z', linewidth = 1.5, color = 'k')
        ax.plot(Tieline2x, Tieline2y, zs = -1, zdir ='z', linewidth = 1.5, color = 'k')
        
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+self.Compounds[0].capitalize(), fontsize = 14)
        matplotlib.pyplot.ylabel(r'Mole Fraction of '+self.Compounds[1].capitalize(), fontsize = 14)
        #matplotlib.pyplot.zlabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$', fontsize = 14)
        matplotlib.pyplot.title(r"Predicted Gibbs Surface ", fontsize=14) 
        fig.savefig('Results/'+self.Name+'/'+Model+'/PredictedGibbsWireframe.png')

        #matplotlib.pyplot.xlabel(r'Mole Fraction of '+self.Compounds[0].capitalize(), fontsize = 14)
        #matplotlib.pyplot.ylabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$', fontsize = 14)
        #matplotlib.pyplot.title(r'\textbf{Predicted $\displaystyle\frac{\Delta G_{mix}}{RT}$ %3.2f}'%T, fontsize = 14)
        #ax = fig.add_axes([0,0,1,1])
        #ax.text(0,0, 'Model: %s, Params: %0.4E, %0.4E, %.2f, %.2f, %.0f'%(Model, Params[0], Params[1], Cell_s[0], Cell_s[1], Converged), fontsize=12, transform=ax.transAxes)
        #ax.set_axis_off()
        
        show()
        
        #matplotlib.pyplot.close()
        
        #if not(path.exists('Results/'+self.Name+'/'+ Model)):
        #    mkdir('Results/'+self.Name+'/'+ Model)

#        if not(path.exists('Results/'+self.Name+'/'+ Model+'/'+self.Name +'.h5')):
#            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+self.Name +'.h5', 'w', "Optimization Outputs")
#            table = h5file.createTable("/", "Outputs", ResultsFile1, "Optimal model parameters, predicted phase equilibrium, errors etc")
#        else:
#            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+self.Name +'.h5', 'r+')
#            table = h5file.root.Outputs

        #table.row['T'] = T
        #table.row['ModelParams'] = reshape(append(Params, Cell_s),-1)
        #table.row['Actual'] = Actual
        #table.row['PureCompParams']  = reshape(array(c),-1)
        #table.row['Converged'] = Converged
        #table.row.append()
        #table.flush()
        h5file.close()
           
    def ParamCalc(self, ModelInstance, InitParamsLimit, Bounds, R, System, Cell_s, MixtureDataDir):

        if not(path.exists('Results/'+self.Name+'/'+ Model)):
            mkdir('Results/'+self.Name+'/'+ Model)

        if path.exists('Results/'+self.Name+'/'+Model+'/'):
            fileList = listdir('Results/'+self.Name+'/'+Model+'/')
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model+'/', fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model+'/')
                        
        InitParams = (InitParamsLimit[1]-InitParamsLimit[0])*random(6) + reshape(array([InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0]]), -1)
        InitParamj = append(InitParams, array([-2, -2., 0., -2, -2., 0.])+ (4., 4., -1., 4., 4., -1.)*random(6))

        for T in self.M['Temperatures']:
            
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]

            h5file = tables.openFile(MixtureDataDir+'/'+ self.Name +'.h5', 'r')
            exec('ActualPhase1 = h5file.root.ExperimentalData.TielineData.T' + str(int(where(self.M['Temperatures']==T)[0]) +1) +'.Phase1.read()')
            exec('ActualPhase2 = h5file.root.ExperimentalData.TielineData.T' + str(int(where(self.M['Temperatures']==T)[0]) +1) +'.Phase2.read()')

                #random_integers(size(self.M['Temperatures'])

            Tieline1 = 0
            Tieline2 = size(ActualPhase1[:,0])-1

            Actual_a1 = ActualPhase1[Tieline1, :]
            Actual_b1 = ActualPhase2[Tieline1, :]

            Actual_a2 = ActualPhase1[Tieline2, :]
            Actual_b2 = ActualPhase2[Tieline2, :]    

            Actual1 = append(Actual_a1[:2],Actual_b1[:2])
            Actual2 = append(Actual_a2[:2],Actual_b2[:2])
            Actual = reshape(append(Actual1, Actual2), (-1, size(Actual1)))

            #print Actual
       
            MaxFEval = 1000
            TolParam = 1e-6
            Converge = 3
            NStarts = 1
            MaxNStarts = 10000
            
            while (Converge!=1 and NStarts<MaxNStarts):
                                                           
                [Paramj, infodict, Converge, Mesg] = scipy.optimize.fsolve(System.SystemEquations, InitParamj, (c, Actual, R, T), System.SystemEquationsJac, 1, 0, TolParam, MaxFEval, None, 0.0, 100, None)
                print Mesg
                NFCalls = infodict['nfev']
                NJCalls = infodict['njev']
                NStarts = NStarts + 1
                                                  
                if (Paramj[8]>0.01) or (Paramj[11]>0.01): #or not(all(Bounds[0]<Paramj[0:6], 0) and all(Paramj[0:6]<Bounds[1], 0)): 
                    Converge = 3

                InitParams = (InitParamsLimit[1]-InitParamsLimit[0])*random(6) + reshape(array([InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0]]), -1)
                InitParamj = append(InitParams, array([-2, -2., 0., -2, -2., 0.])+ (4., 4., -1., 4., 4., -1.)*random(6))
                                                                                                        
            # print Mesg
            ##print NStarts
            
            DWPMParams = Paramj[0:6]
            Plane1 = Paramj[6:9]
            Plane2 = Paramj[9:]
            
            ParamInstance = ModelInstance(append(DWPMParams, Cell_s))
            self.Plotter(DWPMParams, Plane1, Plane2, Cell_s, Model, ParamInstance, Actual, c, T, Converge)
             
#    def PhaseDiagramCalc(self, Model, System, Params, R, T):

                                   
##=============================================================##
Model = 'DWPM' 
ModelInstance = GibbsClasses.DWPM
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
Bounds = (-100, 0.00)
InitParamsLimit =(-0.1500, 0.00)
alpha = 0.2
z = 10
R = 8.314

if not(path.exists('Results/')):
    mkdir('Results/')   

for file in listdir(MixtureDataDir):

    h5file = tables.openFile(MixtureDataDir+'/'+file, 'r')
    Compounds = h5file.root.Compounds.read()
    h5file.close()

    Optimization = Mixture(Compounds, MixtureDataDir, PureDataDir)

    if not(path.exists('Results/'+Optimization.Name)):
        mkdir('Results/'+Optimization.Name)
    
    Cell_s = (0.5, 0.5, 0.5)
    System = AnalyticSystemsClasses.SystemDWPM(Cell_s)
   # PhaseSystem = AnalyticPhaseCompSystemsClasses.SystemDWPM(Cell_s)
    
    Optimization.ParamCalc(ModelInstance, InitParamsLimit, Bounds, R, System, Cell_s, MixtureDataDir)

        
    
