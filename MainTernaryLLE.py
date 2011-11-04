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
import AnalyticBinarySystems



class ResultsFile1(tables.IsDescription):
    T = tables.Float64Col()
    PureCompParams = tables.Float64Col(shape = (3))
    ModelParams = tables.Float64Col(shape = (6)) 
    SParameters = tables.Float64Col(shape = (3))
    Converged = tables.Float64Col()

class TielineData(tables.IsDescription):
    TielineCompositions = tables.Float64Col(shape=(4))
    PerpVectorMN = tables.Float64Col(shape=(4))
    TangentPlane = tables.Float64Col(shape=(3)) 

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
        self.BinaryLLE  = h5file.root.BinaryLLE.read()
        h5file.close()
        self.vdWaalsInstance = vanDerWaals.Slope(self.Data, Compounds)
        [Slope, VPData] = self.vdWaalsInstance.BestFit()
                    
    def Plotter(self, DWPMParams, Plane1, Plane2, Cell_s, Model, ModelInstance, Actual, c, T, Converged):

        print "Plotter"

        XRange = arange(0.01, 1.00, 0.01)
        GibbsSurface = array([ModelInstance.deltaGmix(array([x1, x2]), c) for x1 in XRange for x2 in XRange])
        PlotGmix = transpose(reshape(GibbsSurface, (size(XRange), -1)))
        OffSet =nanmin(PlotGmix)-0.1
        HessianProjection = array([ModelInstance.Hessian(array([x1, x2]), c) for x1 in XRange for x2 in XRange])
        HessianProjection[where(HessianProjection == 0)] = NaN
        HessianProjection[where(HessianProjection == 1)] = OffSet
        PlotProjection = transpose(reshape(HessianProjection, (size(XRange), -1)))
        x_1, x_2 = meshgrid(XRange, XRange)
        
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        fig.suptitle(r'Predicted'+r'$\displaystyle\frac{\Delta G_{mix}}{RT}$'+'Surface', fontsize=14) 
        ax1 = fig.add_subplot(121)#CS, projection= '3d')
        ax2 = fig.add_subplot(122)#CS, projection= '3d')
        
        ax1.plot_wireframe(x_1, x_2, PlotGmix, rstride=1, cstride=1, linewidth = 0.15, antialiased =True)        
        ax1.plot_surface(x_1, x_2, PlotProjection, rstride =1, cstride=1, linewidth = 0.1, antialiased =True, color = "r", alpha = 0.8)
        
        Tielines1x1 = [Actual[0, 0], Actual[0, 2]]
        Tielines2x1 = [Actual[1, 0], Actual[1, 2]]
        Tielines1x2 = [Actual[0, 1], Actual[0, 3]]
        Tielines2x2 = [Actual[1, 1], Actual[1, 3]]

        ax1.plot(Tielines1x1,
                 Tielines1x2,
                 #zs=OffSet, zdir='z',
                 linewidth = 1.0, color = 'k')
        ax1.plot(Tielines2x1, Tielines2x2,
                 #zs=OffSet, zdir='z',
                 linewidth = 1.0, color = 'k')
        
        ax1.set_xlabel(r'Mole Fraction '+self.Compounds[0].capitalize(), fontsize = 14)
        ax1.set_ylabel(r'Mole Fraction '+self.Compounds[1].capitalize(), fontsize = 14)

        ax2.plot_wireframe(x_1, x_2, PlotGmix, rstride=1, cstride=1, linewidth = 0.15, antialiased =True)        
        ax2.plot_surface(x_1, x_2, PlotProjection, rstride =1, cstride=1, linewidth = 0.1, antialiased =True, color = "r", alpha = 0.8)
        
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        CS = ax2.contour(x_1, x_2, PlotGmix, 10,
                         #zdir='z', offset=OffSet,
                         colors = "k")
        matplotlib.pyplot.clabel(CS, fontsize=14, inline=1)

        ax2.set_xlabel(r'Mole Fraction '+self.Compounds[0].capitalize(), fontsize = 14)
        ax2.set_ylabel(r'Mole Fraction '+self.Compounds[1].capitalize(), fontsize = 14)
       
        fig.savefig('Results/'+self.Name+'/'+Model+'/PredictedGibbsWireframe.png')

        #ax.text(0,0, 'Model: %s, Params: %0.4E, %0.4E, %.2f, %.2f, %.0f'%(Model, Params[0], Params[1], Cell_s[0], Cell_s[1], Converged), fontsize=12, transform=ax.transAxes)
        #ax.set_axis_off()
                        
        matplotlib.pyplot.close()
        
        if not(path.exists('Results/'+self.Name+'/'+ Model)):
            mkdir('Results/'+self.Name+'/'+ Model)

        if not(path.exists('Results/'+self.Name+'/'+ Model+'/'+self.Name +'.h5')):
            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+self.Name +'.h5', 'w', "Optimization Outputs")
            table = h5file.createTable("/", "Outputs", ResultsFile1, "Calculated model parameters")
        else:
            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+self.Name +'.h5', 'r+')
            table = h5file.root.Outputs
        
        table.row['T'] = T
        table.row['PureCompParams'] = c
        table.row['ModelParams'] = DWPMParams
        table.row['SParameters'] = Cell_s
        table.row['Converged'] = Converged
        table.row.append()
        table.flush()
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
            h5file.close()
                
            Tieline1 = 0
            Tieline2 = size(ActualPhase1[:,0])-1

            Actual_a1 = ActualPhase1[Tieline1, :]
            Actual_b1 = ActualPhase2[Tieline1, :]

            Actual_a2 = ActualPhase1[Tieline2, :]
            Actual_b2 = ActualPhase2[Tieline2, :]    

            Actual1 = append(Actual_a1[:2],Actual_b1[:2])
            Actual2 = append(Actual_a2[:2],Actual_b2[:2])
            Actual = reshape(append(Actual1, Actual2), (-1, size(Actual1)))
                  
            MaxFEval = 1000
            TolParam = 1e-6
            Converge = 3
            NStarts = 1
            MaxNStarts = 10000
            
            while (Converge!=1 and NStarts<MaxNStarts):
                                                           
                [Paramj, infodict, Converge, Mesg] = scipy.optimize.fsolve(System.SystemEquations, InitParamj, (c, Actual, R, T), System.SystemEquationsJac, 1, 0, TolParam, MaxFEval, None, 0.0, 100, None)
                NFCalls = infodict['nfev']
                NJCalls = infodict['njev']
                NStarts = NStarts + 1
                                                  
                if (Paramj[8]>0.01) or (Paramj[11]>0.01): #or not(all(Bounds[0]<Paramj[0:6], 0) and all(Paramj[0:6]<Bounds[1], 0)): 
                    Converge = 3

                InitParams = (InitParamsLimit[1]-InitParamsLimit[0])*random(6) + reshape(array([InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0], InitParamsLimit[0]]), -1)
                InitParamj = append(InitParams, array([-2, -2., 0., -2, -2., 0.])+ (4., 4., -1., 4., 4., -1.)*random(6))
            
            DWPMParams = Paramj[0:6]
            Plane1 = Paramj[6:9]
            Plane2 = Paramj[9:]
            
            ParamInstance = ModelInstance(append(DWPMParams, Cell_s))
            self.Plotter(DWPMParams, Plane1, Plane2, Cell_s, Model, ParamInstance, Actual, c, T, Converge)
             
    def PhaseDiagramCalc(self, ModelInstance, Params, PureParams, Cell_s, R, T):
        
        h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/PhaseDiagramData.h5', 'w', "Calculated Tielines")
                           
        MixturePairs = ['12', '23', '31']
        
        for Pair in arange(size(self.BinaryLLE)):
        #BinaryLLE -->Boolean array with pairs[1-2, 2-3, 3-1]
            if self.BinaryLLE[Pair]:
                
                table = h5file.createTable("/", "Tielines", TielineData, "Calculated tielines starting from Binary Pair "+MixturePairs[Pair])
                
                BinaryCell_s = Cell_s[Pair]
                
                if Pair<2:
                    BinaryPureCompParams = PureParams[Pair:Pair+2]
                    BinaryParams = Params[Pair*4:Pair*4+2]/BinaryPureCompParams
                else:
                    BinaryPureCompParams = PureParams[0:3:2]
                    BinaryParams = Params[Pair:Pair+2]/BinaryPureCompParams
                    BinaryParams = BinaryParams[::-1]
                
                BinaryPhaseSystem = AnalyticBinarySystems.SystemDWPM(BinaryCell_s)
            
                MaxFEval = 10000
                TolPhase = 1e-3
                ConvergePhase = 3
                NStarts = 1
                MaxNStarts = 10000
        
                while (ConvergePhase!=1 and NStarts<MaxNStarts):
                    
                    InitPhase = append(array([0.1, -0.1])*random(2)+array([0.0, 1.0]), array([-2., -3]))
                    
                    [Phase, infodictPhase, ConvergePhase, MesgPhase] = scipy.optimize.fsolve(BinaryPhaseSystem.SystemEquations, InitPhase, (BinaryParams, R, T), BinaryPhaseSystem.SystemEquationsJac, 1, 0, TolPhase, MaxFEval, None, 0.0, 100, None)
            
                    NFCalls = infodictPhase['nfev']
                    NJCalls = infodictPhase['njev']
                    NStarts = NStarts + 1
                    
                    if (Phase[3]>0.001) or ((Phase[2]*1+Phase[3])>0.001) or not(0.0 <= Phase[0]<=1.0) or not(0.0 <= Phase[1]<= 1.0) or (abs(Phase[0]-Phase[1])<0.1):
                        ConvergePhase = 3
                             
                print Phase
                                                              
                System = AnalyticSystemsClasses.SystemDWPM(Params, Cell_s, PureParams)
        
                Tieline = zeros(6)
                if Pair==0:
                    x_alpha = array([Phase[0], 1-Phase[0]])
                    x_Beta = array([Phase[1], 1-Phase[1]])        
                elif Pair==1:
                    x_alpha = array([0.0001 ,Phase[0]])
                    x_Beta = array([0.0001, Phase[1]])        
                else:
                    x_alpha = array([1-Phase[0], 0.0001])
                    x_Beta = array([1-Phase[1], 0.0001])
            
                nInHessianRegion = True
                ParamInstance = ModelInstance(append(Params, Cell_s))
                TangentPlaneArray = zeros(3)
                
                
                while (linalg.norm(x_Beta-x_alpha)>0.01) and nInHessianRegion:
                    print 'Alpha_0:', x_alpha
                    print 'Beta_0:', x_Beta
                    print

                    m = x_alpha + (x_Beta - x_alpha)/2.0
                    stepsize = 0.05
                     # Check direction of step for initial binary mixture 2-3 => might have to be -?
                    n = System.NewN(x_alpha, x_Beta, m, stepsize)

                    table.row['TielineCompositions'] = append(x_alpha, x_Beta)
                    table.row['PerpVectorMN'] = append(m, n)
                    table.row['TangentPlane'] = TangentPlaneArray
                    table.row.append()
                    table.flush()
                    
                    if not(ParamInstance.Hessian(n, PureParams)):
                        nInHessianRegion = False
                        print 'Initial point outside Hessian region'

                             
                    Equilibrium = False
                                        
                    while (not(Equilibrium)) and nInHessianRegion:
                        
                        #InitComps = zeros(4)
                        #InitComps[0:3:2] = random(2)
                        #InitComps[1:4:2] = 1-InitComps[0:3:2]*random(2)
                        
                        InitComps = append(x_alpha, x_Beta)
                        
                        InitPhaseSplit = append(InitComps, array([-2, -2., 0.]) + (4., 4., -1.)*random(3))
                        [PhaseSplit, infodict3D, Converge3D, Mesg3D] = scipy.optimize.fsolve(System.SystemEquations, InitPhaseSplit, (n,), System.SystemEquationsJac, 1, 0, TolPhase, MaxFEval, None, 0.0, 100, None)                        

                        #print PhaseSplit
                        x_alpha_t = PhaseSplit[0:2]
                        x_Beta_t = PhaseSplit[2:4]

                        #t = (n[1] - x_alpha_t[1])/(x_Beta_t[1]-x_alpha_t[1])
                        #g_m = 
                        #g_alpha = 
                        #g_Beta = 
                        #LineTest = abs(f_alpha +(f_Beta - f_alpha)*t - f_m)
                        
                        if (Converge3D == 1) \
                           and (sum(x_alpha_t) <= 1.0) \
                           and (sum(x_Beta_t) <= 1.0) \
                           and (all(PhaseSplit[0:4] >= 0)) \
                           and (linalg.norm(x_alpha_t - x_Beta_t) > 0.01):
                            Equilibrium = True
                            x_alpha = x_alpha_t
                            x_Beta = x_Beta_t
                            TangentPlaneArray = PhaseSplit[4:7]
                           
                                
        print linalg.norm(x_Beta-x_alpha)
        h5file.close()    

        XRange = arange(0.001, 1.00, 0.01)
        GibbsSurface = array([ParamInstance.deltaGmix(array([x1, x2]), PureParams) for x1 in XRange for x2 in XRange])
        PlotGmix = transpose(reshape(GibbsSurface, (size(XRange), -1)))
        OffSet =nanmin(PlotGmix)-0.1
        #HessianProjection = array([ParamInstance.Hessian(array([x1, x2]), PureParams) for x1 in XRange for x2 in XRange])
        #HessianProjection[where(HessianProjection == 0)] = NaN
        #HessianProjection[where(HessianProjection == 1)] = OffSet
        #PlotProjection = transpose(reshape(HessianProjection, (size(XRange), -1)))
        x_1, x_2 = meshgrid(XRange, XRange)
        
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        fig.suptitle(r'Calculated Tielines', fontsize=14) 
        ax = fig.add_subplot(111)#,projection= '3d')
        
        #ax.plot_wireframe(x_1, x_2, PlotGmix, rstride=1, cstride=1, linewidth = 0.15, antialiased =True)        
        #ax.plot_surface(x_1, x_2, PlotProjection, rstride =1, cstride=1, linewidth = 0.1, antialiased =True, color = "r", alpha = 0.8)

        h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/PhaseDiagramData.h5', 'r')

        TielinesPlot = array([row['TielineCompositions'] for row in h5file.root.Tielines.iterrows()])
        PerpVectorsPlot = array([row['PerpVectorMN'] for row in h5file.root.Tielines.iterrows()])
        TangentPlanes = array([row['TangentPlane'] for row in h5file.root.Tielines.iterrows()])
        
        for row in arange(TielinesPlot.shape[0]):
            #ax.plot(TielinesPlot[row,0:3:2], TielinesPlot[row,1:4:2], array([ParamInstance.deltaGmix(TielinesPlot[row,0:2], PureParams),ParamInstance.deltaGmix(TielinesPlot[row,2:4], PureParams)]), zdir='z', linewidth = 1.0, color = 'k')
       
            
            ax.plot(TielinesPlot[row,0:3:2],
                    TielinesPlot[row,1:4:2],
                    linewidth = 1.0, color = 'k')
            ax.axis('equal')
            ax.plot(PerpVectorsPlot[row,0:3:2],
                    PerpVectorsPlot[row,1:4:2],
                    linewidth = 1.0, color = 'r')
               
        ax.set_xlabel(r'Mole Fraction '+self.Compounds[0].capitalize(), fontsize = 14)
        ax.set_ylabel(r'Mole Fraction '+self.Compounds[1].capitalize(), fontsize = 14)
        fig.savefig('Results/'+self.Name+'/'+Model+'/CalculatedPhaseDiagram.png')
        matplotlib.pyplot.show()
        
                
        h5file.close()
        

                    

            
                                   
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
    
    ### Parameter Calculation
       
    Cell_s = (0.5, 0.5, 0.5)
    #System = AnalyticSystemsClasses.SystemDWPM(Cell_s)
       
    #Optimization.ParamCalc(ModelInstance, InitParamsLimit, Bounds, R, System, Cell_s, MixtureDataDir)

    ### Phase Diagram Calcualtion
    
    h5file = tables.openFile('Results/'+Optimization.Name+'/'+Model+'/'+Optimization.Name +'.h5', 'r')
    Temps = array([row['T'] for row in h5file.root.Outputs.iterrows()])
    ModelParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
    PureParams = array([row['PureCompParams'] for row in h5file.root.Outputs.iterrows()])
    SParameters = array([row['SParameters'] for row in h5file.root.Outputs.iterrows()])
    h5file.close()        
    
    for i in arange(size(Temps)):
       
        Optimization.PhaseDiagramCalc(ModelInstance, ModelParams[i], PureParams[i], SParameters[i], R, Temps[i])
        

        
    
