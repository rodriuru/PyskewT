# -*- coding: utf-8 -*-
"""
Created on Thu May 20 17:04:20 2021

@author: rodri
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cte as C


class SkewT:

    
    def __init__(self,id=1):
        mpl.rcParams.update({'font.size': 4})
        self.fig, self.ax = plt.subplots()

    def Cp(self,w):
        return 1004*(1+0.87*w)
        
    def Lv(self,T):
        pol = np.poly1d(np.array([8.17422613e-04, -6.00650658e-02,  9.57297089e-01, -2.36500484e+03,2.50090732e+06]))
        return pol(T)
    
    def e_sat(self,T):
            return np.exp((self.Lv(T)/C.Rv)*(1/C.T_lt - 1/(T+C.T_lt)))*C.e_lt

    def w(self, T,p):
        return C.epsi*self.e_sat(T)/p
    
    def pdeTheta(self,T,Tp):
        return np.log(1e5) + C.cp*np.log(T/Tp)/C.Rd 
    
    def dibujar(self):
        ax = self.ax
        p=C.p
        eti = []
        val = []
        for i in range(6,len(p),10):
            ax.axhline(y=np.log(p[i]),linewidth=0.5, color='#ba967f')
            eti.append('{:6.0f}'.format(p[i]//100))
            val.append(np.log(p[i]))
        ax.set(ylim=(11.542,9.2103),xlim=(-40, 40),yticks=val, yticklabels=eti)
        ax.tick_params(labelright=True,labelleft=True)
        ax.yaxis.set_ticks_position('both')
        ax.set_title('SkewT')
        ax.set_xticks(np.linspace(-40,40,17))
        ax.set_xlabel('T(°C) en 1000hPa', fontsize=5)
        ax.set_ylabel('Presión(hPa)', fontsize=3)
        ax.set_aspect(aspect=40)
        plt.show()

    def T(self):
        offset = 0
        p_0 = 1e5
        T = C.T
        N = C.N
        beta=C.beta
        p = C.p
      
        for i in range(offset,N,10):
            Tp = T[i] + beta*np.log(p/p_0)
            lp = np.log(p)
            self.ax.plot(Tp,lp, color='#fffdf0', linewidth=9)    
        for i in range(offset,N,5):
            Tp = T[i] + beta*np.log(p/p_0)
            self.ax.plot(Tp,lp, color='#efffe3', linewidth=1)   
        for i in range(offset,N,1):
            Tp = T[i] + beta*np.log(p/p_0)
            self.ax.plot(Tp,lp, color='orange', linewidth=0.05)
        for i in range(offset,N,5):
            Tp = T[i] + beta*np.log(p/p_0)
            self.ax.plot(Tp,lp, color='orange', linewidth=0.1)  
        for i in range(offset,N,10):
            Tp = T[i] + beta*np.log(p/p_0)
            self.ax.plot(Tp,lp, color='orange', linewidth=0.4)
            if (T[i]>=-80)and(T[i]<=-20):
                plt.text(Tp[162],np.log(p[162]), '{:2d}'.format(int(T[i])),
                rotation=49,color='#ed8600', fontsize=4)  
                
    def ws(self,w_s=np.array([])):
        T = C.T
        if w_s.size == 0:
            w_s = np.array([0.001,0.002,0.005,0.01,0.02,0.03,0.06,0.08,0.1,0.2,
                            0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.6,2.0, 2.5,3,4,
                            5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,44,48,52])*1e-3
        
        e_satT=self.e_sat(T)
        for i in range(0,len(w_s),2):
            
            p_ws = (C.epsi/w_s[i]+1)*e_satT
          
            Tp = T + C.beta*np.log(p_ws/1e5)
            plog=np.log(p_ws)
            self.ax.plot(Tp,plog, color=C.ws_color, linewidth=0.3,linestyle='--')
            indi = np.where(np.round(100*plog)>=1145)
            indi2 = np.where(np.round(100*plog)>=935)
            if (w_s[i] > 4e-4)and(w_s[i] < 40e-3):
                indi = indi[0][0]
                plt.text(Tp[indi],plog[indi], str(int(w_s[i]*100000)/100),
                         rotation=63,color=C.ws_color, fontsize=2.5)
            if w_s[i] <= 4e-4:
                indi2=indi2[0][0]
                plt.text(Tp[indi2],plog[indi2], '{:6.2e}'.format(w_s[i]*1000),
                         rotation=56,color=C.ws_color, fontsize=2.5)
                
    def TP(self,Tpot=np.array([]), p0=1e5):
        p0 = 1e5
        # cp = self.Cp(w)
        # Rd = cp - self.Cv(w)

        if Tpot.size == 0:
            Tpot=np.arange(C.T_min,C.T_max+120,1)
        for i in range(0,len(Tpot),10):
            ppot = (C.cp/C.Rd)*np.log((C.T+C.T_lt)/(Tpot[i]+C.T_lt))+np.log(p0)
            Tp = C.T + C.beta*(ppot-np.log(1e5))
            self.ax.plot(Tp,ppot, color='#ffee00',linewidth=0.5)
        for i in range(0,len(Tpot),1):
            ppot = (C.cp/C.Rd)*np.log((C.T+C.T_lt)/(Tpot[i]+C.T_lt))+np.log(p0)
            Tp = C.T + C.beta*(ppot-np.log(1e5))
            self.ax.plot(Tp,ppot, color='#ffee00',linewidth=0.05)
            

        
    def tpe(self,unT=np.array([]),unp=1e5):
        def dpdt(p,T):
            Lv = self.Lv(T-C.T_lt)
            return (C.cp-Lv*self.w(T-C.T_lt,p)/T +
                    C.epsi*Lv*Lv*self.e_sat(T-C.T_lt)/p/C.Rv/T/T)*p/T/C.Rd
        def dpdt2(p,T):
            b=1
            w = self.w(T-C.T_lt,p)
            Lv= self.Lv(T-C.T_lt)
            A = b*C.epsi*Lv*Lv*w/C.Rd/T/T
            return p*(A+C.cp)/(Lv*w +C.Rd*T)/b
        
        def integrar(p0,T0):
            N = 500
            Tdom = np.linspace(T0,100,N)
            tp = np.empty((N,2))
            tp[0,0] = T0
            tp[0,1] = p0
            for i in range(0,N-1):
                dif = dpdt2(tp[i,1],tp[i,0])
                tp[i+1,:] = [Tdom[i+1],dif*(Tdom[i+1]-Tdom[i])+tp[i,1]]
            return tp
        
        if unT > 0:

            pyT =  integrar(unp,unT[0])
            plog = np.log(pyT[:,1])
            T = pyT[:,0]-C.T_lt + C.beta*(plog-np.log(1e5))
            self.ax.plot(T,plog, color='#008526',linewidth=1.0)
    
        TPot = np.linspace(253,313,12+1)
        for TT in TPot:
            pyT =  integrar(1e5,TT)
            plog = np.log(pyT[:,1])
            T = pyT[:,0]-C.T_lt + C.beta*(plog-np.log(1e5))
            self.ax.plot(T,plog, color='#00f044',linewidth=0.7)
            
        return plog
            

    def graficar(self,mat):
        Tw = mat[:,3]
        p = mat[:,0]*100
        TT = mat[:,2]
        Tp = TT + C.beta*np.log(p/1e5)
        Tpw = Tw + C.beta*np.log(p/1e5)
        i = np.where(p == 1e5)
        pp = self.pdeTheta(Tw[i]+273,TT[i]+273)
        self.tpe(Tw[i]+273,np.exp(pp))
        print(np.exp(pp))
        self.ax.axhline(y=pp,linewidth=0.5, color='r',linestyle='--')
        self.ax.plot(Tpw,np.log(p),linewidth=0.5, color='k')
        self.ax.plot(Tp,np.log(p),linewidth=0.5, color='k')        

    # def TPE(self,pot=np.array([]), p0=1e5):
    #     p0 = 1e5
        
    #     # cp = self.Cp(w)
    #     # Rd = cp - self.Cv(w)
    #     if pot.size == 0:
    #         pote=(C.T+C.T_lt)*np.exp(C.epsi*self.Lv(C.T)*self.e_sat(C.T)/C.cp/(C.T+C.T_lt)/p0)
    #     aux = self.Lv(C.T)*self.e_sat(C.T)/C.Rv/(C.T+C.T_lt)/p0
    #     #lnT = np.log(pote/(C.T+C.T_lt))
    #     for i in range(0,len(pote),2):
            
    #         #print()
    #         fac=1
    #         pp = (fac*(C.cp/C.Rd)*np.log((C.T+C.T_lt)/pote[i])+aux*(1+np.log(p0))+ np.log(p0))/(1+aux)
    #         #pp = np.log(C.epsi*self.Lv(C.T)*self.e_sat(C.T)/C.cp/(C.T+C.T_lt)/lnT)
    #         Tp = C.T + C.beta*(pp-np.log(1e5))
    #         self.ax.plot(Tp,pp, color='#00f044',linewidth=1)


            

m = SkewT(1)
m.T()
archivo = 'sondeo3.txt'
sondeo = np.loadtxt(archivo, skiprows=4)
m.ws()
m.TP()
m.graficar(sondeo)


#aver = m.tpe()

m.dibujar()

m.fig.savefig('aa.png', dpi=800)