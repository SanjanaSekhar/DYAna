from HiggsAnalysis.CombinedLimit.PhysicsModel import *
 
#this version measures AFB and gauge boson polarization (alpha) and leptoquark coupling
#includes LQ
class DY_AFB(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6, -1.0,1.0]");
        self.modelBuilder.doVar("A0[0.05, -1.0, 1.0]");
        self.modelBuilder.doVar("yLQ[0.01, -5.0, 5.0]");
        self.modelBuilder.doSet("POI","A0,Afb,yLQ")

      
        self.modelBuilder.factory_('expr::Alph("2.0*@0/(2.0-@0)",A0)')
        self.modelBuilder.factory_('expr::Norm("3.0/4.0/(2.0+@0)",Alph)')
        self.modelBuilder.factory_('expr::RAlph("@0*@1",Alph,Norm)')
        self.modelBuilder.factory_('expr::Rpl("(@0+@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::Rmn("(@0-@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::yLQ2("(@0*@0)",yLQ)')
        self.modelBuilder.factory_('expr::yLQ4("(@0*@0*@0*@0)",yLQ)')


 
 
 
 
    def getYieldScale(self,bin,process):
        if 'alpha' in process: return "RAlph"
        elif 'pl' in process: return "Rpl"
        elif 'mn' in process : return "Rmn"
        elif 'LQpure_d' in process: return "yLQ4"
        elif 'LQint_d' in process: return "yLQ2"
        elif 'LQpure_u' in process: return "yLQ4"
        elif 'LQint_u' in process: return "yLQ2"
        else:
            #print("Didnt find process %s bin %s in specifications \n" % (process, bin))
            return 1

#This version includes constraining fakes between ss and os regions
class DY_AFB_ss(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.7,0.7]");
        self.modelBuilder.doVar("A0[0.05, -1.0, 1.0]");
        self.modelBuilder.doSet("POI","Afb,A0")

        # ss templates
        self.modelBuilder.doVar("R_ee_os_fakes[0.6,0.0,1.0]");
        self.modelBuilder.doVar("ee16_fakes_norm[1.0, 0.01, 10.]");
        self.modelBuilder.doVar("ee17_fakes_norm[1.0, 0.01, 10.]");
        self.modelBuilder.doVar("ee18_fakes_norm[1.0, 0.01, 10.]");
        #Remember, cant use spaces in these formulas!
        #self.modelBuilder.options.verbose = 10
        self.modelBuilder.factory_('expr::R_ee16_qcd_os("@0*@1",ee16_fakes_norm,R_ee_os_fakes)')
        self.modelBuilder.factory_('expr::R_ee17_qcd_os("@0*@1",ee17_fakes_norm,R_ee_os_fakes)')
        self.modelBuilder.factory_('expr::R_ee18_qcd_os("@0*@1",ee18_fakes_norm,R_ee_os_fakes)')
        self.modelBuilder.factory_('expr::R_ee16_qcd_ss("@0*(1.0-@1)",ee16_fakes_norm,R_ee_os_fakes)')
        self.modelBuilder.factory_('expr::R_ee17_qcd_ss("@0*(1.0-@1)",ee17_fakes_norm,R_ee_os_fakes)')
        self.modelBuilder.factory_('expr::R_ee18_qcd_ss("@0*(1.0-@1)",ee18_fakes_norm,R_ee_os_fakes)')
      
        self.modelBuilder.factory_('expr::Alph("2.0*@0/(2.0-@0)",A0)')
        self.modelBuilder.factory_('expr::Norm("3.0/4.0/(2.0+@0)",Alph)')
        self.modelBuilder.factory_('expr::RAlph("@0*@1",Alph,Norm)')
        self.modelBuilder.factory_('expr::Rpl("(@0+@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::Rmn("(@0-@1)",Norm,Afb)')



 
 
 
 
    def getYieldScale(self,bin,process):
        if 'alpha' in process: return "RAlph"
        elif 'pl' in process: return "Rpl"
        elif 'mn' in process : return "Rmn"

        elif 'qcd' in process and 'ee16_ss' in bin: return "R_ee16_qcd_ss"
        elif 'qcd' in process and 'ee17_ss' in bin: return "R_ee17_qcd_ss"
        elif 'qcd' in process and 'ee18_ss' in bin: return "R_ee18_qcd_ss"
        elif 'qcd' in process and 'ee16' in bin: return "R_ee16_qcd_os"
        elif 'qcd' in process and 'ee16' in bin: return "R_ee17_qcd_os"
        elif 'qcd' in process and 'ee16' in bin: return "R_ee18_qcd_os"
        else:
            #print("Didnt find process %s bin %s in specifications \n" % (process, bin))
            return 1
 
### This is class measures AFB only
class DY_AFB_LO(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.70,0.70]");
        self.modelBuilder.doSet("POI","Afb")

        # ss templates
        self.modelBuilder.doVar("Rdy_mumu_ss[1.0,0.0,10.0]");
        self.modelBuilder.doVar("Rdy_ee_ss[1.0,0.0,10.0]");
      
        self.modelBuilder.factory_('expr::Rpl("(1.+@0)",Afb)')
        self.modelBuilder.factory_('expr::Rmn("(1.-@0)",Afb)')


 
 
 
 
    def getYieldScale(self,bin,process):
        if 'pl' in process: return "Rpl"
        if 'mn' in process : return "Rmn"
        if 'dy' in process and 'ee_ss' in bin: return "Rdy_ee_ss"
        if 'dy' in process and 'mumu_ss' in bin: return "Rdy_mumu_ss"
        else:
            return 1

class DY_AFB_dilu(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.75,0.75]");
        self.modelBuilder.doSet("POI","Afb")

        # ss templates
        self.modelBuilder.doVar("Dilu_ratio[1.0,0.0,10.0]");
        self.modelBuilder.doVar("Rdy_mumu_ss[1.0,0.0,10.0]");
        self.modelBuilder.doVar("Rdy_ee_ss[1.0,0.0,10.0]");
      
        self.modelBuilder.factory_('expr::Rpl("0.5*(1.+@0*@1)",Afb, Dilu_ratio)')
        self.modelBuilder.factory_('expr::Rmn("0.5*(1.-@0*@1)",Afb, Dilu_ratio)')


 
 
 
 
    def getYieldScale(self,bin,process):
        if 'pl' in process: return "Rpl"
        if 'mn' in process : return "Rmn"
        if 'dy' in process and 'ee_ss' in bin: return "Rdy_ee_ss"
        if 'dy' in process and 'mumu_ss' in bin: return "Rdy_mumu_ss"
        else:
            return 1

class Samesign(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Rdy[1.,0.0,10.0]");
        self.modelBuilder.doVar("Rqcd[1,0.0,10.0]");
        self.modelBuilder.doSet("POI","Rdy,Rqcd")
      
 
 
 
    def getYieldScale(self,bin,process):
        if 'dy' in process : return "Rdy"
        if 'qcd' in process : return "Rqcd"
        else:
            return 1

class EMu(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Rdy[1.,0.0,10.0]");
        self.modelBuilder.doVar("Rbk[1.,0.0,10.0]");
        self.modelBuilder.doVar("Rqcd_emu[1,0.0,10.0]");
        self.modelBuilder.doSet("POI","Rbk,Rdy,Rqcd_emu")
      
 
 
 
    def getYieldScale(self,bin,process):
        if 'bk' in process : return "Rbk"
        if 'dy' in process : return "Rdy"
        if 'qcd' in process : return "Rqcd_emu"
        else:
            return 1
samesign = Samesign()
emu = EMu()
dy_AFB = DY_AFB() 
dy_AFB_LO = DY_AFB_LO() 
dy_AFB_dilu = DY_AFB_dilu() 
