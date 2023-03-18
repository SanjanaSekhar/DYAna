from HiggsAnalysis.CombinedLimit.PhysicsModel import *

 
#This version has YLQ squared as a fundamental parameter, YLQ computed from that
class LQ_YLQ_SQ(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("A4[1.6, 0, 2.0]");
        self.modelBuilder.doVar("A0[0.05, 0.0, 2.0]");
        self.modelBuilder.doVar("yLQ2[0.001, -1.5, 1.5]");
        self.modelBuilder.doSet("POI","yLQ2")
        self.modelBuilder.doVar('expr::yLQ("((TMath::Abs(@0))**(0.5))",yLQ2)')
        #self.modelBuilder.doSet("POI","yLQ")

      
        self.modelBuilder.doVar('expr::Afb("3.0*@0/8.0",A4)');
        self.modelBuilder.factory_('expr::Alph("2.0*@0/(2.0-@0)",A0)')
        self.modelBuilder.factory_('expr::Norm("3.0/4.0/(2.0+@0)",Alph)')
        self.modelBuilder.factory_('expr::RAlph("@0*@1",Alph,Norm)')
        self.modelBuilder.factory_('expr::Rpl("(@0+@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::Rmn("(@0-@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::yLQ4("(@0*@0)",yLQ2)')
	self.modelBuilder.factory_('expr::minusyLQ2("(-1.0*@0)",yLQ2)')
 
 
 
 
    def getYieldScale(self,bin,process):
        if 'LQpure' in process: return "yLQ4"
        elif 'LQint_u' in process or 'LQint_c' in process: return "yLQ2"
	elif 'LQint_d_vec' in process or 'LQint_s_vec' in process: return "minusyLQ2"
	elif 'LQint_d' in process or 'LQint_s' in process: return "yLQ2"
        elif 'alpha' in process: return "RAlph"
        elif 'pl' in process: return "Rpl"
        elif 'mn' in process : return "Rmn"
        else:
            #print("Didnt find process %s bin %s in specifications \n" % (process, bin))
            return 1


 
#This version has YLQ as a fundamental parameter
class LQ_YLQ(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("A4[1.6, -2.0, 2.0]");
        self.modelBuilder.doVar("A0[0.05, -1.0, 1.0]");
        self.modelBuilder.doVar("yLQ[0.001, -4.0, 4.0]");
        self.modelBuilder.doSet("POI","yLQ")

      
        self.modelBuilder.doVar('expr::Afb("3.0*@0/8.0",A4)');
        self.modelBuilder.factory_('expr::Alph("2.0*@0/(2.0-@0)",A0)')
        self.modelBuilder.factory_('expr::Norm("3.0/4.0/(2.0+@0)",Alph)')
        self.modelBuilder.factory_('expr::RAlph("@0*@1",Alph,Norm)')
        self.modelBuilder.factory_('expr::Rpl("(@0+@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::Rmn("(@0-@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::yLQ2("(@0*@0)",yLQ)')
        self.modelBuilder.factory_('expr::yLQ4("(@0*@0)",yLQ2)')

 
 
    def getYieldScale(self,bin,process):
        if 'LQpure' in process: return "yLQ4"
        elif 'LQint' in process: return "yLQ2"
        elif 'alpha' in process: return "RAlph"
        elif 'pl' in process: return "Rpl"
        elif 'mn' in process : return "Rmn"
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
 
lq_ylq = LQ_YLQ()
lq_ylq_sq = LQ_YLQ_SQ()
