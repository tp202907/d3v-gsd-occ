from optbase import *
from optlib_scipy import ScipyOptimizationAlgorithm
from hullmoddir.michell_adaptor import michell_resitance

#varijable: pomicanje 4 tocke po x i y, treba li po z?
class form_analysis_model():
    def __init__(self, Hullform = None, optimise_aft = False, speed = 3.5): #speed in knots
        self.Hullform = Hullform
        self.speed = speed
        knots_to_add = [0, 0, 0]
        self.Hullform.make_form_ffd_cage(knots_to_add)
        self.optimise_aft = optimise_aft
        self.resistance_calc = michell_resitance(self.Hullform)
        self.Cw = self.resistance_calc.wave_resistance(self.speed)[1]
        self.original_Cw = copy(self.Cw)
        self.Hullform.calc_stab()
        self.original_displacement = copy(self.Hullform.displacement)
        self.original_displacementCGx = copy(self.Hullform.displacementCG[0])
        # print(self.qq)
        self.b_mo_x1 = 0.1
        self.b_mo_x2 = 0.1
        self.a_mo_x1 = -0.1
        self.a_mo_x2 = -0.1
        print("start:")
        print(self.get_b_mo_x1(),self.get_b_mo_x2(),self.get_a_mo_x1(),self.get_a_mo_x2())
        print(self.get_displacement())
        print(self.get_displacementCG_x())
        print(self.get_resistance())



    #bow functions:
    def set_b_mo_x1(self, mo_x):    #ffd displacement of bow row1, the one closest to 0,0
        self.b_mo_x1 = mo_x
        # for ind in self.brow1:
        #     self.ffd.array_mu_x[ind] = self.b_mo_x1

    def set_b_mo_x2(self, mo_x):    #ffd displacement of brow2
        self.b_mo_x2 = mo_x
        # for ind in self.brow2:
        #     self.ffd.array_mu_x[ind] = self.b_mo_x2

    def get_b_mo_x1(self):
        return self.b_mo_x1

    def get_b_mo_x2(self):
        return self.b_mo_x2

    #aft functions
    def set_a_mo_x1(self, mo_x):  # ffd displacement of brow1
        self.a_mo_x1 = mo_x
        # for ind in self.arow1:
        #     self.ffd.array_mu_x[ind] = self.a_mo_x1

    def set_a_mo_x2(self, mo_x):  # ffd displacement of brow1
        self.a_mo_x2 = mo_x
        # for ind in self.arow2:
        #     self.ffd.array_mu_x[ind] = self.a_mo_x2

    def get_a_mo_x1(self):
        return self.a_mo_x1

    def get_a_mo_x2(self):
        return self.a_mo_x2

    # dodaj za y i z ako treba

    def calc_resistance(self):
        self.Cw = self.resistance_calc.wave_resistance(self.speed)[1]

    def get_resistance(self):
        return self.Cw

    def get_displacement(self):
        return self.Hullform.displacement     #koje su granice volumena?

    def get_displacementCG_x(self):
        return self.Hullform.displacementCG[0]


    def analyze(self):
        self.Hullform._surfaces = self.Hullform.original_clean_surface    #reset to original surface
        self.Hullform.move_form_ffd_row(4, self.b_mo_x1)
        self.Hullform.move_form_ffd_row(5, self.b_mo_x2)


        if self.optimise_aft:
            self.Hullform.move_form_ffd_row(1, self.a_mo_x1)
            self.Hullform.move_form_ffd_row(2, self.a_mo_x2)


        deformed_surfaces = self.Hullform.ffd_deform_surfaces()
        self.Hullform._surfaces = deformed_surfaces
        # self.Hullform.visualise_surface()
        self.Hullform.regenerateHullHorm()
        self.resistance_calc = michell_resitance(self.Hullform)     #re initialise resist calc with original clean surfaces
        self.calc_resistance()
        self.Hullform.calc_stab()
        print(self.get_b_mo_x1(),self.get_b_mo_x2(),self.get_a_mo_x1(),self.get_a_mo_x2())
        print(self.get_displacement())
        print(self.get_displacementCG_x())
        print(self.get_resistance())
        return AnalysisResultType.OK

class form_OptimizationProblem(OptimizationProblem):
    def __init__(self,Hullform, name='', optimise_aft = True, speed = 3.5):
        if name == '':
            name = 'form'
        super().__init__(name)
        self.speed = speed
        am = form_analysis_model(Hullform = Hullform, optimise_aft = optimise_aft, speed = self.speed)
        self.deform_min = -0.1
        self.deform_max = 0.1

        b_mo_x1 = DesignVariable('b_mo_x1', CallbackGetSetConnector(am.get_b_mo_x1, am.set_b_mo_x1),self.deform_min,self.deform_max)
        b_mo_x2 = DesignVariable('b_mo_x2', CallbackGetSetConnector(am.get_b_mo_x2, am.set_b_mo_x2),self.deform_min,self.deform_max)
        self.add_design_variable(b_mo_x1)
        self.add_design_variable(b_mo_x2)

        #krma  varijable:
        if optimise_aft:
            a_mo_x1 = DesignVariable('a_mo_x1', CallbackGetSetConnector(am.get_a_mo_x1, am.set_a_mo_x1),self.deform_min,self.deform_max)
            a_mo_x2 = DesignVariable('a_mo_x2', CallbackGetSetConnector(am.get_a_mo_x2, am.set_a_mo_x2),self.deform_min,self.deform_max)
            self.add_design_variable(a_mo_x1)
            self.add_design_variable(a_mo_x2)


        self.add_objective(DesignObjective('Cw', CallbackGetConnector(am.get_resistance)))

        self.add_constraint(DesignConstraint('displacement_upper', CallbackGetConnector(am.get_displacement), am.original_displacement*1.02 , ConstrType.LT))
        self.add_constraint(DesignConstraint('displacement_lower', CallbackGetConnector(am.get_displacement), am.original_displacement*0.98, ConstrType.GT))

        # self.add_constraint(DesignConstraint('displacementCGx_upper', CallbackGetConnector(am.get_displacementCG_x), am.original_displacementCGx*1.4 , ConstrType.LT))
        # self.add_constraint(DesignConstraint('displacementCGx_lower', CallbackGetConnector(am.get_displacementCG_x), am.original_displacementCGx*0.6, ConstrType.GT))
        # self.add_constraint(DesignConstraint('displacementCGx_upper', CallbackGetConnector(am.get_displacementCG_x), 0.02 , ConstrType.LT))
        # self.add_constraint(DesignConstraint('displacementCGx_lower', CallbackGetConnector(am.get_displacementCG_x), -0.02, ConstrType.GT))

        self.add_analysis_executor(am)
