import os.path
from typing import List

from PySide6.QtCore import Slot
from PySide6.QtWidgets import QApplication, QMenu, QWidget
from PySide6.QtWidgets import QMenuBar
from hullformdir.shipstability import ShipStability
from optlib_scipy import ScipyOptimizationAlgorithm
try:
    # d3v imports
    from signals import Signals
    from commands import Command
    from iohandlers import IOHandler
    from core import geometry_manager as manager
    from core.geometry import Geometry
    from mainwin import MainFrame
    from glwin import GlWin
    from application import App
    # d3v-gsd
    from hullformdir.hullform import *
    from hullformdir.hullgeneratorform import HullGeneratorForm
    from hydrostatics_gui import HydrostaticsGUI
    from ship_stability_gui import ShipStabilityGUI
    from hullmoddir.occhullform import OCCHullform
    from hullform_command import HullFormCommand
    import hullmoddir.optocchullform as opthf
    # pygem
    from hullmoddir.pygemmenus import CreateFFDBox_menu, DeformFFDBox_menu
    from hullmoddir.pygemoptimisation import form_OptimizationProblem
    from michell_dir.michell import michell_resitance

except BaseException as error:
    print('An exception occurred: {}'.format(error))
except:
    print('Unknown exception occurred during signals connection')


def get_menu_by_name(mb: QMenuBar, name: str) -> QMenu:
    actions = mb.actions()
    for action in actions:
        if action.menu():
            if action.text() == name:
                return action
    return None


class HullmodCommand(Command):
    def __init__(self):
        super().__init__()
        self._app = QApplication.instance()
        self.app.registerIOHandler(HullmodImporter())
        self.hfcom: HullmodCommand = None
        for com in self.app.commands:
            if isinstance(com, HullFormCommand):
                self.hfcom = com
                break
        try:
            manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            manager.geometry_created.connect(self.onGeometryCreated)
            manager.geometry_removed.connect(self.onGeometryRemoved)
            manager.visible_geometry_changed.connect(self.onVisibleGeometryChanged)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')
        self.menuOCCForm = QMenu("&OCC Form")
        self.hfcom.menuMain.insertMenu(self.hfcom.menuHullGenForm.menuAction(), self.menuOCCForm)

        menuHullgenModify = self.menuOCCForm.addAction("&Modify form")
        menuHullgenModify.triggered.connect(self.onModifyHullgenForm)

        menuScale1D = self.menuOCCForm.addAction("&Scale 1-D")
        menuScale1D.triggered.connect(self.onScale1D)

        menuFullScale = self.menuOCCForm.addAction("&Full Scale")
        menuFullScale.triggered.connect(self.onFullScale)

        menuHullsections = self.menuOCCForm.addAction("&Cross section curves")
        menuHullsections.triggered.connect(self.onSectionCurves)

        menuHullextrudeSide = self.menuOCCForm.addAction("&Extrude hull side shell")
        menuHullextrudeSide.triggered.connect(self.onExtrudeSide)

        menuHullcloseTransom = self.menuOCCForm.addAction("&Close transom")
        menuHullcloseTransom.triggered.connect(self.onCloseTransom)

        menuHullcloseCover = self.menuOCCForm.addAction("&Close cover")
        menuHullcloseCover.triggered.connect(self.onCloseCover)

        menu_hullmod_opt = self.menuOCCForm.addAction("&Optimize form")
        menu_hullmod_opt.triggered.connect(self.on_hulmod_opt)

        menu_PrintHydroData = self.menuOCCForm.addAction("&Hydrostatic data")
        menu_PrintHydroData.triggered.connect(self.onPrintHydroData)
        #PyGem:

        menu_CreateFFDBox = self.menuOCCForm.addAction("&Create FFD Box")
        menu_CreateFFDBox.triggered.connect(self.onCreateFDDBox)

        menu_FFDDeform = self.menuOCCForm.addAction("&FDD Deform")
        menu_FFDDeform.triggered.connect(self.onFDDDeform)

        menu_CalcResist = self.menuOCCForm.addAction("&Calc Wave Resistance")
        menu_CalcResist.triggered.connect(self.on_calc_resistance)

        menu_CalcStab = self.menuOCCForm.addAction("&Calc Stab")
        menu_CalcStab.triggered.connect(self.on_calc_stab)

        menu_OptForm = self.menuOCCForm.addAction("&Optimise Form")
        menu_OptForm.triggered.connect(self.on_Opt_Form)
    #

    def on_Opt_Form(self):
        self.opt = form_OptimizationProblem(Hullform=self.hfcom.active_hull_form)
        #
        opt_ctrl = {}
        self.opt.opt_algorithm = ScipyOptimizationAlgorithm('SLSQP_mi=1000', 'SLSQP', opt_ctrl)
        res = self.opt.optimize()
        self.opt.print_output()
        self.hfcom.active_hull_form.visualise_surface()


        # hullform = self.hfcom.active_hull_form
        # calc = michell_resitance(self.hfcom.active_hull_form)
        # data1 = calc.wave_resistance(3.5)
        # print(data1)                                                                                    #(5.357551808216177, 0.0013889292068666674)
        # #                                                                                         # 0.17743882099421268 [-4.58088043e-04  2.34463456e-05  1.56304175e-01]
        # self.hfcom.active_hull_form.calc_stab()
        # print(self.hfcom.active_hull_form.displacement,self.hfcom.active_hull_form.displacementCG)
        #
        #
        # hullform.make_form_ffd_cage([0,0,0])
        # hullform.move_form_ffd_row(5,0)
        # hullform.move_form_ffd_row(4,0)
        #
        # deformed_surfaces = self.hfcom.active_hull_form.ffd_deform_surfaces()
        # self.hfcom.active_hull_form._surfaces = deformed_surfaces
        # self.hfcom.active_hull_form.regenerateHullHorm()
        #
        # calc = michell_resitance(self.hfcom.active_hull_form)
        # data1 = calc.wave_resistance(3.5)
        # print(data1)
        #
        # self.hfcom.active_hull_form.calc_stab()
        # print(self.hfcom.active_hull_form.displacement,self.hfcom.active_hull_form.displacementCG)
        #
        # self.hfcom.active_hull_form.visualise_surface()



        # self.hfcom.active_hull_form._surfaces =  cleaned_surfaces
        # calc = michell_resitance(self.hfcom.active_hull_form, True)
        # data1 = calc.wave_resistance(3.5)
        # print(data1)  # (5.357551808216177, 0.0013889292068666674)
        # self.hfcom.active_hull_form.calc_stab()
        # print(self.hfcom.active_hull_form.displacement,self.hfcom.active_hull_form.displacementCG)

        # deformed_surfaces = self.hfcom.active_hull_form.ffd_deform_surfaces(cleaned_surfaces)


        # self.resist = michell_resitance(hullform)
        # all_surfaces = hullform._surfaces
        # cleaned_surfaces = self.resist.mirrored_surfaces
        # one_side_surfaces = self.resist.one_side_surfaces


        # deformed_surfaces = self.hfcom.active_hull_form.ffd_deform_surfaces(cleaned_surfaces)
        # self.hfcom.active_hull_form.visualise_surface()
        # self.hfcom.active_hull_form.visualise_surface()
        # self.emit_geometries_rebuild()
        # self.opt = form_OptimizationProblem(Hullform = self.hfcom.active_hull_form)
        # print(self.opt.original_Cw, self.opt.original_displacement ,self.opt.original_displacementCGx)
        # print(res)

    def on_calc_stab(self):         #zrcalo i stvara novi mesh????
        self.hfcom.active_hull_form.calc_stab()
        print(self.hfcom.active_hull_form.displacement, self.hfcom.active_hull_form.displacementCG)
    #     # points = self.hfcom.active_hull_form.mesh.points()
    #     # self.hfcom.active_hull_form.visualise_surface()
    #     # z = points[:, 2]
    #     # print(min(z), max(z),self.hfcom.active_hull_form.T)
    #     # print(self.hfcom.active_hull_form.L, self.hfcom.active_hull_form.H, self.hfcom.active_hull_form.T)
    #     sscalc = ShipStability(self.hfcom.active_hull_form, (self.hfcom.active_hull_form.H-10))
    #     sscalc.wl.set_plane_point_z(self.hfcom.active_hull_form.T)
    #     self.displacement, self.displacementCG, new_fvs, new_pts = sscalc.calculate_displacement_and_displacementCG()
    #     # print('displacement, m3', displacement)
    #     # print('displacement, t', displacement)
    #     # print('displacement CG', displacementCG)

    def on_calc_resistance(self):
        calc = michell_resitance(self.hfcom.active_hull_form)
        # calc.set_new_michell()
        print(calc.wave_resistance(3.5))
        # calc.set_old_michell()
        # print(calc.wave_resistance(3.5))
        # calc.visualise()
        # print(calc.n_calls)
        # calc.test_intersection()


    def onCreateFDDBox(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.Box_creation_menu = CreateFFDBox_menu(self.hfcom.active_hull_form)
            self.Box_creation_menu.show_window()

    def onFDDDeform(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.Box_deform_menu = DeformFFDBox_menu(self.hfcom.active_hull_form)
            self.Box_deform_menu.show_window()
    def onPrintHydroData(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.printing_values()

    def onModifyHullgenForm(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.modify_form()

    def onScale1D(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.scale_1D()

    def onFullScale(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.full_scale()

    def onSectionCurves(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.cross_section_curves()

    def onExtrudeSide(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.extrude_side_shell()

    def onCloseTransom(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.close_transom()

    def onCloseCover(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.close_cover()

    def on_hulmod_opt(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            opthf.surface_through_curves_fit(self.hfcom.active_hull_form)

    @Slot()
    def onVisibleGeometryChanged(self, visible: List[Geometry], loaded: List[Geometry], selected: List[Geometry]):
        for g in visible:
            pass

    @Slot()
    def onSelectedGeometryChanged(self, visible: List[Geometry], loaded: List[Geometry], selected: List[Geometry]):
        pass

    @Slot()
    def onGeometryCreated(self, geometries: List[Geometry]):
        pass

    @Slot()
    def onGeometryRemoved(self, geometries: List[Geometry]):
        for g in geometries:
            if isinstance(g, OCCHullform):
                pass

    @property
    def app(self) -> App:
        return self._app

    @property
    def mainwin(self) -> MainFrame:
        return self.app.mainFrame

    @property
    def glwin(self) -> GlWin:
        return self.mainwin.glWin


class HullmodImporter(IOHandler):
    def __init__(self, force_import=False):
        super().__init__()

    def import_geometry(self, fileName):
        if len(fileName) < 1:
            return
        filename_no_ext, file_extension = os.path.splitext(os.path.basename(fileName))
        hf = None
        if file_extension in self.getImportFormats():
            hf = OCCHullform(fileName, filename_no_ext)
        if hf is not None:
            return hf

    def export_geometry(self, fileName, geometry2export):
        if isinstance(geometry2export, HullForm):
            geometry2export.exportGeometry(fileName)
        om.write_mesh(geometry2export.mesh, fileName)
        pass

    def getExportFormats(self):
        return ".hoc", ".igs"

    def getImportFormats(self):
        return ".hoc", ".igs"


def createCommand():
    return HullmodCommand()
