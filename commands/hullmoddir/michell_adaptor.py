from michell_dir.michell import michell
try:
    import numpy as np
except:
    pass
from scipy import integrate
import matplotlib.pyplot as plt
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Lin, gp_Pln
from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector
from sys import maxsize
from copy import copy
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Extend.TopologyUtils import *
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCC.Core.BRep import BRep_Tool, BRep_Tool_Curve
from OCC.Core.GeomConvert import geomconvert_CurveToBSplineCurve, GeomConvert_CompCurveToBSplineCurve
import time











class NoIntersection(Exception):
    pass



class michell_resitance(michell):
    def __init__(self, Hullform, remove_symmetry=True):  # speed in knots, later L,T,B in hullform
        self.Hullform = Hullform
        self.draft = self.Hullform.T
        self.n_calls = 0
        self.visualisation_points = []
        self.visualisation_sections = []
        self.visualisation_y = []
        self.density = 1000
        self.gravity = 9.81

        if self.Hullform is not None:
            self.length = self.Hullform.L
            self.beam = self.Hullform.B
            self.H = self.Hullform.H
            self.ship_start = self.Hullform.bbox.minCoord[0]  # aft
            self.ship_end = self.Hullform.bbox.maxCoord[0]  # bow
            self.lower_z_integration_limit = self.Hullform.bbox.minCoord[2]
            self.upper_z_integration_limit = self.lower_z_integration_limit + self.draft
            self.half_beam = self.interpolate_surface
            self.calc_wetted_area()
            self.clean_surfaces(remove_symmetry)
            self.prepare_hullform_approximation()


        else:  # rucno calc_wetted_area!!!!
            self.length = 0
            self.beam = 0
            self.H = 0
            self.ship_start = 0  # aft
            self.ship_end = 0  # bow
            self.lower_z_integration_limit = 0
            self.upper_z_integration_limit = 0
            self.half_beam = None



    def set_new_michell(self):
        self.length = self.Hullform.L
        self.beam = self.Hullform.B
        self.H = self.Hullform.H
        self.ship_start = self.Hullform.bbox.minCoord[0]  # aft
        self.ship_end = self.Hullform.bbox.maxCoord[0]  # bow
        self.lower_z_integration_limit = self.Hullform.bbox.minCoord[2]
        self.upper_z_integration_limit = self.lower_z_integration_limit + self.draft
        self.half_beam = self.interpolate_surface
        self.calc_wetted_area()
        # self.prepare_hullform_aproximation()
        xes = [1]
        zes = np.array(
            [0.00397997, 0.02049611, 0.04832857, 0.08446832, 0.125, 0.16553168, 0.20167143, 0.22950389, 0.24602003])
        print("new")
        for x in xes:
            print(self.half_beam(x, zes))
        # print(self.wave_resistance(3.5))

    def set_old_michell(self):
        self.draft = 0.25
        self.length = 4

        self.beam = 0.4

        # self.H = 0
        self.ship_start = 0  # aft
        self.ship_end = 4  # bow
        self.lower_z_integration_limit = 0
        self.upper_z_integration_limit = 0.25
        self.half_beam = self.wigly_half_beam
        self.calc_wetted_area()
        xes = [3]
        zes = np.array(
            [0.00397997, 0.02049611, 0.04832857, 0.08446832, 0.125, 0.16553168, 0.20167143, 0.22950389, 0.24602003])
        print("old")
        for x in xes:
            print(self.half_beam(x, zes))


    def prepare_hullform_approximation(self):
        # self.Fr = 0.51444*self.speed/(self.gravity*self.length)**0.5
        self._nurb_curves = {}      #a dictionary with key as x_cut, and value as [k, c, z] needed for interpolation
        # self.clean_surfaces()

        bow_n_curves = 20  #how many intervals are there in bow_range
        aft_n_curves = 20
        center_n_curves = 20
        bow_range = 0.15       # % of L that is considered part of bow
        aft_range = 0.15
        n_z_samples = 12    #how many points will it interpolate on the curve



        L = self.ship_end - self.ship_start
        aft_dx = L*aft_range/(aft_n_curves + 0.5)       #pocinnje na ship start - 0.5 aft dx kako bi vrh mogao dobro aproksimirati, vrh je uvijek vertikalan i vrijednosti su y = 0
        bow_dx = L*bow_range/(bow_n_curves + 0.5)
        center_dx = L*(1 - bow_range - aft_range)/(center_n_curves + 1)    #centar pocinje od L aft +dx, granica je dio afta

        aft_curve_xes = np.arange(aft_n_curves+2)*aft_dx + self.ship_start -0.5*aft_dx       #one extra curve for out of bounds 0, and for common one with center
        # print(aft_curve_xes[-1] + center_dx)
        center_curve_exes = np.arange(center_n_curves)*center_dx + aft_curve_xes[-1] + center_dx
        bow_curve_exes = np.arange(bow_n_curves+2)*bow_dx + center_curve_exes[-1] + center_dx

        self.xes = np.concatenate([aft_curve_xes, center_curve_exes, bow_curve_exes])

        t1 = time.time()
        for x_cut in self.xes.tolist():
            self.get_form_nurbs(x_cut, n_z_samples)
        t2 = time.time()
        print(f"Form nurbs generation done in: {t2-t1}")



    def visualise(self, points = np.array([])):
        from OCC.Display.SimpleGui import init_display
        display, start_display, add_menu, add_function_to_menu = init_display()
        display.DisplayShape(self._clean_surface_compound, update=True)
        for n in self.visualisation_sections:
            display.DisplayShape(n, update=True)

        start_display()

    def clean_surfaces(self, remove_symmetry = True):
        if remove_symmetry:
            self.remove_form_symmetry()
        else:                                       #make compound out of symmetrical form
            compound_builder = BRep_Builder()
            compound = TopoDS_Compound()
            compound_builder.MakeCompound(compound)
            for face in self.Hullform._surfaces:
                compound_builder.Add(compound, face)
            self._clean_surface_compound = compound



    def remove_form_symmetry(self):
        compound_builder = BRep_Builder()
        compound = TopoDS_Compound()
        compound_builder.MakeCompound(compound)
        # remove all surfaces from y- side of ship for faster calculation
        # find a point for every surface (u=0.5, v = 0.5) and check if that point is on right side of plane
        for face in self.Hullform._surfaces:
            bsrf = BRepAdaptor_Surface(face, True)  # .BSpline()  # .locateu, locatev

            ustart = bsrf.FirstUParameter()
            uend = bsrf.LastUParameter()
            vstart = bsrf.FirstVParameter()
            vend = bsrf.LastVParameter()
            midpoint = bsrf.Value(ustart + (uend - ustart) / 2, vstart + (vend - vstart) / 2)

            if midpoint.Y() < 0.0:
                continue


            compound_builder.Add(compound, face)
        self._clean_surface_compound = compound

    def get_form_nurbs(self, x_cut, n = 15, tolerance = 0.01):
        # cutting plane
        # print(x_cut)
        plane_point = gp_Pnt(x_cut, 0, self.H / 2)
        plane_normal = gp_Dir(1, 0, 0)
        plane = gp_Pln(plane_point, plane_normal)

        section = BRepAlgoAPI_Section(self._clean_surface_compound, plane, True)
        if section.SectionEdges().Size() == 0:  #if there are no edges (if plane was outside of form)
            self._nurb_curves[x_cut] = None     #if none, interpolation will return y = 0
        else:

            # display.DisplayShape(section.Shape(), update=True)

            explorer = TopExp_Explorer(section.Shape(), TopAbs_EDGE)
            composite_curve_builder = GeomConvert_CompCurveToBSplineCurve()

            # check for duplictes in wires: in order for them to be duplicates they need to both end and start points as same? for simplicity, calculate only for x coordinate
            # edges = np.empty(0, dtype=object)
            edges = np.empty(0, dtype=object)
            edge_start = np.empty((0, 3))
            edge_end = np.empty((0, 3))

            while explorer.More():
                edge = topods_Edge(explorer.Current())
                # print(edge)
                # display.DisplayShape(edge, update=True)
                # edge can be joined only if it is not degenerated (zero length)
                if BRep_Tool.Degenerated(edge):
                    print("Degenerate edge encountered when joining!")
                    explorer.Next()
                    continue

                # the edge must be converted to Nurbs edge
                nurbs_converter = BRepBuilderAPI_NurbsConvert(edge)
                nurbs_converter.Perform(edge)
                nurbs_edge = topods_Edge(nurbs_converter.Shape())

                # here we extract the underlying curve from the Nurbs edge
                nurbs_curve = BRep_Tool_Curve(nurbs_edge)[0]

                edges = np.append(edges, nurbs_curve)
                ustart = nurbs_curve.FirstParameter()
                uend = nurbs_curve.LastParameter()
                start = nurbs_curve.Value(ustart)
                end = nurbs_curve.Value(uend)
                npstart = np.array([[start.X(), start.Y(), start.Z()]])
                npend = np.array([[end.X(), end.Y(), end.Z()]])
                edge_start = np.append(edge_start, npstart, 0)
                edge_end = np.append(edge_end, npend, 0)
                explorer.Next()


            # clean wires from coincident instances:
            edges_arr = np.concatenate([np.expand_dims(edge_start, 1), np.expand_dims(edge_end, 1)], 1)
            # print(edges)

            # bool_arr = np.isclose(np.expand_dims(edges_arr, -4), np.expand_dims(edges_arr, -3))    # rtol=1e-01, atol=1e-01
            bool_arr = np.isclose(edges_arr.reshape(1, -1, 1, 2, 3), edges_arr.reshape(-1, 1, 2, 1, 3)).all(-1).any(-1).all(-1)  # rtol=1e-01, atol=1e-01
            bool_arr = np.tril(bool_arr)
            bool_arr[np.diag_indices(bool_arr.shape[0])] = False

            bad_edges_i = np.where(bool_arr)[0]  # dont know if always correct
            edges = np.delete(edges, bad_edges_i)

            for nurbs_curve in edges.tolist():

                bspline_curve = geomconvert_CurveToBSplineCurve(nurbs_curve)


                composite_curve_builder.Add(bspline_curve, tolerance)  # needs large tolerance?

            nurbs = composite_curve_builder.BSplineCurve()
            # self.visualisation_sections.append(nurbs)           #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # get parameters k and c for linear interpolation:

            u_end = nurbs.LastParameter()
            u_start = nurbs.FirstParameter()
            u_step = (u_end - u_start) / (n + 1)
            startpoint = nurbs.Value(u_start)
            endpoint = nurbs.Value(u_end)
            # nurbs_orientation = endpoint.Z() - startpoint.Z()

            # print(u_start,u_end, u_step)
            u_current = copy(u_start)
            z = []
            y = []

            for i in range(n + 2):  #+2 jer pocetnji i zadnji

                point = nurbs.Value(u_current)
                # nurbs.D0(u_current, point)
                u_current += u_step
                # self.visualisation_points.append(point)         # !!!!!!!!!!!!!!!!!!!!!!!!
                z.append(point.Z())
                y.append(point.Y())
                # tangents.append(tangent)

            z = np.asarray(z)
            y = np.asarray(y)
            z = z[::-1]
            y = y[::-1]

            y1 = y[:-1]
            y2 = y[1:]
            z1 = z[:-1]
            z2 = z[1:]

            k = (y2 - y1) / (z2 - z1)  # nagibi pravaca, k[i] vrijedi za tocku sa x-om izmedju x[i] i x[i+1]
            c = (y1 * z2 - z1 * y2) / (z2 - z1)


            self._nurb_curves[x_cut] = [k, c, z, y]



    def interpolate_surface(self,x_cut, z_dash):
        # print(x_cut, z_dash)
        bool = np.isclose(x_cut, self.xes)
        # print(bool)
        if bool.any():      #if x_cut is exactly on one of the curves, interpolate on that curve

            y = self.interpolate_curve(self.xes[bool][0], z_dash)
            self.visualisation_y.append([x_cut, y, z_dash])
            return y[::-1]


        bool = ~((x_cut <= self.xes[:-1]) | (x_cut >= self.xes[1:]))

        if bool.any():

            i = np.where(bool)[0]
            x1 = self.xes[i][0]        #const
            x2 = self.xes[i+1][0]
            y1 = self.interpolate_curve(x1, z_dash) #arr on same z
            y2 = self.interpolate_curve(x2, z_dash)

            k = (y2-y1)/(x2-x1)        #nagibi pravaca, k[i] vrijedi za tocku sa x-om izmedju x[i] i x[i+1]
            c = (y1*x2-x1*y2)/(x2-x1)

            yfinal = k * x_cut + c
            self.visualisation_y.append([x_cut, yfinal, z_dash])
            return yfinal[::-1]

        else:   #if x is not in any interval, return 0
            # print("3333333333333333")
            y = np.zeros(z_dash.shape)
            self.visualisation_y.append([x_cut,y,z_dash])
            return y[::-1]



    def interpolate_curve(self, x_cut, z_dash):
        data = self._nurb_curves[x_cut]
        if data is None:
            y = np.zeros(z_dash.shape)

            return y
        k = data[0]
        c = data[1]
        z = data[2]
        y_orig = data[3]

        z_dashr = np.expand_dims(z_dash, 0)
        zr = np.expand_dims(z, 1)
        y_dash_dummy = np.zeros(z_dash.shape, dtype=np.float64)
        on_point_ind = np.where(z_dash == zr)   #0 is z_ind, 1 = z_dash

        if z[0] - z[-1] > 0.0:      #check the orientation of z
            bool = ~((z_dashr <= zr[1:]) | (z_dashr >= zr[:-1]))

        elif z[0] - z[-1] < 0.0:
            bool = ~((z_dashr <= zr[:-1]) | (z_dashr >= zr[1:]))
        else:
            return y_dash_dummy



        bool[:, on_point_ind[1]] = False
        ind = np.where(bool)
        I_ind = ind[0]
        z_ind = ind[1]
        z_to_inter = z_dash[z_ind]

        k_ind = k[I_ind]
        c_ind = c[I_ind]

        y_dash = k_ind * z_to_inter + c_ind

        y_dash_dummy[on_point_ind[1]] = y_orig[on_point_ind[0]] #keep onpoint same
        y_dash_dummy[z_ind] = y_dash    #replace with interpolated values


        return y_dash_dummy



