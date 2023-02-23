from OCC.Core.BRep import BRep_Tool_Pnt, BRep_Tool_Surface, BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform, BRepBuilderAPI_Transform, BRepBuilderAPI_NurbsConvert
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.GProp import GProp_GProps
from OCC.Core.GeomConvert import geomconvert_SurfaceToBSplineSurface
from OCC.Core.IGESControl import IGESControl_Reader, IGESControl_Writer, IGESControl_Controller_Init
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_VERTEX, TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.gp import gp_Vec, gp_Dir, gp_Pnt, gp_Ax2, gp_GTrsf, gp_Ax3, gp_Trsf, gp_Ax1
from OCC.Core.TopoDS import topods_Face
from OCC.Display.SimpleGui import init_display
from hullmoddir import occhullmanipulation as occhm

import scipy as sci
T = True
F = False
file_path_sur = "C:\\PycharmProjects\\MasterThesis\\d3v-gsd-occ\\examples\\ROPAX_105M16_2M.igs"
reader_sur = IGESControl_Reader()
# Reading the file, make sure you update path. Note in Windows slash needs to be used
reader_sur.ReadFile(file_path_sur)
# Prepareing the surface
reader_sur.TransferRoots()
shape_sur = reader_sur.Shape()
bbox = Bnd_Box()
bbox.SetGap(1e-5)
brepbndlib_Add(shape_sur, bbox, False)
xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
surf = BRep_Tool.Surface(shape_sur)
nurbs_surf = geomconvert_SurfaceToBSplineSurface(surf)


IGESControl_Controller_Init()
writer = IGESControl_Writer()
writer.AddGeom(nurbs_surf)
writer.Write("test_surf.igs")

if False:
    display, start_display, add_menu, add_function_to_menu = init_display()
    # display.DisplayShape(shape_sur, update=True, color='black')
    display.DisplayShape(nurbs_surf, update=True, color='red')
    start_display()

#
# cover_high = 10000
# edges = TopExp_Explorer(shape_sur, TopAbs_EDGE)
# z_ed = []
# while edges.More():
#     vertices = TopExp_Explorer(edges.Current(), TopAbs_VERTEX)
#     zmin = BRep_Tool_Pnt(vertices.Current()).Z()
#     vertices.Next()
#     zmax = BRep_Tool_Pnt(vertices.Current()).Z()
#     z_ed.append([zmin, zmax])
#     if z_ed[-1][0] <= zmin and z_ed[-1][1] <= zmax:
#         border = edges.Current()
#     edges.Next()
# pris = BRepPrimAPI_MakePrism(border, gp_Vec(0, 0, cover_high), False)
# pris.Build()
# pris_surf = BRep_Tool_Surface(pris.Shape())
#
# adap = BRepAdaptor_Surface(pris.Shape())
# print(pris.LastShape())
# sec = BRepAlgoAPI_Section(pris_surf, pris.LastShape())
# print(sec.Shape())
# geomconvert_SurfaceToBSplineSurface(pris_surf)
# shape_sur = BRepAlgoAPI_Fuse(shape_sur, pris).Shape()
# bbox = Bnd_Box()
# bbox.SetGap(1e-5)
# brepbndlib_Add(shape_sur, bbox, False)
# xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
#
# face_explorer = TopExp_Explorer(shape_sur, TopAbs_FACE)
#
# nurbs_faces = []
# while face_explorer.More():
#     nurbs_faces.append(face_explorer.Current())
#     face_explorer.Next()
# surfaces = [BRep_Tool_Surface(nurbs_face) for nurbs_face in nurbs_faces]
#
#
# # surf_adp = BRepAdaptor_Surface(nurbs_faces[0])
# nurbs_surfs = [geomconvert_SurfaceToBSplineSurface(surface) for surface in surfaces]
# print(nurbs_surfs)

