import OCC
import numpy as np
import openmesh as om
from pygem.cad import *
from copy import copy
from sys import maxsize
from OCC.Core.gp import gp_Trsf, gp_Lin
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.Tesselator import *
# from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Extend.TopologyUtils import *
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_BSplineSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from hullformdir.shipstability import ShipStability
from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Core.GeomFill import (
    GeomFill_BSplineCurves,
    GeomFill_StretchStyle,
    GeomFill_CoonsStyle,
    GeomFill_CurvedStyle,
)
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.Geom import Geom_BSplineCurve

from OCC.Display.SimpleGui import init_display

from OCC.Extend.ShapeFactory import point_list_to_TColgp_Array1OfPnt, make_face

def make_block(block_dims=np.array([20, 6, 3]), move_vector=np.array([0, 0, 0])):
    mesh = om.TriMesh()
    axes = []
    # stvara 2 tocke na svakoj osi
    for dim in block_dims:
        axes.append(np.linspace(0, dim, 2))

    block_corners = np.asarray(np.meshgrid(*axes)).T.reshape(8, 3)
    block_corners += move_vector

    # shrink block to try to avoid errors in boolean operations
    # block_corners[[4,5,6,7]] += np.array([0,0,-0.00005])			#z+ face is lowered
    # block_corners[[0,1,2,3]] += np.array([0,0,+0.00005])			#z- face is increased

    corner_vertices = []
    for corner in block_corners:
        corner_vertices.append(mesh.add_vertex(corner))

    # x+face
    mesh.add_face(corner_vertices[2], corner_vertices[3], corner_vertices[6])
    mesh.add_face(corner_vertices[3], corner_vertices[7], corner_vertices[6])

    # x-face
    mesh.add_face(corner_vertices[0], corner_vertices[4], corner_vertices[1])
    mesh.add_face(corner_vertices[1], corner_vertices[4], corner_vertices[5])

    # y+face
    mesh.add_face(corner_vertices[3], corner_vertices[1], corner_vertices[5])
    mesh.add_face(corner_vertices[3], corner_vertices[5], corner_vertices[7])

    # y-face
    mesh.add_face(corner_vertices[0], corner_vertices[2], corner_vertices[4])
    mesh.add_face(corner_vertices[2], corner_vertices[6], corner_vertices[4])

    # z+face
    mesh.add_face(corner_vertices[4], corner_vertices[6], corner_vertices[5])
    mesh.add_face(corner_vertices[6], corner_vertices[7], corner_vertices[5])

    # z-face
    mesh.add_face(corner_vertices[2], corner_vertices[0], corner_vertices[1])
    mesh.add_face(corner_vertices[2], corner_vertices[1], corner_vertices[3])

    return mesh
def soft_merge_meshes(meshes, vh_idx_to_sync_list=None):  # meshes je lista sa meshevima, vh_idx_to_sync_list sa lista isog lena ko meshes, svaka sadržava array sa vh_idx koji zelimo syncat
    points = np.empty((0, 3))
    merged_fvi = np.empty((0, 3))

    if vh_idx_to_sync_list is None:
        for mesh in meshes:
            mesh_fvi = mesh.face_vertex_indices()
            if mesh_fvi.size == 0:
                continue
            merged_fvi = np.append(merged_fvi, mesh_fvi + points.shape[0],
                                   axis=0)  # +points.shape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
            points = np.append(points, mesh.points(), axis=0)

        return om.TriMesh(points, merged_fvi)

    else:
        synced_vh_idx = []
        for i in range(len(meshes)):
            mesh = meshes[i]
            mesh_fvi = mesh.face_vertex_indices()
            merged_fvi = np.append(merged_fvi, mesh_fvi + points.shape[0],
                                   axis=0)  # +points.shape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
            synced_vh_idx.append(vh_idx_to_sync_list[i] + points.shape[0])
            points = np.append(points, mesh.points(), axis=0)

        return (om.TriMesh(points, merged_fvi), synced_vh_idx)

def hard_merge_meshes(meshes, vh_idx_to_sync_list=None):  # vh_idx_to_sync_list is list with numpy arrays
    if vh_idx_to_sync_list is None:
        merged_mesh = soft_merge_meshes(meshes)
        merged_mesh_fvi = merged_mesh.face_vertex_indices()
        merged_mesh_points = merged_mesh.points()

        bool = np.isclose(np.expand_dims(merged_mesh_points, 0), np.expand_dims(merged_mesh_points, 1)).all(
            -1)  # check for coincident points
        # clean diag and lower triangle matrix
        bool[np.diag_indices(bool.shape[0])] = False
        bool = np.triu(bool)
        ind = np.asarray(np.where(bool))
        # remove duplicates incase 3+ vh idx on same point
        data = np.unique(ind[1], return_index=True)  # [0] unique values, [1] their indices in orig array,
        ind = ind[:, data[1]]
        # delete duplicate points, replace duplicate vh_idx in fvi
        # duplicate vh_idx reduction:
        fvi_ind = np.where(np.expand_dims(merged_mesh_fvi, 0) == ind[1].reshape(-1, 1, 1))
        merged_mesh_fvi[fvi_ind[1:3]] = ind[0][fvi_ind[
            0]]  # slice fvi ind because [0] is indice of what vh_idx the fvi were compared to, 1,2 are actual indices of fvi to be replaced
        # syncing fvi after deleting duplicate points:
        vh_to_delete = np.unique(ind[1])
        vh_to_keep = np.delete(np.arange(merged_mesh_points.shape[0]), vh_to_delete, 0)
        merged_mesh_points = np.delete(merged_mesh_points, vh_to_delete, 0)
        fvi_ind = np.where(np.expand_dims(merged_mesh_fvi, 0) == vh_to_keep.reshape(-1, 1, 1))
        merged_mesh_fvi[fvi_ind[1:3]] = fvi_ind[
            0]  # slice fvi ind because [0] is indice of what vh_idx the fvi were compared to, 1,2 are actual indices of fvi to be replaced

        return om.TriMesh(merged_mesh_points, merged_mesh_fvi)

    else:
        data = soft_merge_meshes(meshes, vh_idx_to_sync_list)
        merged_mesh = data[0]
        vh_idx_to_sync_list = data[1]
        merged_mesh_fvi = merged_mesh.face_vertex_indices()
        merged_mesh_points = merged_mesh.points()

        bool = np.isclose(np.expand_dims(merged_mesh_points, 0), np.expand_dims(merged_mesh_points, 1)).all(-1)
        # clean diag and lower triangle matrix
        bool[np.diag_indices(bool.shape[0])] = False
        bool = np.triu(bool)
        ind = np.asarray(np.where(bool))
        # remove duplicates incase 3+ vh idx on same point
        data = np.unique(ind[1], return_index=True)  # [0] unique values, [1] their indices in orig array,
        ind = ind[:, data[1]]  # ind[0] new vh idx for duplicates, ind[1] old vh idx for duplicates
        # delete duplicate points, replace duplicate vh_idx in fvi
        # duplicate vh_idx reduction:
        fvi_ind = np.where(np.expand_dims(merged_mesh_fvi, 0) == ind[1].reshape(-1, 1, 1))
        merged_mesh_fvi[fvi_ind[1:3]] = ind[0][fvi_ind[
            0]]  # slice fvi ind because [0] is indice of what vh_idx the fvi were compared to, 1,2 are actual indices of fvi to be replaced
        # syncing fvi afrer deleting duplicate points:
        vh_to_delete = np.unique(ind[1])
        vh_to_keep = np.delete(np.arange(merged_mesh_points.shape[0]), vh_to_delete, 0)
        merged_mesh_points = np.delete(merged_mesh_points, vh_to_delete, 0)
        fvi_ind = np.where(np.expand_dims(merged_mesh_fvi, 0) == vh_to_keep.reshape(-1, 1, 1))
        merged_mesh_fvi[fvi_ind[1:3]] = fvi_ind[
            0]  # slice fvi ind because [0] is indice of what vh_idx the fvi were compared to, 1,2 are actual indices of fvi to be replaced

        # sync vh idx:
        synced_vh_idx_list = []
        for vh_idx_to_sync in vh_idx_to_sync_list:
            data = np.intersect1d(vh_idx_to_sync, ind[1], return_indices=True)
            vh_idx_to_sync[data[1]] = ind[0][data[2]]

            # syncing vi after deleting duplicate points:
            data = np.where(np.expand_dims(vh_idx_to_sync, 0) == np.expand_dims(vh_to_keep, 1))
            vh_idx_to_sync[data[1:3]] = data[0]

            synced_vh_idx_list.append(vh_idx_to_sync)
        return (om.TriMesh(merged_mesh_points, merged_mesh_fvi), synced_vh_idx_list)


def move_mesh(mesh, move_vector):
    return om.TriMesh(mesh.points() + move_vector, mesh.face_vertex_indices())

def is_mesh_closed(mesh):  # for use in dbb.py
    ffi = mesh.face_face_indices()
    if ffi.shape[1] <= 2:  # if all faces are only connected to 2 or less other faces, mesh i not closed by default
        return False
    if (ffi == -1).any():
        return False
    else:
        return True

class ffd_maker():
    def __init__(self):
        self._ffd_volumes = {}
        self._current_volume_id = 0

    #makes
    def make_form_ffd_cage(self, knots_to_add = [0,0,0]):
        self._ffd_volumes = {}  #reset any ffd volumes
        self._current_volume_id = 0

        box5_dims = np.array([self.L*1.01, self.B*1.01, self.H*1.01])
        box5_center = np.array([self.ship_start + self.L/2, 0, self.H / 2])
        self.make_ffd_volume(box5_center, box5_dims, n_control_points= [7,2,2], knots_to_add=knots_to_add)

        # print(self._ffd_volumes)

    def move_form_ffd_row(self,x_ind, mu_x,): #dx is absolute len by which to modify boxes #x_ind is index of row to move 0 is aft, 6 iis bow. mu_x is ffd deformation in x direction.
        row_indices = [[x_ind,0,0],[x_ind,1,0],[x_ind,0,1],[x_ind,1,1]]

        for ind in row_indices:
            self.move_ffd_pole(0, ind, [mu_x,0,0])
            current_ffd_volume = self._ffd_volumes[0]
            current_ffd_volume.array_mu_x[ind[0], ind[1], ind[2]] = mu_x


    def make_ffd_volume(self, box_center, box_length, n_control_points = [2,2,2], box_rotation_angles = [0,0,0], knots_to_add = [0,0,0]):
        ffd_volume = FFD(n_control_points, knots_to_add[0],knots_to_add[1],knots_to_add[2])
        # ffd_volume = FFD(n_control_points)
        ffd_volume.box_length = box_length
        ffd_volume.box_origin = (box_center - box_length/2)
        ffd_volume.rot_angle = box_rotation_angles
        # _ffd_volume.control_points()
        self._ffd_volumes[self._current_volume_id] = ffd_volume
        self._current_volume_id += 1

    def move_ffd_pole(self, ffd_id, pole_id:list, move_vector):
        current_ffd_volume = self._ffd_volumes[ffd_id]
        # print(dir(current_ffd_volume))
        current_ffd_volume.array_mu_x[pole_id[0], pole_id[1], pole_id[2]] = move_vector[0]
        current_ffd_volume.array_mu_y[pole_id[0], pole_id[1], pole_id[2]] = move_vector[1]
        current_ffd_volume.array_mu_z[pole_id[0], pole_id[1], pole_id[2]] = move_vector[2]
        # current_ffd_volume.array_mu_x[[*pole_id,]] = move_vector[0]
        # current_ffd_volume.array_mu_y[[*pole_id,]] = move_vector[1]
        # current_ffd_volume.array_mu_z[[*pole_id,]] = move_vector[2]

    def make_ffd_box_mesh(self):    #moza ih ne radi na tocno pravom mjestu pogledaj odakle je origin make_boxa
        pole_size = self.B/12
        pole = make_block(np.array([pole_size,]*3), np.array([0, 0, 0]))
        self.ffd_volume_mesh = []      # reset box mesh
        for ffd_id, ffd_volume in self._ffd_volumes.items():
            ffd_control_points = ffd_volume.control_points().tolist()

            for cpoint in ffd_control_points:
                pole_mesh = copy(pole)
                pole_mesh = move_mesh(pole_mesh, cpoint)
                self.ffd_volume_mesh.append(pole_mesh)


        # self.regenerateHullHorm()
        self.mesh = soft_merge_meshes(self.ffd_volume_mesh+[self.mesh,])
        # print(self.mesh)
        self.emit_geometries_rebuild()

class PyGemHullform(ffd_maker):
    def __init__(self):
        super().__init__()

    def position_form(self,x_offset = 0):    #positions the mesh and surfaces so that the x=0 is roughly at aft perpendicular and the min point is at z = 0, x offset is percentage L of how much to move the form in -x direction
        points = self.mesh.points()
        x = points[:,0]
        z = points[:,2]
        xmin = np.min(x)
        zmin = np.min(z)

        #move the sufaces:
        transformation = gp_Trsf()
        translation_vector = gp_Vec(-xmin - self.L*x_offset, 0, -zmin)
        transformation.SetTranslation(translation_vector)
        for i in range(len(self._surfaces)):
            self._surfaces[i] = BRepBuilderAPI_Transform(self._surfaces[i], transformation).Shape()
        #move mesh:
        self.mesh = move_mesh(self.mesh, np.array([-xmin - self.L*x_offset, 0, -zmin]))



    def matplotlib_visualise_mesh(self):
        import matplotlib.pyplot as plt
        ax = plt.figure().add_subplot(projection='3d')
        # plt.xlim(-self.L/2, self.L/2)
        # plt.ylim(-self.L/2, self.L/2)
        # plt.autoscale(False)
        mesh = self.mesh
        points = mesh.points()
        fvi = mesh.face_vertex_indices()
        tp = points[fvi]

        test_tp = np.array([[[0,1,2],[3,4,5],[6,7,8]],[[0,1,6],[3,4,5],[6,7,8]],[[0,1,2],[3,4,5],[6,7,8]]])

        # tp_bool = np.isclose(np.expand_dims(tp,0), np.expand_dims(tp,1)).all(-1).all(-1)
        # tp_bool = np.tril(tp_bool)
        # diag_ind = np.diag_indices(tp_bool.shape[0])
        # tp_bool[diag_ind] = False


        # print(tp_bool.any())      #nema conincidentnih tocaka


        normals = np.cross(tp[:, 1] - tp[:, 0], tp[:, 2] - tp[:, 0])
        normal_len = np.expand_dims(((normals ** 2).sum(-1)) ** 0.5, -1)  # normal lenghts
        normals = (normals / normal_len) * 0.2  # unit vectors
        # print(np.expand_dims(tp[:,0],0))
        tp_centroid = tp.sum(1)/3
        plot_tp = np.append(tp, np.expand_dims(tp[:,0],1),1)
        # print(plot_tp)


        #
        # tp_centroid = tp_centroid[0:10]
        # ax.scatter(tp_centroid[:,0],tp_centroid[:,1],tp_centroid[:,2], c = "red")

        # rc = np.array([self.L/2,self.B/2,self.H/2])
        rc = np.array([0,0,self.H/2])
        ax.scatter(rc[0],rc[1],rc[2], c = "green")
        x = tp_centroid - rc.reshape(1,-1)
        signed_normals = (x*normals).sum(-1)
        bool = signed_normals > 0

        normals = np.append(np.expand_dims(tp_centroid, 1), np.expand_dims(tp_centroid + normals, 1), 1)

        good_normals = copy(normals)
        good_normals[~bool] = np.array([[0,0,0],[0,0,0]])

        bad_normals = copy(normals)
        bad_normals[bool] = np.array([[0, 0, 0], [0, 0, 0]])

        # print(good_normals, bad_normals)
        #
        # ax.plot(plot_good_normals[:, 0], plot_good_normals[:, 1], plot_good_normals[:, 2], color="blue")
        # ax.plot(plot_bad_normals[:, 0], plot_bad_normals[:, 1], plot_bad_normals[:, 2], color="red")


        # ima losih faceva!!!!!!!!!!
        ns = 900
        ne = 1200

        # plot_tp = np.delete(plot_tp, ~bool, 0)
        # normals = np.delete(normals, ~bool, 0)
        # new_fvi = np.delete(fvi, ~bool, 0)
        # new_mesh = om.TriMesh(points, new_fvi)
        # self._mesh = new_mesh
        # print("new mesh !!!!!!!!!")
        # print(is_mesh_closed(mesh))

        # for t in list(plot_tp[ns:ne]):
        #     ax.plot(t[:,0], t[:,1], t[:,2], color = "blue", linewidth=0.25)
        #
        # for t in list(normals[ns:ne]):
        #     ax.plot(t[:,0], t[:,1], t[:,2], color = "green", linewidth=0.5)


        # for t in list(bad_normals[ns:ne]):
        #     ax.plot(t[:,0], t[:,1], t[:,2], color = "red", linewidth=1)
        #
        # ns = 95
        # ne = 100
        #
        for t in list(plot_tp):
            ax.plot(t[:,0], t[:,1], t[:,2], color = "blue", linewidth=0.25)
        #
        # for t in list(good_normals[ns:ne]):
        #     ax.plot(t[:,0], t[:,1], t[:,2], color = "green", linewidth=1)

        for t in list(bad_normals):
            ax.plot(t[:,0], t[:,1], t[:,2], color = "red", linewidth=1)



        # ax.plot(x, y, z, label='parametric curve')

        # xlim = plt.xlim()
        # max_xlim = max(map(abs, xlim))
        # plt.xlim((-max_xlim, max_xlim))
        # ylim = plt.ylim()
        # max_ylim = max(map(abs, ylim))
        # plt.ylim((-max_ylim, max_ylim))
        ax.legend()
        plt.show()


    def visualise_surface(self, surfaces = None):
        display, start_display, add_menu, add_function_to_menu = init_display()
        if surfaces is None:
            for i in self._surfaces:
                display.DisplayShape(i, update=True)
        else:
            for i in surfaces:
                display.DisplayShape(i, update=True)

        # for i in self.clean_surfaces:
        #     display.DisplayShape(i, update=True)
        # display.DisplayShape(self._clean_surface, update=True)
        # boxes = []
        for ffd_volume in self._ffd_volumes.values():
            cntpnts = ffd_volume.control_points(False)
            for cp in cntpnts:
                point = gp_Pnt(cp[0],cp[1],cp[2])
                display.DisplayShape(point, update=True)
                display.DisplayMessage(point, f"U")

        for ffd_volume in self._ffd_volumes.values():
            cntpnts = ffd_volume.control_points(True)
            for cp in cntpnts:
                point = gp_Pnt(cp[0],cp[1],cp[2])
                display.DisplayShape(point, update=True)
                display.DisplayMessage(point, f"D")
        # for box in boxes:
        #     display.DisplayShape(box, update=True)

        start_display()

    def ffd_deform_surfaces(self,surfaces = None):
        deformed_surfaces = []
        if surfaces is None:
            surfaces = self._surfaces

        for surface in surfaces:
            for ffd in self._ffd_volumes.values():
                deformed_surface =  ffd(surface)       #vraca compound?
                deformed_surface = self.regenerate_surface(deformed_surface)
            deformed_surfaces.append(deformed_surface)

        return deformed_surfaces


    def regenerate_surface(self, surface):
        if isinstance(surface, OCC.Core.TopoDS.TopoDS_Compound):      #nekada vraca compound a nekada surface...
            _compound = surface
            _expl = TopologyExplorer(_compound)
            for _face in _expl.faces():   #uvijek bi trebao biti samo jedan surface u compoundu
                _face_to_regenerate = _face
        else:
            _face_to_regenerate = surface

        surf = BRepAdaptor_Surface(_face_to_regenerate, True)
        bsrf = surf.BSpline()
        face = make_face(bsrf, 1e-6)
        return face

    def calc_stab(self):
        sscalc = ShipStability(self, self.H - 0.05)
        sscalc.wl.set_plane_point_z(self.T)
        self.displacement, self.displacementCG, new_fvs, new_pts = sscalc.calculate_displacement_and_displacementCG()
        # print('displacement, m3', displacement)
        # print('displacement, t', displacement)
        # print('displacement CG', displacementCG)

    def remove_form_deck_and_aft(self, remove_top_deck = True, remove_aft_surface = False):
        # makni zrcalo i palubu
        ray_origin = gp_Pnt(self.ship_start + self.L / 2, 0, self.H - 0.01)
        ray_dir_zpos = gp_Dir(0, 0, 1)
        ray_dir_xneg = gp_Dir(-1, 0, 0)
        ray_zpos = gp_Lin(ray_origin, ray_dir_zpos)
        ray_xneg = gp_Lin(ray_origin, ray_dir_xneg)
        surfs = []
        # clean surface:
        for face in self._surfaces:
            intersector = IntCurvesFace_ShapeIntersector()
            intersector.Load(face, 1e-6)
            if remove_top_deck:
                intersector.Perform(ray_zpos, 0, maxsize)  # check if surface is deck, if is, skip face
                # print(intersector.IsDone() and intersector.NbPnt() > 0)
                if intersector.IsDone() and intersector.NbPnt() > 0:
                    continue
            if remove_aft_surface:
                intersector.Perform(ray_xneg, 0, maxsize)  # check if surface is zrcalo, if is, skip face
                if intersector.IsDone() and intersector.NbPnt() > 0:
                    continue

            surfs.append(face)
            self._surfaces = surfs


if __name__ == "__main__":
    maker = ffd_maker()
    maker.make_ffd_volume(np.array([0,0,0]), np.array([1,1,1]))
    # print(maker.ffd_volumes[0].control_points())






















