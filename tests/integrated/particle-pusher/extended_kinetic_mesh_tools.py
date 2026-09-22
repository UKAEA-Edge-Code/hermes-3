from boututils.run_wrapper import shell, launch_safe
from netCDF4 import Dataset
import numpy as np


def kinetic_mesh_path():
    return "extended_kinetic_mesh_data/BOUT.dmp.vantage.0.kinetic.msh"


def extended_kinetic_mesh_test_input(
    bout_grid_file,
    msh_file,
    dt,
    nsteps,
    iz_rate,
    rec_rate,
    remove_threshold,
    merge_threshold,
):

    input_file_string = f"""
    nout = 1
    timestep = 1

    [mesh]
    file="{bout_grid_file}"
    extrapolate_y=false
    extrapolate_x=false

    [dmplex]
    test_dmplex_cell_volumes = true
    test_dmplex_cell_centres = false
    use_external_msh = true
    msh_file = "{msh_file}"

    [solver]
    type = pvode

    [hermes]
    components = (d+, e, vantage)
    Nnorm = 1e19
    normalise_metric=false

    [d+]
    type = evolve_density
    AA = 1
    charge = 1

    [Nd+]
    function = 1

    [e]
    type = quasineutral
    AA = 1/1836
    charge = -1

    [vantage]
    dt = {dt}
    nsteps = {nsteps}
    test_mass_conservation = true
    initial_neutral_pressure = 1
    initial_neutral_temperature = 1
    npart_per_cell = 20
    background_ion_density = 1e19
    background_ion_temperature = 50
    background_ion_Vy = 1
    remove_threshold = {remove_threshold}
    merge_threshold = {merge_threshold}
    iz_rate_override = {iz_rate}
    rec_rate_override = {rec_rate}
    """
    return input_file_string


def generate_BOUT_grid_data(base_grid_dir, kinetic_nc_file_path, verbose=True):
    # make the directory for the basic hermes-3 run which
    # makes a slab grid
    cmd = f"mkdir {base_grid_dir}"
    shell(cmd)
    # this input file grid resolutions here
    # cannot be modified without changing
    # the mesh variables
    #   map_RZ_to_itriangle_0
    #   map_RZ_to_itriangle_1
    #   vertices
    #   tri_cell_vertices
    input_file_string = """
    nout = 0
    timestep = 1

    [mesh]
    J = 1

    nx = 8  # X grid size
    ny = 4  # Y grid size

    dx = 1.0/(nx-4)  # X mesh spacing
    dy = 2*pi/ny  # Y mesh spacing
    dz = 1  # Unity in toroidal direction

    Rxy = x
    Rxy_corners = x - 0.5*dx
    Rxy_lower_right_corners = x + 0.5*dx
    Rxy_upper_left_corners = x - 0.5*dx
    Rxy_upper_right_corners = x + 0.5*dx
    Zxy = y
    Zxy_corners = y - 0.5*dy
    Zxy_lower_right_corners = y - 0.5*dy
    Zxy_upper_left_corners = y + 0.5*dy
    Zxy_upper_right_corners = y + 0.5*dy

    [dmplex]
    use_cxx_ivertex=true
    test_dmplex_cell_volumes = true
    test_dmplex_cell_centres = true

    [solver]
    type = pvode

    [hermes]
    components = (d+, e, vantage)
    Nnorm = 1e19

    [d+]
    type = evolve_density
    AA = 1
    charge = 1

    [Nd+]
    function = 1

    [e]
    type = quasineutral
    AA = 1/1836
    charge = -1

    [vantage]
    nsteps = 0
    test_mass_conservation = true
    npart_per_cell = 1
    """
    file = f"{base_grid_dir}/BOUT.inp"
    with open(file, "w") as file:
        file.write(input_file_string)
    # Command to run
    cmd = f"./hermes-3 -d {base_grid_dir}"
    nproc = 1
    # Launch using MPI, with OMP_NUM_THREADS=1
    if verbose:
        print(f"execute: {cmd}")
    s, out = launch_safe(cmd, nproc=nproc, mthread=1, pipe=True)

    # copy the file to be used as a grid file, and insert the necessary
    # mesh variables that should be computed in preprocessing by the gridding/meshing workflow
    cmd = f"cp {base_grid_dir}/BOUT.dmp.vantage.0.nc {kinetic_nc_file_path}"
    if verbose:
        print(f"execute: {cmd}")
    shell(cmd)

    map_RZ_to_itriangle_0 = np.array(
        [
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [4.0, 6.0, 0.0, 2.0, 4.0, 6.0, 0.0, 2.0],
            [12.0, 14.0, 8.0, 10.0, 12.0, 14.0, 8.0, 10.0],
            [20.0, 22.0, 16.0, 18.0, 20.0, 22.0, 16.0, 18.0],
            [28.0, 30.0, 24.0, 26.0, 28.0, 30.0, 24.0, 26.0],
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        ]
    )
    map_RZ_to_itriangle_1 = np.array(
        [
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [5.0, 7.0, 1.0, 3.0, 5.0, 7.0, 1.0, 3.0],
            [13.0, 15.0, 9.0, 11.0, 13.0, 15.0, 9.0, 11.0],
            [21.0, 23.0, 17.0, 19.0, 21.0, 23.0, 17.0, 19.0],
            [29.0, 31.0, 25.0, 27.0, 29.0, 31.0, 25.0, 27.0],
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        ]
    )
    vertices = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.5707963267948966],
            [0.0, 3.141592653589793],
            [0.0, 4.71238898038469],
            [0.25, 0.0],
            [0.25, 1.5707963267948966],
            [0.25, 3.141592653589793],
            [0.25, 4.71238898038469],
            [0.5, 0.0],
            [0.5, 1.5707963267948966],
            [0.5, 3.141592653589793],
            [0.5, 4.71238898038469],
            [0.75, 0.0],
            [0.75, 1.5707963267948966],
            [0.75, 3.141592653589793],
            [0.75, 4.71238898038469],
            [1.0, 0.0],
            [1.0, 1.5707963267948966],
            [1.0, 3.141592653589793],
            [1.0, 4.71238898038469],
            [0.25, 6.283185307179586],
            [0.5, 6.283185307179586],
            [0.75, 6.283185307179586],
            [1.0, 6.283185307179586],
            [0.0, 6.283185307179586],
            [-1.0, 6.283185307179586],
            [-1.0, 0.0],
            [2.0, 0.0],
            [2.0, 6.283185307179586],
            [-0.4999999999986942, 6.283185307179586],
            [-1.0, 5.799863360474505],
            [-1.0, 5.316541413769544],
            [-1.0, 4.833219467065129],
            [-1.0, 4.349897520360367],
            [-1.0, 3.866575573654959],
            [-1.0, 3.383253626950145],
            [-1.0, 2.899931680244037],
            [-1.0, 2.41660973353672],
            [-1.0, 1.933287786829305],
            [-1.0, 1.449965840122101],
            [-1.0, 0.9666438934146182],
            [-1.0, 0.4833219467074823],
            [-0.500000000002059, 0.0],
            [1.5, 0.0],
            [2.0, 0.4833219467051673],
            [2.0, 0.9666438934102256],
            [2.0, 1.44996584011498],
            [2.0, 1.933287786819744],
            [2.0, 2.416609733524625],
            [2.0, 2.899931680229439],
            [2.0, 3.383253626935547],
            [2.0, 3.866575573642864],
            [2.0, 4.349897520350281],
            [2.0, 4.833219467057484],
            [2.0, 5.316541413764967],
            [2.0, 5.799863360472104],
            [1.5, 6.283185307179586],
            [-0.6071522086797074, 0.784761011285081],
            [-0.607152208680398, 5.498424295897261],
            [-0.5271436386483882, 2.235364003519466],
            [-0.5271436386490111, 4.047821303667518],
            [-0.6377166249249197, 1.29275877701848],
            [-0.6377166249255449, 4.990426530165874],
            [-0.581430915947327, 4.591558493712746],
            [-0.6217149109190194, 2.767018144897421],
            [-0.5814309159448909, 1.691626813475703],
            [-0.621714910919144, 3.516167162291901],
            [-0.6486859643676327, 3.14159265359466],
            [-0.6960731294749422, 5.95330761966039],
            [-0.3642360307921816, 5.790607309349078],
            [-0.696073129475413, 0.3298776875203251],
            [-0.3642360307925395, 0.4925779978312037],
            [1.607152208679471, 5.498424295894448],
            [1.60715220868047, 0.7847610112824089],
            [1.527143638648388, 4.047821303660118],
            [1.527143638649029, 2.235364003512242],
            [1.637716624924872, 4.990426530161094],
            [1.63771662492552, 1.292758777013975],
            [1.581430915947132, 1.691626813467362],
            [1.621714910919019, 3.516167162282163],
            [1.581430915944891, 4.59155849370388],
            [1.621714910919148, 2.767018144887718],
            [1.648685964367633, 3.141592653584932],
            [1.696073129475332, 0.3298776875191727],
            [1.364236030792596, 0.4925779978304361],
            [1.696073129474763, 5.953307619659157],
            [1.364236030791846, 5.790607309348226],
        ]
    )
    tri_cell_vertices = np.array(
        [
            [0, 4, 5],
            [0, 5, 1],
            [1, 5, 6],
            [1, 6, 2],
            [2, 6, 7],
            [2, 7, 3],
            [3, 7, 20],
            [3, 20, 24],
            [4, 8, 9],
            [4, 9, 5],
            [5, 9, 10],
            [5, 10, 6],
            [6, 10, 11],
            [6, 11, 7],
            [7, 11, 21],
            [7, 21, 20],
            [8, 12, 13],
            [8, 13, 9],
            [9, 13, 14],
            [9, 14, 10],
            [10, 14, 15],
            [10, 15, 11],
            [11, 15, 22],
            [11, 22, 21],
            [12, 16, 17],
            [12, 17, 13],
            [13, 17, 18],
            [13, 18, 14],
            [14, 18, 19],
            [14, 19, 15],
            [15, 19, 23],
            [15, 23, 22],
            [1, 59, 65],
            [60, 3, 63],
            [40, 41, 57],
            [30, 31, 58],
            [39, 40, 61],
            [31, 32, 62],
            [29, 25, 68],
            [26, 42, 70],
            [38, 39, 65],
            [32, 33, 63],
            [37, 38, 59],
            [33, 34, 60],
            [36, 37, 64],
            [35, 36, 67],
            [34, 35, 66],
            [41, 26, 70],
            [25, 30, 68],
            [57, 41, 70],
            [30, 58, 68],
            [40, 57, 61],
            [58, 31, 62],
            [2, 66, 67],
            [64, 2, 67],
            [39, 61, 65],
            [62, 32, 63],
            [59, 38, 65],
            [33, 60, 63],
            [3, 62, 63],
            [61, 1, 65],
            [37, 59, 64],
            [60, 34, 66],
            [36, 64, 67],
            [66, 35, 67],
            [57, 70, 71],
            [68, 58, 69],
            [42, 0, 71],
            [24, 29, 69],
            [70, 42, 71],
            [29, 68, 69],
            [0, 1, 71],
            [1, 2, 59],
            [1, 57, 71],
            [59, 2, 64],
            [2, 60, 66],
            [58, 3, 69],
            [57, 1, 61],
            [2, 3, 60],
            [3, 24, 69],
            [3, 58, 62],
            [19, 74, 80],
            [75, 17, 78],
            [54, 55, 72],
            [44, 45, 73],
            [53, 54, 76],
            [45, 46, 77],
            [43, 27, 83],
            [28, 56, 85],
            [52, 53, 80],
            [46, 47, 78],
            [51, 52, 74],
            [47, 48, 75],
            [50, 51, 79],
            [49, 50, 82],
            [48, 49, 81],
            [55, 28, 85],
            [27, 44, 83],
            [72, 55, 85],
            [44, 73, 83],
            [54, 72, 76],
            [73, 45, 77],
            [18, 81, 82],
            [79, 18, 82],
            [53, 76, 80],
            [77, 46, 78],
            [74, 52, 80],
            [47, 75, 78],
            [17, 77, 78],
            [76, 19, 80],
            [51, 74, 79],
            [75, 48, 81],
            [50, 79, 82],
            [81, 49, 82],
            [72, 85, 86],
            [83, 73, 84],
            [16, 43, 84],
            [56, 23, 86],
            [85, 56, 86],
            [43, 83, 84],
            [23, 19, 86],
            [19, 18, 74],
            [19, 72, 86],
            [74, 18, 79],
            [18, 75, 81],
            [73, 17, 84],
            [72, 19, 76],
            [18, 17, 75],
            [17, 16, 84],
            [17, 73, 77],
        ]
    )

    # open the kinetic nc file and append mesh data
    with Dataset(kinetic_nc_file_path, mode="a") as ncdataset:
        if verbose:
            print(f"append: map_RZ_to_itriangle_0 to {kinetic_nc_file_path}")
        ptr_map_RZ_to_itriangle_0 = ncdataset.createVariable(
            "map_RZ_to_itriangle_0", "f8", ("x", "y")
        )
        ptr_map_RZ_to_itriangle_0[:] = map_RZ_to_itriangle_0

        if verbose:
            print(f"append: map_RZ_to_itriangle_1 to {kinetic_nc_file_path}")
        ptr_map_RZ_to_itriangle_1 = ncdataset.createVariable(
            "map_RZ_to_itriangle_1", "f8", ("x", "y")
        )
        ptr_map_RZ_to_itriangle_1[:] = map_RZ_to_itriangle_1

        if verbose:
            print(f"append: vertices to {kinetic_nc_file_path}")
        nvertices, vertexdim = np.shape(vertices)
        ncdataset.createDimension("nvertices", nvertices)
        ncdataset.createDimension("vertexdim", vertexdim)
        ptr_vertices = ncdataset.createVariable(
            "vertices", "f8", ("nvertices", "vertexdim")
        )
        ptr_vertices[:] = vertices

        if verbose:
            print(f"append: tri_cell_vertices to {kinetic_nc_file_path}")
        ntriangle, tricorners = np.shape(tri_cell_vertices)
        ncdataset.createDimension("ntriangle", ntriangle)
        ncdataset.createDimension("tricorners", tricorners)
        ptr_tri_cell_vertices = ncdataset.createVariable(
            "tri_cell_vertices", "i4", ("ntriangle", "tricorners")
        )
        ptr_tri_cell_vertices[:] = tri_cell_vertices

    return None
