import argparse
import meshio

if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(description='Convert gmsh-generated msh mesh to xdmf.')
    arg_parser.add_argument('--msh_mesh', type=str, help='name of msh mesh file')
    arg_parser.add_argument('--save_prefix', type=str, nargs='?', help='prefix to add to xdmf files')
    arg_parser.add_argument('--dim', type=str, help='dimesionality of mesh')

    args = arg_parser.parse_args()

    msh_name = args.msh_mesh

    prefix = msh_name[:-4] # drop .msh extension
    if args.save_prefix:
        prefix = args.save_prefix

    mesh_dim = args.dim
    volumetric, planar = {
        '3d': ('tetra', 'triangle'),
        '2d': ('triangle', 'line')
    }[mesh_dim]

    msh_mesh = meshio.read(msh_name)

    meshio.write(f"{prefix}_volume.xdmf", meshio.Mesh(
                                          points=msh_mesh.points, 
                                          cells=[(volumetric, msh_mesh.cells_dict[volumetric])],
                                          cell_data={"volumes": [msh_mesh.cell_data_dict["gmsh:physical"][volumetric],]},
                                          field_data=msh_mesh.field_data
                                            )
                )



    meshio.write(f"{prefix}_boundary.xdmf", meshio.Mesh(
                                        points=msh_mesh.points, 
                                        cells=[(planar, msh_mesh.cells_dict[planar])],
                                        cell_data={"boundaries": [msh_mesh.cell_data_dict["gmsh:physical"][planar],]},
                                        field_data=msh_mesh.field_data
                                            )
                )


