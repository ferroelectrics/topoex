from argparse import ArgumentParser
import copy
import math
import random

import dolfin as dlf
from dolfin import (
    Constant,
    DirichletBC,
    Expression,
    FacetNormal,
    FiniteElement,
    Function,
    FunctionSpace,
    LinearVariationalProblem,
    Measure,
    Mesh,
    MeshFunction,
    MeshValueCollection,
    MixedElement,
    NonlinearProblem,
    PETScKrylovSolver,
    PETScOptions,
    PETScSNESSolver,
    SubDomain,
    TestFunction,
    TrialFunction,
    VectorElement,
    VectorFunctionSpace,
    UserExpression,
    XDMFFile,

    as_backend_type,
    assemble,
    curl,
    derivative,
    div,
    dot,
    errornorm,
    grad,
    inner,
    interpolate,
    lhs,
    near,
    norm,
    project,
    rhs,
    split
)

from equations import (
    FreeEnergy,
    FreeEnergyF
)

import ufl

dlf.parameters["form_compiler"]["optimize"]     = True
dlf.parameters["form_compiler"]["cpp_optimize"] = True
dlf.parameters["form_compiler"]["representation"] = "uflacs"
dlf.parameters["form_compiler"]["quadrature_degree"] = 2
dlf.parameters["form_compiler"]["cpp_optimize_flags"] = "-O2"
dlf.parameters["std_out_all_processes"] = False

random.seed()


class RandomDistribution(UserExpression):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def eval_cell(self, values, x, cell):
        values[0] = random.uniform(-1e-6, 1e-6)
        values[1] = random.uniform(-1e-6, 1e-6)
        values[2] = random.uniform(-1e-6, 1e-6)

    def value_shape(self):
        return (3,)


class NLProblem(NonlinearProblem):
    def __init__(self, a, L, bcs):
        NonlinearProblem.__init__(self)
        self.a = a
        self.L = L
        self.bcs = bcs

    def F(self, b, x):
        assemble(self.L, tensor=b)
        for bc in self.bcs:
            bc.apply(b, x)

    def J(self, A, x):
        assemble(self.a, tensor=A, keep_diagonal=True)
        for bc in self.bcs:
            bc.apply(A)
        A.ident_zeros()


def create_nlsolver():
    solver = dlf.PETScSNESSolver()

    PETScOptions.set('ksp_type', 'gmres')
    PETScOptions.set('pc_type', 'gasm')
    PETScOptions.set('snes_type', 'newtonls')
    PETScOptions.set('snes_linesearch_type', 'bt')
    solver.set_from_options()

    return solver


def create_kspsolver():
    solver = dlf.PETScKrylovSolver()

    prefix = 'linselec_'
    solver.set_options_prefix(prefix)

    PETScOptions.set(f'{prefix}ksp_type', 'gmres')
    PETScOptions.set(f'{prefix}pc_type', 'gasm')
    solver.set_from_options()

    return solver


def read_mesh(name3d, name2d):
    mesh = Mesh(dlf.MPI.comm_world)
    with XDMFFile(name3d) as infile:
        infile.read(mesh)

    subdomains_collection = MeshValueCollection('size_t', mesh, mesh.topology().dim())
    with XDMFFile(name3d) as infile:
        infile.read(subdomains_collection, "volumes")
    subdomains = MeshFunction('size_t', mesh, subdomains_collection)

    boundaries_collection = MeshValueCollection('size_t', mesh, mesh.topology().dim()-1)
    with XDMFFile(name2d) as infile:
        infile.read(boundaries_collection, "boundaries")
    boundaries = MeshFunction('size_t', mesh, boundaries_collection)

    return mesh, subdomains, boundaries


def create_measures(mesh, subdomains, boundaries):
    dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    dS = Measure('dS', domain=mesh, subdomain_data=boundaries)

    return dx, ds, dS


def read_msh_physical_names(msh_file):
    names = {}

    with open(msh_file) as msh:
        readline = lambda msh=msh: msh.readline().strip()
        while True:
            if readline() == '$PhysicalNames':
                nnames = int(readline())
                for i in range(nnames):
                    _, num, name = readline().split()
                    names[name.strip('"\'')] = int(num)

                break

    return names


def normal_rate(cur, prev):
    return (cur-prev) / cur


def dynamic_rate(cur, prev, dt):
    return abs((abs(cur) - abs(prev)) / abs(cur)) / dt


def write_to_xdmf(mesh, fname, f, filename, t, append=False):
    with XDMFFile(mesh.mpi_comm(), filename) as out:
        out.parameters['rewrite_function_mesh'] = False
        out.parameters['functions_share_mesh'] = True
        out.parameters['flush_output'] = True

        out.write_checkpoint(f, fname, t, append=append)


def read_from_xdmf(mesh, fname, f, filename, t=-1):
    with XDMFFile(mesh.mpi_comm(), filename) as input:
        input.read_checkpoint(f, fname, t)


def simulation(config):

    init_polar_name = config['init_polar_name']
    sim_name = config['sim_name']
    basename_msh = config['basename_msh']

    mesh, volumes, boundaries = read_mesh(f'{basename_msh}_volume.xdmf', f'{basename_msh}_boundary.xdmf')
    dx, ds, dS = create_measures(mesh, volumes, boundaries)
    entities = read_msh_physical_names(f'{basename_msh}.msh')

    SFE = FiniteElement('CG', mesh.ufl_cell(), 1)
    SFS = FunctionSpace(mesh, SFE)

    VFE = VectorElement('CG', mesh.ufl_cell(), 1)
    VFS = FunctionSpace(mesh, VFE)

    VFE_u = VectorElement('CG', mesh.ufl_cell(), 1)
    VFS_u = FunctionSpace(mesh, VFE_u)

    EpsFE = FiniteElement('DG', mesh.ufl_cell(), 0)
    EpsFS = FunctionSpace(mesh, EpsFE)

    ElFE = FiniteElement('DG', mesh.ufl_cell(), 0)
    ElFS = FunctionSpace(mesh, EpsFE)

    p = [Function(VFS) for _ in range(2)]
    pv = TestFunction(VFS)
    pd = TrialFunction(VFS)
    phi = Function(SFS)
    phiv = TestFunction(SFS)
    phid = TrialFunction(SFS)

    alphas = {  '1': config['a1'],
               '11': config['a11'],
               '12': config['a12'],
              '111': config['a111'],
              '112': config['a112'],
              '123': config['a123'] }

    gs = { '11':     config['g11'],
           '12':     config['g12'],
           '44':     config['g44'],
           '44prim': config['g44p'] }

    if not init_polar_name:
        p[1].interpolate(RandomDistribution(degree=2))
    else:
        read_from_xdmf(mesh, 'polar', p[2], init_polar_name)
    p[0].assign(p[1])

    dt = Constant(1e-3)

    Equations = FreeEnergyF(alphas, gs, p[0], pv, phi)
    TS_BE = dot(p[0]-p[1], pv) / dt
    F_BE = TS_BE*dx() + Equations*dx()
    J_BE = sum( derivative(F_BE, fun, dfun) for fun, dfun in zip(split(p[0]), split(pd)) )

    bcs_nls = []

    nlp_BE = NLProblem(J_BE, F_BE, bcs_nls)
    nl_solver = create_nlsolver()

    eps_0 = config['eps0']
    eps_i = config['epsi']

    a_ = a_Poisson(eps_0, eps_i, phid, phiv)*dx()
    L_ = L_Poisson(p[0], phiv)*dx()
    bcs = [ DirichletBC(SFS, Constant(0), boundaries, entities['top']) ]

    ksp_solver = create_kspsolver()

    A_P = assemble(a_)
    b_P = assemble(L_)
    for bc in bcs:
        bc.apply(A_P, b_P)
    ksp_solver.set_operator(A_P)

    Energy_C = FreeEnergy(alphas, gs, p[0], phi)

    current_iteration = 0
    T = 0.0
    current_energy = assemble(Energy_C*dx())
    previous_energy = current_energy
    delta_energy = 0.0

    save_prefix = f'{sim_name}'
    
    ksp_solver.solve(phi.vector(), b_P)
    niter, converged = nl_solver.solve(nlp_BE, p[0].vector())

    p[1].assign(p[0])
    b_P = assemble(L_)
    for bc in bcs:
        bc.apply(b_P)

    stop_solve = False
    while not stop_solve:
        ksp_solver.solve(phi.vector(), b_P)
        niter, converged = nl_solver.solve(nlp_BE, p[0].vector())

        T += float(dt)
        dlf.info(f'Time: {T}')
        dlf.info(f'Dt: {float(dt)}')

        current_energy = assemble( Energy_C*dx() )
        delta_energy = dynamic_rate(current_energy, previous_energy, float(dt))
        dlf.info(f'Full energy: {current_energy}')
        dlf.info(f'Delta energy: {delta_energy}')
    
        if delta_energy < 1e-6:
            write_to_xdmf(mesh, 'polar', p[0], f'{save_prefix}_polar_out.xdmf', T, append=False)
            stop_solve = True
            break

        current_iteration += 1
        previous_energy = current_energy
        p[1].assign(p[0])

        b_P = assemble(L_)
        for bc in bcs:
            bc.apply(b_P)


if __name__ == '__main__':

    parser = ArgumentParser(
                prog='Ferroelectric simulation',
                description='Simulation of topological excitations in ferrroic materials')

    parser.add_argument('-a1', '--alpha1', type=float, required=True)
    parser.add_argument('-a11', '--alpha11', type=float, required=True)
    parser.add_argument('-a12', '--alpha12', type=float, required=True)
    parser.add_argument('-a111', '--alpha111', type=float, required=True)
    parser.add_argument('-a112', '--alpha112', type=float, required=True)
    parser.add_argument('-a123', '--alpha123', type=float, required=True)

    parser.add_argument('-g11', '--grad11', type=float, required=True)
    parser.add_argument('-g12', '--grad12', type=float, required=True)
    parser.add_argument('-g44', '--grad44', type=float, required=True)
    parser.add_argument('-g44p', '--grad44prim', type=float, required=True)

    parser.add_argument('-e0', '--epsilon0', type=float, required=True)
    parser.add_argument('-ei', '--epsiloni', type=float, required=True)

    args = parser.parse_args()

    config = {
        'init_polar_name': None,
        'sim_name': 'initial_run',
        'basename_msh': 'box',
        'a1': args.alpha1,
        'a11': args.alpha11,
        'a12': args.alpha12,
        'a111': args.alpha111,
        'a112': args.alpha112,
        'a123': args.alpha123,
        'g11': args.grad11,
        'g12': args.grad12,
        'g44': args.grad44,
        'g44p': args.grad44prim,
        'eps0': args.epsilon0,
        'epsi': args.epsiloni,
    }

    simulation(config)

