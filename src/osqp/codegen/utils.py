"""
Utilities to generate embedded C code from OSQP sources
"""
# Compatibility with Python 2
from __future__ import print_function
from builtins import range

# Path of osqp module
import os.path
import osqp
files_to_generate_path = os.path.join(osqp.__path__[0],
                                      'codegen', 'files_to_generate')

# Timestamp
import datetime


def write_vec(f, vec, name, vec_type):
    """
    Write vector to file
    """
    if len(vec) > 0:

        f.write('%s %s[%d] = {\n' % (vec_type, name, len(vec)))

        # Write vector elements
        for i in range(len(vec)):
            if vec_type == 'c_float':
                f.write('(c_float)%.20f,\n' % vec[i])
            else:
                f.write('%i,\n' % vec[i])

        f.write('};\n')


def write_vec_extern(f, vec, name, vec_type):
    """
    Write vector prototype to file
    """
    if len(vec) > 0:
        f.write("extern %s %s[%d];\n" % (vec_type, name, len(vec)))


def write_mat(f, mat, name):
    """
    Write scipy sparse matrix in CSC form to file
    """
    write_vec(f, mat['p'], name + '_p', 'c_int')
    if len(mat['x']) > 0:
        write_vec(f, mat['i'], name + '_i', 'c_int')
        write_vec(f, mat['x'], name + '_x', 'c_float')

    f.write("csc %s = {" % name)
    f.write("%d, " % mat['nzmax'])
    f.write("%d, " % mat['m'])
    f.write("%d, " % mat['n'])
    f.write("%s_p, " % name)
    if len(mat['x']) > 0:
        f.write("%s_i, " % name)
        f.write("%s_x, " % name)
    else:
        f.write("0, 0, ")
    f.write("%d};\n" % mat['nz'])


def write_mat_extern(f, mat, name):
    """
    Write matrix prototype to file
    """
    f.write("extern csc %s;\n" % name)


def write_data_src(f, data, namespace):
    """
    Write data structure to file
    """
    f.write("// Define data structure\n")
    Pdata = f'{namespace}Pdata'
    Adata = f'{namespace}Adata'
    ldata = f'{namespace}ldata'
    qdata = f'{namespace}qdata'
    udata = f'{namespace}udata'
    data_name = f'{namespace}data'
    # Define matrix P
    write_mat(f, data['P'], Pdata)

    # Define matrix A
    write_mat(f, data['A'], Adata)

    # Define other data vectors
    write_vec(f, data['q'], qdata, 'c_float')
    write_vec(f, data['l'], ldata, 'c_float')
    write_vec(f, data['u'], udata, 'c_float')

    # Define data structure
    f.write(" ".join(("OSQPData", data_name, "=", "{")))
    f.write("%d, " % data['n'])
    f.write("%d, " % data['m'])
    f.write(f"&{Pdata}, &{Adata}, {qdata}, {ldata}, {udata}")
    f.write("};\n\n")


def write_data_inc(f, data, namespace):
    """
    Write data structure prototypes to file
    """
    f.write("// Data structure prototypes\n")

    # Define matrix P
    write_mat_extern(f, data['P'], f'{namespace}Pdata')

    # Define matrix A
    write_mat_extern(f, data['A'], f'{namespace}Adata')

    # Define other data vectors
    write_vec_extern(f, data['q'], f'{namespace}qdata', 'c_float')
    write_vec_extern(f, data['l'], f'{namespace}ldata', 'c_float')
    write_vec_extern(f, data['u'], f'{namespace}udata', 'c_float')

    # Define data structure
    f.write(f'extern OSQPData {namespace}data;\n\n')


def write_settings_src(f, settings, embedded_flag, namespace):
    """
    Write settings structure to file
    """
    f.write("// Define settings structure\n")
    f.write(" ".join(("OSQPSettings", f"{namespace}settings", "= {")))
    f.write("(c_float)%.20f, " % settings['rho'])
    f.write("(c_float)%.20f, " % settings['sigma'])
    f.write("%d, " % settings['scaling'])

    if embedded_flag != 1:
        f.write("%d, " % settings['adaptive_rho'])
        f.write("%d, " % settings['adaptive_rho_interval'])
        f.write("(c_float)%.20f, " % settings['adaptive_rho_tolerance'])

    f.write("%d, " % settings['max_iter'])
    f.write("(c_float)%.20f, " % settings['eps_abs'])
    f.write("(c_float)%.20f, " % settings['eps_rel'])
    f.write("(c_float)%.20f, " % settings['eps_prim_inf'])
    f.write("(c_float)%.20f, " % settings['eps_dual_inf'])
    f.write("(c_float)%.20f, " % settings['alpha'])
    f.write("(enum linsys_solver_type) LINSYS_SOLVER, ")

    f.write("%d, " % settings['scaled_termination'])
    f.write("%d, " % settings['check_termination'])
    f.write("%d, " % settings['warm_start'])

    f.write("};\n\n")


def write_settings_inc(f, settings, embedded_flag, namespace):
    """
    Write prototype for settings structure to file
    """
    f.write("// Settings structure prototype\n")
    f.write(f"extern OSQPSettings {namespace}settings;\n\n")


def write_scaling_src(f, scaling, namespace):
    """
    Write scaling structure to file
    """
    f.write("// Define scaling structure\n")
    if scaling is not None:
        write_vec(f, scaling['D'],    f'{namespace}Dscaling',    'c_float')
        write_vec(f, scaling['Dinv'], f'{namespace}Dinvscaling', 'c_float')
        write_vec(f, scaling['E'],    f'{namespace}Escaling',    'c_float')
        write_vec(f, scaling['Einv'], f'{namespace}Einvscaling', 'c_float')
        f.write(" ".join(("OSQPScaling", f"{namespace}scaling", "=", "{")))
        f.write("(c_float)%.20f, " % scaling['c'])
        f.write(f"{namespace}Dscaling, {namespace}Escaling, ")
        f.write("(c_float)%.20f, " % scaling['cinv'])
        f.write(" ".join((f"{namespace}Dinvscaling,", f"{namespace}Einvscaling", "};\n\n")))
    else:
        f.write(f"OSQPScaling {namespace}scaling;\n\n")


def write_scaling_inc(f, scaling, namespace):
    """
    Write prototypes for the scaling structure to file
    """
    f.write("// Scaling structure prototypes\n")

    if scaling is not None:
        write_vec_extern(f, scaling['D'],    f'{namespace}Dscaling',    'c_float')
        write_vec_extern(f, scaling['Dinv'], f'{namespace}Dinvscaling', 'c_float')
        write_vec_extern(f, scaling['E'],    f'{namespace}Escaling',    'c_float')
        write_vec_extern(f, scaling['Einv'], f'{namespace}Einvscaling', 'c_float')

    f.write(f"extern OSQPScaling {namespace}scaling;\n\n")


def write_linsys_solver_src(f, linsys_solver, embedded_flag, namespace):
    """
    Write linsys_solver structure to file
    """
    
    f.write("// Define linsys_solver structure\n")
    write_mat(f, linsys_solver['L'],            f'{namespace}linsys_solver_L')
    write_vec(f, linsys_solver['Dinv'],         f'{namespace}linsys_solver_Dinv',           'c_float')
    write_vec(f, linsys_solver['P'],            f'{namespace}linsys_solver_P',              'c_int')
    f.write(f"c_float {namespace}linsys_solver_bp[%d];\n"  % (len(linsys_solver['bp'])))
    f.write(f"c_float {namespace}linsys_solver_sol[%d];\n" % (len(linsys_solver['sol'])))
    write_vec(f, linsys_solver['rho_inv_vec'],  f'{namespace}linsys_solver_rho_inv_vec',    'c_float')

    if embedded_flag != 1:
        write_vec(f, linsys_solver['Pdiag_idx'], f'{namespace}linsys_solver_Pdiag_idx', 'c_int')
        write_mat(f, linsys_solver['KKT'],       f'{namespace}linsys_solver_KKT')
        write_vec(f, linsys_solver['PtoKKT'],    f'{namespace}linsys_solver_PtoKKT',    'c_int')
        write_vec(f, linsys_solver['AtoKKT'],    f'{namespace}linsys_solver_AtoKKT',    'c_int')
        write_vec(f, linsys_solver['rhotoKKT'],  f'{namespace}linsys_solver_rhotoKKT',  'c_int')
        write_vec(f, linsys_solver['D'],         f'{namespace}linsys_solver_D',         'QDLDL_float')
        write_vec(f, linsys_solver['etree'],     f'{namespace}linsys_solver_etree',     'QDLDL_int')
        write_vec(f, linsys_solver['Lnz'],       f'{namespace}linsys_solver_Lnz',       'QDLDL_int')
        f.write(f"QDLDL_int   {namespace}linsys_solver_iwork[{len(linsys_solver['iwork'])}];\n")
        f.write(f"QDLDL_bool  {namespace}linsys_solver_bwork[{len(linsys_solver['bwork'])}];\n")
        f.write(f"QDLDL_float {namespace}linsys_solver_fwork[{len(linsys_solver['fwork'])}];\n")

    f.write(f"qdldl_solver {namespace}linsys_solver = ")
    f.write(" ".join(("{QDLDL_SOLVER,", f"&solve_linsys_qdldl, ")))

    if embedded_flag != 1:
        f.write("&update_linsys_solver_matrices_qdldl, &update_linsys_solver_rho_vec_qdldl, ")

    f.write(f"&{namespace}linsys_solver_L, {namespace}linsys_solver_Dinv, "
            f"{namespace}linsys_solver_P, {namespace}linsys_solver_bp, "
            f"{namespace}linsys_solver_sol, {namespace}linsys_solver_rho_inv_vec, ")
    f.write("(c_float)%.20f, " % linsys_solver['sigma'])
    f.write("%d, " % linsys_solver['n'])
    f.write("%d, " % linsys_solver['m'])
    
    if embedded_flag != 1:
        if len(linsys_solver['Pdiag_idx']) > 0:
            linsys_solver_Pdiag_idx_string = f'{namespace}linsys_solver_Pdiag_idx'
            linsys_solver_PtoKKT_string = f'{namespace}linsys_solver_PtoKKT'
        else:
            linsys_solver_Pdiag_idx_string = '0'
            linsys_solver_PtoKKT_string = '0'
        if len(linsys_solver['AtoKKT']) > 0:
            linsys_solver_AtoKKT_string = f'{namespace}linsys_solver_AtoKKT'
        else:
            linsys_solver_AtoKKT_string = '0'
        f.write("%s, " % linsys_solver_Pdiag_idx_string)
        f.write("%d, " % linsys_solver['Pdiag_n'])
        f.write(f"&{namespace}linsys_solver_KKT, %s, %s, {namespace}linsys_solver_rhotoKKT, "
                % (linsys_solver_PtoKKT_string, linsys_solver_AtoKKT_string) +
                f"{namespace}linsys_solver_D, {namespace}linsys_solver_etree, {namespace}linsys_solver_Lnz, " +
                f"{namespace}linsys_solver_iwork, {namespace}linsys_solver_bwork, {namespace}linsys_solver_fwork, ")
    
    f.write("};\n\n")


def write_linsys_solver_inc(f, linsys_solver, embedded_flag, namespace):
    """
    Write prototypes for linsys_solver structure to file
    """
    f.write("// Prototypes for linsys_solver structure\n")
    write_mat_extern(f, linsys_solver['L'],    f'{namespace}linsys_solver_L')
    write_vec_extern(f, linsys_solver['Dinv'], f'{namespace}linsys_solver_Dinv', 'c_float')
    write_vec_extern(f, linsys_solver['P'],    f'{namespace}linsys_solver_P',    'c_int')
    f.write(f"extern c_float {namespace}linsys_solver_bp[%d];\n"  % len(linsys_solver['bp']))
    f.write(f"extern c_float {namespace}linsys_solver_sol[%d];\n" % len(linsys_solver['sol']))
    write_vec_extern(f, linsys_solver['rho_inv_vec'], f'{namespace}linsys_solver_rho_inv_vec', 'c_float')

    if embedded_flag != 1:
        write_vec_extern(f, linsys_solver['Pdiag_idx'], f'{namespace}linsys_solver_Pdiag_idx', 'c_int')
        write_mat_extern(f, linsys_solver['KKT'],       f'{namespace}linsys_solver_KKT')
        write_vec_extern(f, linsys_solver['PtoKKT'],    f'{namespace}linsys_solver_PtoKKT',    'c_int')
        write_vec_extern(f, linsys_solver['AtoKKT'],    f'{namespace}linsys_solver_AtoKKT',    'c_int')
        write_vec_extern(f, linsys_solver['rhotoKKT'],  f'{namespace}linsys_solver_rhotoKKT',  'c_int')
        write_vec_extern(f, linsys_solver['D'],         f'{namespace}linsys_solver_D',         'QDLDL_float')
        write_vec_extern(f, linsys_solver['etree'],     f'{namespace}linsys_solver_etree',     'QDLDL_int')
        write_vec_extern(f, linsys_solver['Lnz'],       f'{namespace}linsys_solver_Lnz',       'QDLDL_int')
        f.write(f"extern QDLDL_int   {namespace}linsys_solver_iwork[{len(linsys_solver['iwork'])}];\n")
        f.write(f"extern QDLDL_bool  {namespace}linsys_solver_bwork[{len(linsys_solver['iwork'])}];\n")
        f.write(f"extern QDLDL_float {namespace}linsys_solver_fwork[{len(linsys_solver['fwork'])}];\n")

    f.write(f"extern qdldl_solver {namespace}linsys_solver;\n\n")


def write_solution_src(f, data, namespace):
    """
    Preallocate solution vectors
    """
    f.write("// Define solution\n")
    f.write(f"c_float {namespace}xsolution[%d];\n" % data['n'])
    f.write(f"c_float {namespace}ysolution[%d];\n\n" % data['m'])
    f.write(" ".join(("OSQPSolution",  f"{namespace}solution",
                      "=", "{", f"{namespace}xsolution",
                      ",", f"{namespace}ysolution", "};\n\n")))


def write_solution_inc(f, data, namespace):
    """
    Prototypes for solution vectors
    """
    f.write("// Prototypes for solution\n")
    f.write(f"extern c_float {namespace}xsolution[{data['n']}];\n")
    f.write(f"extern c_float {namespace}ysolution[{data['m']}];\n\n" )
    f.write(f"extern OSQPSolution {namespace}solution;\n\n")


def write_info_src(f, namespace):
    """
    Preallocate info structure
    """
    f.write("// Define info\n")
    f.write(" ". join(("OSQPInfo",
                       f"{namespace}info",
                       "=",
                       '{0, "Unsolved", OSQP_UNSOLVED, 0.0, 0.0, 0.0};\n\n')))

def write_info_inc(f, namespace):
    """
    Prototype for info structure
    """
    f.write("// Prototype for info structure\n")
    f.write(f"extern OSQPInfo {namespace}info;\n\n")


def write_workspace_src(f, n, m, rho_vectors, embedded_flag, namespace):
    """
    Preallocate workspace structure and populate rho vectors
    """

    f.write("// Define workspace\n")

    write_vec(f, rho_vectors['rho_vec'],     f'{namespace}work_rho_vec',     'c_float')
    write_vec(f, rho_vectors['rho_inv_vec'], f'{namespace}work_rho_inv_vec', 'c_float')
    if embedded_flag != 1:
        write_vec(f, rho_vectors['constr_type'], f'{namespace}work_constr_type', 'c_int')

    f.write(f"c_float {namespace}work_x[{n}];\n")
    f.write(f"c_float {namespace}work_y[{m}];\n")
    f.write(f"c_float {namespace}work_z[{m}];\n")
    f.write(f"c_float {namespace}work_xz_tilde[%d];\n" % (m + n))
    f.write(f"c_float {namespace}work_x_prev[%d];\n" % n)
    f.write(f"c_float {namespace}work_z_prev[%d];\n" % m)
    f.write(f"c_float {namespace}work_Ax[%d];\n" % m)
    f.write(f"c_float {namespace}work_Px[%d];\n" % n)
    f.write(f"c_float {namespace}work_Aty[%d];\n" % n)
    f.write(f"c_float {namespace}work_delta_y[%d];\n" % m)
    f.write(f"c_float {namespace}work_Atdelta_y[%d];\n" % n)
    f.write(f"c_float {namespace}work_delta_x[%d];\n" % n)
    f.write(f"c_float {namespace}work_Pdelta_x[%d];\n" % n)
    f.write(f"c_float {namespace}work_Adelta_x[%d];\n" % m)
    f.write(f"c_float {namespace}work_D_temp[%d];\n" % n)
    f.write(f"c_float {namespace}work_D_temp_A[%d];\n" % n)
    f.write(f"c_float {namespace}work_E_temp[%d];\n\n" % m)

    f.write(" ".join(("OSQPWorkspace", f"{namespace}workspace", "= {\n")))
    f.write(f"&{namespace}data, (LinSysSolver *)&{namespace}linsys_solver,\n")
    f.write(f"{namespace}work_rho_vec, {namespace}work_rho_inv_vec,\n")
    if embedded_flag != 1:
        f.write(f"{namespace}work_constr_type,\n")

    f.write(f"{namespace}work_x, {namespace}work_y, {namespace}work_z, {namespace}work_xz_tilde,\n")
    f.write(f"{namespace}work_x_prev, {namespace}work_z_prev,\n")
    f.write(f"{namespace}work_Ax, {namespace}work_Px, {namespace}work_Aty,\n")
    f.write(f"{namespace}work_delta_y, {namespace}work_Atdelta_y,\n")
    f.write(f"{namespace}work_delta_x, {namespace}work_Pdelta_x, {namespace}work_Adelta_x,\n")
    f.write(f"{namespace}work_D_temp, {namespace}work_D_temp_A, {namespace}work_E_temp,\n")
    f.write(f"&{namespace}settings, &{namespace}scaling, &{namespace}solution, &{namespace}info" + "};\n\n")


def write_workspace_inc(f, n, m, rho_vectors, embedded_flag, namespace):
    """
    Prototypes for the workspace structure and rho_vectors
    """
    f.write("// Prototypes for the workspace\n")
    write_vec_extern(f, rho_vectors['rho_vec'],     f'{namespace}work_rho_vec',     'c_float')
    write_vec_extern(f, rho_vectors['rho_inv_vec'], f'{namespace}work_rho_inv_vec', 'c_float')
    if embedded_flag != 1:
        write_vec_extern(f, rho_vectors['constr_type'], f'{namespace}work_constr_type', 'c_int')

    f.write(f"extern c_float {namespace}work_x[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_y[%d];\n" % m)
    f.write(f"extern c_float {namespace}work_z[%d];\n" % m)
    f.write(f"extern c_float {namespace}work_xz_tilde[%d];\n" % (m + n))
    f.write(f"extern c_float {namespace}work_x_prev[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_z_prev[%d];\n" % m)
    f.write(f"extern c_float {namespace}work_Ax[%d];\n" % m)
    f.write(f"extern c_float {namespace}work_Px[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_Aty[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_delta_y[%d];\n" % m)
    f.write(f"extern c_float {namespace}work_Atdelta_y[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_delta_x[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_Pdelta_x[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_Adelta_x[%d];\n" % m)
    f.write(f"extern c_float {namespace}work_D_temp[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_D_temp_A[%d];\n" % n)
    f.write(f"extern c_float {namespace}work_E_temp[%d];\n\n" % m)

    f.write(f"extern OSQPWorkspace {namespace}workspace;\n\n")


def render_workspace(variables, hfname, cfname):
    """
    Print workspace dimensions
    """

    rho_vectors = variables['rho_vectors']
    data = variables['data']
    linsys_solver = variables['linsys_solver']
    scaling = variables['scaling']
    settings = variables['settings']
    embedded_flag = variables['embedded_flag']
    namespace: str = variables['namespace']

    n = data['n']
    m = data['m']

    # Open output file
    incFile = open(hfname, 'w')
    srcFile = open(cfname, 'w')

    # Add an include-guard statement
    fname = os.path.splitext(os.path.basename(hfname))[0]
    incGuard = fname.upper() + "_H"
    incFile.write("#ifndef %s\n" % incGuard)
    incFile.write("#define %s\n\n" % incGuard)

    # Print comment headers containing the generation time into the files
    now = datetime.datetime.now()
    daystr = now.strftime("%B %d, %Y")
    timestr = now.strftime("%H:%M:%S")
    incFile.write("/*\n")
    incFile.write(" * This file was autogenerated by OSQP-Python on %s at %s.\n" % (daystr, timestr))
    incFile.write(" * \n")
    incFile.write(" * This file contains the prototypes for all the workspace variables needed\n")
    incFile.write(" * by OSQP. The actual data is contained inside %sworkspace.c.\n" % namespace)
    incFile.write(" */\n\n")

    srcFile.write("/*\n")
    srcFile.write(" * This file was autogenerated by OSQP-Python on %s at %s.\n" % (daystr, timestr))
    srcFile.write(" * \n")
    srcFile.write(" * This file contains the workspace variables needed by OSQP.\n")
    srcFile.write(" */\n\n")

    # Include types, constants and linsys_solver header
    incFile.write("#include \"types.h\"\n")
    incFile.write("#include \"qdldl_interface.h\"\n\n")

    srcFile.write("#include \"types.h\"\n")
    srcFile.write("#include \"qdldl_interface.h\"\n\n")

    # Write data structure
    write_data_src(srcFile, data, namespace)
    write_data_inc(incFile, data, namespace)

    # Write settings structure
    write_settings_src(srcFile, settings, embedded_flag, namespace)
    write_settings_inc(incFile, settings, embedded_flag, namespace)

    # Write scaling structure
    write_scaling_src(srcFile, scaling, namespace)
    write_scaling_inc(incFile, scaling, namespace)

    # Write linsys_solver structure
    write_linsys_solver_src(srcFile, linsys_solver, embedded_flag, namespace)
    write_linsys_solver_inc(incFile, linsys_solver, embedded_flag, namespace)

    # Define empty solution structure
    write_solution_src(srcFile, data, namespace)
    write_solution_inc(incFile, data, namespace)

    # Define info structure
    write_info_src(srcFile, namespace)
    write_info_inc(incFile, namespace)

    # Define workspace structure
    write_workspace_src(srcFile, n, m, rho_vectors, embedded_flag, namespace)
    write_workspace_inc(incFile, n, m, rho_vectors, embedded_flag, namespace)

    # The endif for the include-guard
    incFile.write("#endif // ifndef %s\n" % incGuard)

    incFile.close()
    srcFile.close()


def render_setuppy(variables, output):
    """
    Render setup.py file
    """

    embedded_flag = variables['embedded_flag']
    python_ext_name = variables['python_ext_name']

    f = open(os.path.join(files_to_generate_path, 'setup.py'))
    filedata = f.read()
    f.close()

    filedata = filedata.replace("EMBEDDED_FLAG", str(embedded_flag))
    filedata = filedata.replace("PYTHON_EXT_NAME", str(python_ext_name))

    f = open(output, 'w')
    f.write(filedata)
    f.close()


def render_cmakelists(variables, output, namespace):
    """
    Render CMakeLists file
    """

    embedded_flag = variables['embedded_flag']

    f = open(os.path.join(files_to_generate_path, 'CMakeLists.txt'))
    filedata = f.read()
    f.close()

    filedata = filedata.replace(
        "EMBEDDED_FLAG", str(embedded_flag)
    ).replace("workspace.", f"{namespace}workspace.")

    f = open(output, 'w')
    f.write(filedata)
    f.close()


def render_emosqpmodule(variables, output, namespace):
    """
    Render emosqpmodule.c file
    """

    python_ext_name = variables['python_ext_name']

    f = open(os.path.join(files_to_generate_path, 'emosqpmodule.c'))
    filedata = f.read()
    f.close()

    filedata = filedata.replace(
        "PYTHON_EXT_NAME", str(python_ext_name)
    ).replace("workspace", f'{namespace}workspace')

    f = open(output, 'w')
    f.write(filedata)
    f.close()


def render_example(source, destination, namespace):

    with open(source, 'r') as src_fp:
        filedata = src_fp.read()

    filedata = filedata.replace(
        'workspace.h', f'&{namespace}workspace.h'
    ).replace('&workspace', f'&{namespace}workspace')

    with open(destination, 'w') as dest_fp:
        dest_fp.write(filedata)

