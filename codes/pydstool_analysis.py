#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Limit cycle curve continuations using PyDSTool
# (needs: python3, mat73 [Python modules])
#
# NOTES:
# - The input parameters and initial states are obtained from a pre-computed
#   bifurcation diagram (e.g., after running 'bd_draw_alpha.m').
# - This is not a fully automatic continuation. Knowledge from the pre-computed
#   diagram is needed.

from itertools import product
from mat73 import loadmat
from pathlib import Path
import PyDSTool as dst


def save_lcc(PC, var, headers, out_file, file_sep=','):
    with open(out_file, 'w+') as f:
        var = ['_T'] + var
        headers = ['period'] + headers + ['stability', 'type']
        f.write(file_sep.join(headers))
        for s in PC:
            f.write('\n')
            f.write(file_sep.join([str(s[v]) for v in var]))
            info = s.labels
            f.write(file_sep + info['LC']['stab'])
            if 'LPC' in info:
                f.write(file_sep + 'LPC')
            elif 'PD' in info:
                f.write(file_sep + 'PD')
            elif 'NS' in info:
                f.write(file_sep + 'NS')
            elif 'RG' in info:
                f.write(file_sep + 'RG')
            elif 'BP' in info:
                f.write(file_sep + 'BP')
            else:
                f.write(file_sep + 'LC')
    return None


def cont_lcc(PC):
    """Continue limit cycle curves
    """
    if not 'BP' in PC.LocBifPoints:
        PC.LocBifPoints.append('BP')
    PC.forward()
    return None


def cont_epc(PC):
    """Continue Equilibrium Point Curve with increasing 'y_0' (i.e., 'E_Pyr')
    """
    m = PC.MaxNumPoints
    PC.MaxNumPoints = 2
    PC.forward()
    PC.MaxNumPoints = m
    if (PC.sol['y_0'][1] > PC.sol['y_0'][0]):
        PC.forward()
    else:
        PC.backward()
    return None


def setup_args_lcc(curve_name, curve_name_epc):
    PC_args = dst.args(name=curve_name, type='LC-C')
    PC_args.freepars = ['q']
    PC_args.initpoint = curve_name_epc + ':H1'
    PC_args.LocBifPoints = 'all'
    PC_args.MaxNumPoints = 10000
    PC_args.MaxStepSize = 10**(-1)
    PC_args.MinStepSize = 10**(-5)
    PC_args.NumCollocation = 4
    PC_args.NumIntervals = 24
    PC_args.SaveEigen = False
    PC_args.StepSize = 10**(-4)
    PC_args.SolutionMeasures = ['max', 'min']
    PC_args.UseAuto = True
    PC_args.verbosity = 0

    return PC_args


def setup_args_epc(curve_name):
    PC_args = dst.args(name=curve_name, type='EP-C')
    PC_args.freepars = ['q']
    PC_args.LocBifPoints = 'all'
    PC_args.MaxNumPoints = 140
    PC_args.MaxStepSize = 10**(-1)
    PC_args.MinStepSize = 10**(-3)
    PC_args.SaveEigen = False
    PC_args.StepSize = 10**(-2)
    PC_args.StopAtPoints = ['H']
    PC_args.verbosity = 0

    return PC_args


def setup_system(param):
    """Intended to be used for the diagram (q, .)
    """

    ds_args = dst.args(name='namm')
    ds_args.pars = {**{k: v[()] for (k, v) in param.items()}, **{'q': 0.0}}
    ds_args.fnspecs = {
        'S': (['x', 'theta'], 'nu_max / (1 + exp(r * (theta - x)))'),
        'dS': (['x', 'theta'], 'nu_max * r * exp(r * (theta - x)) / (1 + exp(r * (theta - x)))**2'),
        'Jacobian': (['t', 'y_0', 'y_1', 'y_2', 'y_3', 'y_4', 'y_5'], """[
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [-a * a, A * a * dS(y_1 - y_2, v_0 - v_Glu + v_GABA), -A * a * dS(y_1 - y_2, v_0 - v_Glu + v_GABA), -2 * a, 0.0, 0.0],
            [A * a * C_ExIn_to_Pyr * C_Pyr_to_ExIn * dS(C_Pyr_to_ExIn * y_0, v_0), C_Pyr_to_Pyr * A * a * dS(y_1 - y_2, v_0 - v_Glu + v_GABA) - a * a, -C_Pyr_to_Pyr * A * a * dS(y_1 - y_2, v_0 - v_Glu + v_GABA), 0.0, -2 * a, 0.0],
            [B * b * C_InIn_to_Pyr * C_Pyr_to_InIn * dS(C_Pyr_to_InIn * y_0, v_0 - mu_Glu_InIn_by_Pyr * v_Glu), 0.0, -b * b, 0.0, 0.0, -2 * b],
            ]"""),
    }
    ds_args.varspecs = {
        'y_0': 'y_3',
        'y_1': 'y_4',
        'y_2': 'y_5',
        'y_3': 'a * (A * S(y_1 - y_2, v_0 - v_Glu + v_GABA) - 2 * y_3 - a * y_0)',
        'y_4': 'a * (A * (q + C_Pyr_to_Pyr * S(y_1 - y_2, v_0 - v_Glu + v_GABA) + C_ExIn_to_Pyr * S(C_Pyr_to_ExIn * y_0, v_0)) - 2 * y_4 - a * y_1)',
        'y_5': 'b * (B * C_InIn_to_Pyr * S(C_Pyr_to_InIn * y_0, v_0 - mu_Glu_InIn_by_Pyr * v_Glu) - 2 * y_5 - b * y_2)',
    }
    ds_args.ics = {
        'y_0': 0.0,
        'y_1': 0.0,
        'y_2': 0.0,
        'y_3': 0.0,
        'y_4': 0.0,
        'y_5': 0.0,
    }
    ds_args.xdomain = {
        'y_0': [
            10**(-10),
            ds_args.pars['A'] * ds_args.pars['nu_max'] / ds_args.pars['a'] - 10**(-10)
            ]
    }
    var_python_matlab = [
        ('y_0', 'E_Pyr'),
        ('y_1', 'E_ExIn_and_Pyr'),
        ('y_2', 'E_InIn'),
        ('y_3', 'd_E_Pyr'),
        ('y_4', 'd_E_ExIn_and_Pyr'),
        ('y_5', 'd_E_InIn'),
    ]

    return (dst.Generator.Vode_ODEsystem(ds_args), var_python_matlab)


# Main
if __name__ == "__main__":
    # * Example: limit cycle continuations for the JR-diagram.
    # * The input file name is used below to name outputs.
    # * The link between [NAME].bd.mat and [NAME].bd.pydstool.mat is used below to
    #   name outputs.
    proj_dir = Path('./data+codes_pcbi2025_aliobaibk')
    if not proj_dir.is_dir():
        raise NotADirectoryError(
            'The "proj_dir" variable is not set to an existing directory.\n'
            'Ensure "proj_dir" points to the directory where the repository was cloned.\n'
            'For example, if you cloned the repository into '
            '"./data+codes_pcbi2025_aliobaibk",\n'
            'use: proj_dir = Path("./data+codes_pcbi2025_aliobaibk").'
    )

    in_bd_file = proj_dir / 'data' / 'bifurcation-diagrams' / 'alpha' / \
        'C_Pyr_to_Pyr-0.bd.pydstool.mat'
    BD = loadmat(in_bd_file, use_attrdict=True)['P']


    # * BD is a collection of codimension-1 diagrams with q as parameter.
    # * Different indices (diagrams) reflect different values of a parameter other
    #   than q (e.g., v_Glu).
    # * Loop the code below over these indices if needed (excluding the variable
    #   'hopf_indices_matlab' and adjusting 'args_lcc', which are diagram-specific).
    # * Indexing from 1 to facilitate MATLAB interface.
    bd_index_matlab = 301
    BD = BD[bd_index_matlab - 1][0]


    # * Not all Hopf points need to be continued.
    #   In the JR-diagram, consisting of: saddle-saddle-hopf-hopf-hopf (sshhh),
    #   only the first and last Hopf points (third and fifth bifurcation points) are
    #   relevant for continuation (e.g., check < BD['bif_label']['name'] >).
    # * Indexing from 1 to facilitate MATLAB interface.
    hopf_indices_matlab = [3, 5]


    # * Sets up continuations.
    (ode, var_python_matlab) = setup_system(BD['param'])
    name_epc = 'EQ1'
    args_epc = setup_args_epc(name_epc)
    name_lcc = 'LC1'
    args_lcc = setup_args_lcc(name_lcc, name_epc)


    # * Not all state variables are used or useful to save.
    var_valid = [(n_p, n_m) for (n_p, n_m) in var_python_matlab if n_m in BD]

    var = [n_p for (n_p, _) in var_valid]
    var = ['{}_{}'.format(v, m) for (v, m) in product(var,
        args_lcc.SolutionMeasures)]
    var = args_lcc.freepars + var

    headers = [n_m for (_, n_m) in var_valid]
    headers = ['{}_{}'.format(h, m) for (h, m) in product(headers,
        args_lcc.SolutionMeasures)]
    headers = args_lcc.freepars + headers


    # - Main
    out_dir = in_bd_file.parent / in_bd_file.stem.replace('.pydstool', '_lcc')
    out_dir.mkdir(parents=True, exist_ok=True)

    for h in hopf_indices_matlab:
        out_file = out_dir / '{}-{}.lcc'.format(int(BD['original'][()]), h)

        k = int(BD['bif_label']['k_sing'][h - 1][1][()] - 2)
        ode.set(ics={n_p: BD[n_m][k] for (n_p, n_m) in var_valid})
        ode.set(pars={'q': BD['q'][k]})

        PC = dst.ContClass(ode)

        PC.newCurve(args_epc)
        cont_epc(PC[name_epc])

        PC.newCurve(args_lcc)
        cont_lcc(PC[name_lcc])
        # run the last command < cont_lcc(PC[name_lcc]) > as many times as needed or ...
        # ... adjust < PC[name_lcc].MaxNumPoints >

        save_lcc(PC[name_lcc].sol, var, headers, out_file)
