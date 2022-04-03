import sys
import importlib
from types import SimpleNamespace
import warnings
import numpy as np
import scipy.sparse as spa
import qdldl
from osqp import constant
from osqp import algebra_available, default_algebra
from osqp.interface import constant, _ALGEBRA_MODULES
import osqp.utils as utils
import osqp.codegen as cg
import pdb
import time


class OSQP:
    def __init__(self, *args, **kwargs):
        self.m = None
        self.n = None
        self.P = None
        self.q = None
        self.A = None
        self.l = None
        self.u = None

        self.algebra = kwargs.pop('algebra', default_algebra())
        if not algebra_available(self.algebra):
            raise RuntimeError(f'Algebra {self.algebra} not available')
        self.ext = importlib.import_module(_ALGEBRA_MODULES[self.algebra])

        self.settings = self.ext.OSQPSettings()
        self.ext.osqp_set_default_settings(self.settings)

        self._dtype = np.float32 if self.ext.OSQP_DFLOAT == 1 else np.float64
        self._itype = np.int64 if self.ext.OSQP_DLONG == 1 else np.int32

        # The following attributes are populated on setup()
        self._solver = None
        self._derivative_cache = {}

    def __str__(self):
        return f'OSQP with algebra={self.algebra}'

    @property
    def solver_type(self):
        return 'direct' if self.settings.linsys_solver == self.ext.osqp_linsys_solver_type.OSQP_DIRECT_SOLVER else 'indirect'

    @solver_type.setter
    def solver_type(self, value):
        assert value in ('direct', 'indirect')
        self.settings.linsys_solver = self.ext.osqp_linsys_solver_type.OSQP_DIRECT_SOLVER if value == 'direct' else self.ext.osqp_linsys_solver_type.OSQP_INDIRECT_SOLVER

    def constant(self, which):
        return constant(which, algebra=self.algebra)

    def update_settings(self, **kwargs):

        # Some setting names have changed. Support the old names for now, but warn the caller.
        renamed_settings = {'polish': 'polishing', 'warm_start': 'warm_starting'}
        for k, v in renamed_settings.items():
            if k in kwargs:
                warnings.warn(f'"{k}" is deprecated. Please use "{v}" instead.', DeprecationWarning)
                kwargs[v] = kwargs[k]
                del kwargs[k]

        new_settings = self.ext.OSQPSettings()
        for k in self.ext.OSQPSettings.__dict__:
            if not k.startswith('__'):
                if k in kwargs:
                    setattr(new_settings, k, kwargs[k])
                else:
                    setattr(new_settings, k, getattr(self.settings, k))

        if self._solver is not None:
            if 'rho' in kwargs:
                self._solver.update_rho(kwargs.pop('rho'))
            if kwargs:
                self._solver.update_settings(new_settings)
            self.settings = self._solver.get_settings()  # TODO: Why isn't this just an attribute?
        else:
            self.settings = new_settings

    def update(self, **kwargs):
        # TODO: sanity-check on types/dimensions

        q, l, u = kwargs.get('q'), kwargs.get('l'), kwargs.get('u')
        if l is not None:
            l = np.maximum(l, -constant('OSQP_INFTY'))
        if u is not None:
            u = np.minimum(u, constant('OSQP_INFTY'))

        if q is not None or l is not None or u is not None:
            self._solver.update_data_vec(q=q, l=l, u=u)
        if 'Px' in kwargs or 'Px_idx' in kwargs or 'Ax' in kwargs or 'Ax_idx' in kwargs:
            self._solver.update_data_mat(
                P_x=kwargs.get('Px'),
                P_i=kwargs.get('Px_idx'),
                A_x=kwargs.get('Ax'),
                A_i=kwargs.get('Ax_idx'),
            )

        if q is not None:
            self._derivative_cache['q'] = q
        if l is not None:
            self._derivative_cache['l'] = l
        if u is not None:
            self._derivative_cache['u'] = u

        for _var in ('P', 'A'):
            _varx = f'{_var}x'
            if kwargs.get(_varx) is not None:
                if kwargs.get(f'{_varx}_idx') is None:
                    self._derivative_cache[_var].data = kwargs[_varx]
                else:
                    self._derivative_cache[_var].data[kwargs[f'{_varx}_idx']] = kwargs[_varx]

        # delete results from self._derivative_cache to prohibit
        # taking the derivative of unsolved problems
        self._derivative_cache.pop('results', None)

    def setup(self, P, q, A, l, u, **settings):
        self.m = l.shape[0]
        self.n = q.shape[0]
        self.P = self.ext.CSC(spa.triu(P.astype(self._dtype), format='csc'))
        self.q = q.astype(self._dtype)
        self.A = self.ext.CSC(A.astype(self._dtype))
        self.l = l.astype(self._dtype)
        self.u = u.astype(self._dtype)

        self.update_settings(**settings)

        self._solver = self.ext.OSQPSolver(self.P, self.q, self.A, self.l, self.u, self.m, self.n, self.settings)
        self._derivative_cache.update({
            'P': P,
            'q': q,
            'A': A,
            'l': l,
            'u': u
        })

    def warm_start(self, x=None, y=None):
        # TODO: sanity checks on types/dimensions
        return self._solver.warm_start(x, y)

    def solve(self, raise_error=False):
        self._solver.solve()

        info = self._solver.info
        if info.status_val == constant('OSQP_NON_CVX', algebra=self.algebra):
            info.obj_val = np.nan
        # TODO: Handle primal/dual infeasibility

        if info.status_val != constant('OSQP_SOLVED') and raise_error:
            raise ValueError('Problem not solved!')

        # Create a Namespace of OSQPInfo keys and associated values
        _info = SimpleNamespace(**{k: getattr(info, k) for k in info.__class__.__dict__ if not k.startswith('__')})

        # TODO: The following structure is only to maintain backward compatibility, where x/y are attributes
        # directly inside the returned object on solve(). This should be simplified!
        results = SimpleNamespace(
            x=self._solver.solution.x,
            y=self._solver.solution.y,
            info=_info
        )

        self._derivative_cache['results'] = results
        return results

    def codegen(self, folder, project_type='', parameters='vectors', python_ext_name='emosqp', force_rewrite=False,
                FLOAT=False, LONG=True):
        return NotImplementedError

    def derivative_iterative_refinement(self, rhs, max_iter=1000, tol=1e-12):
        M = self._derivative_cache['M']

        # Prefactor
        solver = self._derivative_cache['solver']

        sol = solver.solve(rhs)

        for k in range(max_iter):
            delta_sol = solver.solve(rhs - M @ sol)
            sol = sol + delta_sol

            #  print("norm_iter_ref = %.4e\n" % np.linalg.norm(M @ sol - rhs))
            if np.linalg.norm(M @ sol - rhs) < tol:
                break

        if k == max_iter - 1:
            warnings.warn("max_iter iterative refinement reached.")
        # print('num_iters', k)
        
        return sol

    def adjoint_derivative(self, dx=None, dy_u=None, dy_l=None,
                           P_idx=None, A_idx=None, eps_iter_ref=1e-04):
        """
        Compute adjoint derivative after solve.
        """
        t0 = time.time()
        P, q = self._derivative_cache['P'], self._derivative_cache['q']
        A = self._derivative_cache['A']
        l, u = self._derivative_cache['l'], self._derivative_cache['u']

        try:
            results = self._derivative_cache['results']
        except KeyError:
            raise ValueError("Problem has not been solved. "
                             "You cannot take derivatives. "
                             "Please call the solve function.")

        if results.info.status != "solved":
            raise ValueError("Problem has not been solved to optimality. "
                             "You cannot take derivatives")

        m, n = A.shape
        x = results.x
        y = results.y
        y_u = np.maximum(y, 0)
        y_l = -np.minimum(y, 0)

        if A_idx is None:
            A_idx = A.nonzero()

        if P_idx is None:
            P_idx = P.nonzero()

        if dy_u is None:
            dy_u = np.zeros(m)
        if dy_l is None:
            dy_l = np.zeros(m)

        # identify equality constraints
        eq_indices = np.where(l == u)[0]
        ineq_indices = np.where(l < u)[0]
        num_eq = eq_indices.size
        A_ineq = A[ineq_indices, :]
        l_ineq = l[ineq_indices]
        u_ineq = u[ineq_indices]
        A_eq = A[eq_indices, :]
        b = u[eq_indices]

        # switch to Gx <= h form
        l_non_inf = np.where(l_ineq > -constant('OSQP_INFTY'))[0]
        u_non_inf = np.where(u_ineq < constant('OSQP_INFTY'))[0]
        # pdb.set_trace()
        num_ineq = l_non_inf.size + u_non_inf.size
        G = spa.bmat([
            [-A_ineq[l_non_inf, :]],
            [A_ineq[u_non_inf, :]],
        ])
        h = np.concatenate([-l_ineq[l_non_inf], u_ineq[u_non_inf]])

        nu = y[eq_indices]
        dnu = -dy_l[eq_indices] + dy_u[eq_indices]
        y_ineq = y[ineq_indices].copy()
        y_u_ineq = np.maximum(y_ineq, 0)
        y_l_ineq = -np.minimum(y_ineq, 0)
        lambd = np.concatenate([y_l_ineq[l_non_inf], y_u_ineq[u_non_inf]])
        # lambd = np.concatenate([-y_l[l_non_inf], y_u[u_non_inf]])
        dlambd = np.concatenate([-dy_l[l_non_inf], dy_u[u_non_inf]])

        # compute the derivative
        # Multiply second-third row by diag(y_u)^-1 and diag(y_l)^-1
        # to make the matrix symmetric
        inv_dia_y_u = spa.diags(np.reciprocal(y_u + 1e-20))
        inv_dia_y_l = spa.diags(np.reciprocal(y_l + 1e-20))
        inv_dia_lambda = spa.diags(np.reciprocal(lambd + 1e-20))
        dia_lambda = spa.diags(lambd)
        M = spa.bmat([
            [P, G.T],
            [G, spa.diags(G @ x - h) @ inv_dia_lambda]
        ], format='csc')
        slacks = G @ x - h
        slacks[slacks > -1e-4] = 0
        # M2 = spa.bmat([
        #     [P, G.T],
        #     [dia_lambda @ G, spa.diags(slacks)]
        # ], format='csc')
        M2 = spa.bmat([
            [P, G.T, A_eq.T],
            [dia_lambda @ G, spa.diags(slacks), None],
            [A_eq, None, None]
        ], format='csc')
        delta = spa.bmat([[eps_iter_ref * spa.eye(n), None],
                            [None, -eps_iter_ref * spa.eye(num_ineq)]],
                            format='csc')
        self._derivative_cache['M'] = M
        self._derivative_cache['solver'] = qdldl.Solver(M + delta)
        

        # rhs = - np.concatenate([dx, dlambd])
        rhs = - np.concatenate([dx, dlambd, dnu])
        
        # r_sol = self.derivative_iterative_refinement(rhs)
        # r_x, r_lambda_l, r_lambda_u, r_nu = np.split(r_sol, [n, n + l_non_inf.size, n + num_ineq])

        # try symmetrized
        B = spa.bmat([
            [spa.eye(n + num_ineq + num_eq), M2.T],
            [M2, None]
        ])
        # delta_B = spa.bmat([[eps_iter_ref * spa.eye(n + num_ineq), None],
        #                     [None, -eps_iter_ref * spa.eye(n + num_ineq)]],
        #                     format='csc')
        delta_B = spa.bmat([[eps_iter_ref * spa.eye(n + num_ineq + num_eq), None],
                            [None, -eps_iter_ref * spa.eye(n + num_ineq + num_eq)]],
                            format='csc')
        solver2 = qdldl.Solver(B + delta_B)
        self._derivative_cache['M'] = B
        self._derivative_cache['solver'] = solver2
        rhs_b = np.concatenate([rhs, np.zeros(n + num_ineq + num_eq)])
        
        r_sol_b = self.derivative_iterative_refinement(rhs_b)
        #### r_x_b, r_lambda_l_b, r_lambda_u_b, _ = np.split(r_sol_b, [n + num_ineq, n, n + l_non_inf.size, ])
        dual, primal = np.split(r_sol_b, [n + num_ineq + num_eq])

        # try lsqr
        # out = spa.linalg.lsqr(M2.T, rhs)
        # primal = out[0]

        r_x_b, r_lambda_l_b, r_lambda_u_b, r_nu = np.split(primal, [n, n + l_non_inf.size, n + num_ineq])
        r_x, r_lambda_l, r_lambda_u = r_x_b, r_lambda_l_b, r_lambda_u_b
        # pdb.set_trace()

        
        # # Make sure M matrix exists
        # if 'M' not in self._derivative_cache:
        #     # Multiply second-third row by diag(y_u)^-1 and diag(y_l)^-1
        #     # to make the matrix symmetric
        #     inv_dia_y_u = spa.diags(np.reciprocal(y_u + 1e-20))
        #     inv_dia_y_l = spa.diags(np.reciprocal(y_l + 1e-20))
        #     M = spa.bmat([
        #         [P, A.T, -A.T],
        #         [A, spa.diags(A @ x - u) @ inv_dia_y_u, None],
        #         [-A, None, spa.diags(l - A @ x) @ inv_dia_y_l]
        #     ], format='csc')
        #     delta = spa.bmat([[eps_iter_ref * spa.eye(n), None],
        #                       [None, -eps_iter_ref * spa.eye(2 * m)]],
        #                      format='csc')
        #     self._derivative_cache['M'] = M
        #     self._derivative_cache['solver'] = qdldl.Solver(M + delta)

        # rhs = - np.concatenate([dx, dy_u, dy_l])

        # r_sol = self.derivative_iterative_refinement(rhs)

        # r_x, r_yu, r_yl = np.split(r_sol, [n, n + m])



        # revert back to y form
        r_yu = np.zeros(m)
        r_yl = np.zeros(m)
        r_yu_ineq = np.zeros(m - num_eq)
        r_yl_ineq = np.zeros(m - num_eq)
        r_yu_ineq[u_non_inf] = r_lambda_u
        r_yl_ineq[l_non_inf] = -r_lambda_l

        # pdb.set_trace()
        r_yu[ineq_indices] = r_yu_ineq
        r_yl[ineq_indices] = r_yl_ineq

        # go from (r_nu, r_yu, r_yl) to (r_yu, r_yl)
        r_yu_eq = r_nu.copy() / nu
        r_yu_eq[nu < 0] = 0
        r_yu[eq_indices] = r_yu_eq

        r_yl_eq = -r_nu / (nu)
        r_yl_eq[nu >= 0] = 0
        r_yl[eq_indices] = r_yl_eq

        # Extract derivatives for the constraints
        rows, cols = A_idx
        # dA_vals = (y_u[rows] - y_l[rows]) * r_x[cols] + \
        #           (r_yu[rows] - r_yl[rows]) * x[cols]
        ryu = spa.diags(y_u) @ r_yu
        ryl = -spa.diags(y_l) @ r_yl
        dA_vals = (y_u[rows] - y_l[rows]) * r_x[cols] + \
                  (ryu[rows] - ryl[rows]) * x[cols]
        dA = spa.csc_matrix((dA_vals, (rows, cols)), shape=A.shape)
        # du = - r_yu
        # dl = r_yl
        du = - ryu
        dl = ryl

        # Extract derivatives for the cost (P, q)
        rows, cols = P_idx
        dP_vals = .5 * (r_x[rows] * x[cols] + r_x[cols] * x[rows])
        dP = spa.csc_matrix((dP_vals, P_idx), shape=P.shape)
        dq = r_x
        t1 = time.time()
        # print('derivative time', t1 - t0)
        pdb.set_trace()
        # print('dl', dl)
        return dP, dq, dA, dl, du
