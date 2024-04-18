import matplotlib.pyplot as plt
import matplotlib.colors
from intvalpy import lineqs
import sympy
import numpy as np
import os
import json
from scipy.optimize import linprog
from ipywidgets import interact, widgets
from functools import lru_cache
import time
#plt.rcParams['text.usetex'] = True

def coeff_str(coeff, sym):
    '''
    Represents the expression coeff * sym as a string where coeff is a number and sym is a symbol
    0 * x -> ''
    1 * x -> 'x'
    -1 * x -> '-x'
    n * x -> 'nx'
    -n * x -> '-nx'
    '''
    if coeff==0:
        return ''
    p = '' if coeff >= 0 else '-'
    c = abs(int(coeff)) if coeff%1==0 else abs(coeff)
    return p + (str(c)+'$\cdot$' if c != 1 else '') + sym

def coeffs_str(coeffs, syms):
    s = ''
    for i in range(len(coeffs)):
        if coeffs[i] > 0:
            s += '+'
        s += coeff_str(coeffs[i], syms[i])
    return s

def ineq_str(A, b, i):
    cx, cy =-A[i,0], -A[i,1]
    s = [t for t in [coeff_str(cx, 'x'), coeff_str(cy, 'y')] if t]
    return ' + '.join(s) + ' $\leq$ ' + str(int(-b[i])) if b[i]%1==0 else str(-b[i])

class MILP:
    def __init__(self, n=2):
        self.c = np.zeros((n))
        self.A = np.zeros((0,n))
        self.b = np.zeros((0))
        self.B = np.zeros((0,n))
        self.d = np.zeros((0))
        self.n, self.n1, self.m, self.p = n, n, 0, 0

    def set_objective(self, c):
        self.c = np.array(c, dtype=np.float64)
        assert(self.n == self.c.shape[0])

    def add_inequality_constraint(self, Ai, bi, insert_at_begin=False):
        self.m += 1
        self.A = np.vstack([self.A, Ai]) if not insert_at_begin else np.vstack([Ai, self.A])
        self.b = np.append(self.b, bi) if not insert_at_begin else np.append(bi, self.b)

    def check_dimensions(self):
        assert(self.c.shape[0] == self.n)
        assert(self.A.shape[0] == self.m and self.A.shape[1] == self.n)
        assert(self.n == 2)
        assert(self.p == 0)

    def invert_constraints(self):
        '''
        Multiplies the inequality constraints by -1 on both sides.
        Since Ax >= b <=> -Ax <= -b, this is used to either interpret the inequality constraints with a "smaller-equal" or "greater-equal"
        '''
        self.A = -self.A
        self.b = -self.b

    def load(self, file:str):
        data = json.load(open(file))
        self.A = np.array(data['A'], dtype=np.float64)
        self.b = np.array(data['b'], dtype=np.float64)
        self.c = np.array(data['c'], dtype=np.float64)
        self.B = np.array(data['B'], dtype=np.float64)
        self.d = np.array(data['d'], dtype=np.float64)
        self.n, self.n1, self.m, self.p = data['n'], data['n1'], data['m'], data['p']

    def save(self, file:str):
        data = {}
        data['A'] = self.A.tolist()
        data['b'] = self.b.tolist()
        data['B'] = self.B.tolist()
        data['d'] = self.d.tolist()
        data['c'] = self.c.tolist()
        data['n'], data['n1'], data['m'], data['p'] = self.n, self.n1, self.m, self.p
        with open(file, "w") as f:
            json.dump(data, f)

    def objective_str(self) -> str:
        return r'$\min_{x,y}\:$' + coeffs_str(self.c, ['x', 'y'])

    def __str__(self):
        return "c = {c}\nA =\n{A}\nb =\n{b}".format(c=self.c,A=self.A,b=self.b)

    def save_plot(self, title:str, sol=None, cuts=(np.array([[]]),np.array([])), mark_last_ineq=False, save_file:str=''):
        c,B,d = self.c,self.B,self.d
        A,b = (self.A,self.b) if cuts[1].shape[0]==0 else (np.vstack([self.A, cuts[0]]), np.append(self.b, cuts[1]))

        # Get Extents of Feasible Region
        vs = lineqs(self.A, self.b, show=False)
        xlim, ylim = [-0, 1.1*np.max(vs[:,0])], [-0, 1.1*np.max(vs[:,1])]   

        ratio = (xlim[1]-xlim[0])/(ylim[1]-ylim[0])
        fig = plt.figure(figsize=(10,10/ratio))
        assert(len(c) == 2)
        m = len(b)

        # Generate Meshgrid of Points
        x_vals = np.linspace(xlim[0], xlim[1], 400)
        y_vals = np.linspace(ylim[0], ylim[1], 400)
        X, Y = np.meshgrid(x_vals, y_vals)

        # Plot Objective Function
        Z = c[0]*X + c[1]*Y
        plt.imshow(Z, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], origin='lower', cmap='viridis', alpha=0.3)
        plt.contour(X, Y, Z, levels=20, cmap='viridis', alpha=0.1)

        # Plot Inequality Constraints
        for i in range(m):
            if A[i,1] != 0:
                plt.plot(x_vals, (b[i] - A[i,0]*x_vals) / A[i,1], 'r-' if (mark_last_ineq and i==m-1) else 'g-', lw=2, label=ineq_str(A,b,i))
            elif A[i,0] != 0:
                plt.axvline(x=b[i] / A[i, 0], color='r' if (mark_last_ineq and i==m-1) else 'g', linestyle='-', lw=2, label=ineq_str(A,b,i))

        # Shade Feasible Region
        feasible_cmap = matplotlib.colors.ListedColormap(['white', 'lime'])
        feasible_region = np.all(A[:,0,None,None]*X + A[:,1,None,None]*Y >= b[:,None,None], axis=0)
        plt.contourf(X, Y, feasible_region, cmap=feasible_cmap, alpha=0.2)

        # Plot Integer Points
        x_int = np.arange(np.ceil(xlim[0]), np.floor(xlim[1])+1, 1)
        y_int = np.arange(np.ceil(ylim[0]), np.floor(ylim[1])+1, 1)
        X_int, Y_int = np.meshgrid(x_int, y_int)
        plt.scatter(X_int, Y_int, s=40, color='blue', label='Integers')

        # Plot Optimal Solution
        if sol is not None:
            plt.plot(sol[0],sol[1],'ro', markersize=10, label='(x,y)=(%g, %g)' % (sol[0],sol[1]))

        # Labels and Stuff
        plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(title)
        plt.grid(color='gray', linestyle='-', linewidth=1)
        plt.xticks(np.arange(np.floor(xlim[0]), xlim[1], 1))
        plt.yticks(np.arange(np.floor(ylim[0]), ylim[1], 1))
        plt.xlim(xlim)
        plt.ylim(ylim)

        if save_file:
            plt.savefig(save_file, bbox_inches='tight', transparent=True) # Save Plot
        plt.close(fig) # Don't Show Plot


if __name__ == '__main__':
    pass