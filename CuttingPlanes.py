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

def ineq_str(a, r):
    cx, cy =-a[0], -a[1]
    s = [t for t in [coeff_str(cx, 'x'), coeff_str(cy, 'y')] if t]
    return ' + '.join(s) + ' $\leq$ ' + str(int(-r)) if r%1==0 else str(-r)

class MILP:
    def __init__(self, n, n1):
        self.c = np.zeros((n))
        self.A = np.zeros((0,n))
        self.b = np.zeros((0))
        self.B = np.zeros((0,n))
        self.d = np.zeros((0))
        self.n, self.n1, self.m, self.p = n, n1, 0, 0

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

class Plotter:
    def __init__(self):
        pass

    def init_figure(self, milp:MILP, cuts=(np.array([[]]),np.array([]))) -> plt.Figure:
        A,b = (milp.A,milp.b) if cuts[1].shape[0]==0 else (np.vstack([milp.A, cuts[0]]), np.append(milp.b, cuts[1]))
        vs = lineqs(milp.A, milp.b, show=False)
        self.xmin, self.xmax = [-0, 1.1*np.max(vs[:,0])]
        self.ymin, self.ymax = [-0, 1.1*np.max(vs[:,1])] 
        self.x_vals = np.linspace(self.xmin, self.xmax, 400)
        self.y_vals = np.linspace(self.ymin, self.ymax, 400)
        self.X, self.Y = np.meshgrid(self.x_vals, self.y_vals)

        ratio = (self.xmax - self.xmin) / (self.ymax - self.ymin)
        return plt.figure(figsize=(10,10/ratio))

    def plot_objective_function(self, c, levels=20, origin='lower', cmap='viridis', alpha_shade=0.1, alpha_isoline=0.3):
        '''
        Plots the objective function c[0]x + c[1]y
        '''
        Z = c[0]*self.X + c[1]*self.Y
        plt.imshow(Z, extent=[self.xmin, self.xmax, self.ymin, self.ymax], origin=origin, cmap=cmap, alpha=alpha_isoline)
        plt.contour(self.X, self.Y, Z, levels=levels, cmap=cmap, alpha=alpha_shade)

    def plot_line(self, a, r, color='g', line_width=2, title=''):
        '''
        Plots the line a[0]x + a[1]y = r
        '''
        if a[1] != 0:
            plt.plot(self.x_vals, (r - a[0]*self.x_vals) / a[1], color+'-', lw=line_width, label=title)
        elif a[0] != 0:
            plt.axvline(x=r / a[0], color=color, linestyle='-', lw=line_width, label=title)

    def plot_feasible_region(self, milp:MILP, cuts=(np.array([[]]),np.array([]))):
        '''
        Shades the feasible region A(x,y)' >= b and cuts[0](x,y)' >= cuts[1]
        '''
        A,b = (milp.A,milp.b) if cuts[1].shape[0]==0 else (np.vstack([milp.A, cuts[0]]), np.append(milp.b, cuts[1]))
        feasible_cmap = matplotlib.colors.ListedColormap(['white', 'lime'])
        feasible_region = np.all(A[:,0,None,None]*self.X + A[:,1,None,None]*self.Y >= b[:,None,None], axis=0)
        plt.contourf(self.X, self.Y, feasible_region, cmap=feasible_cmap, alpha=0.2)

    def plot_integer_points(self, n1:int, color='blue'):
        '''
        Highlights the integer constraints.
        If n1=2, x and y must be integer
        If n1=1, only x must be integer
        '''
        if n1 == 2:
            x_int = np.arange(np.ceil(self.xmin), np.floor(self.xmax)+1, 1)
            y_int = np.arange(np.ceil(self.ymin), np.floor(self.ymax)+1, 1)
            X_int, Y_int = np.meshgrid(x_int, y_int)
            plt.scatter(X_int, Y_int, s=40, color=color)
        elif n1 == 1:
            for x in range(int(np.ceil(self.xmin)), int(np.floor(self.xmax)+1)):
                plt.axvline(x=x, color=color, linestyle='-', alpha=1, lw=2)

    def plot_point(self, x, y, color='r'):
        plt.plot(x, y, color+'o', markersize=10, label='(x,y)=(%g, %g)' % (x, y))


    def plot_milp(self, milp:MILP, title:str='MILP', optimal_solution=None, cuts=(np.array([[]]),np.array([])), show_legend=True, save_file:str=''):
        A,b = (milp.A,milp.b) if cuts[1].shape[0]==0 else (np.vstack([milp.A, cuts[0]]), np.append(milp.b, cuts[1]))

        # Initialize figure
        fig = self.init_figure(milp, cuts)

        # Plot objective function
        self.plot_objective_function(milp.c)

        # Plot Inequality Constraints
        for i in range(len(b)):
            self.plot_line(A[i], b[i], color='g', title=ineq_str(A[i], b[i]))

        # Shade Feasible Region
        self.plot_feasible_region(milp, cuts)

        # Plot Integer Points
        self.plot_integer_points(milp.n1)

        # Plot Optimal Solution
        if optimal_solution is not None:
            self.plot_point(optimal_solution[0], optimal_solution[1])

        # Labels and Stuff
        if show_legend:
            plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
            plt.xlabel('x')
            plt.ylabel('y')
        plt.title(title)
        plt.grid(color='gray', linestyle='--', linewidth=1)
        plt.xticks(np.arange(np.floor(self.xmin), self.xmax, 1))
        plt.yticks(np.arange(np.floor(self.ymin), self.ymax, 1))
        plt.xlim((self.xmin, self.xmax))
        plt.ylim((self.ymin, self.ymax))

        if save_file:
            plt.savefig(save_file, bbox_inches='tight', transparent=True) # Save Plot
        plt.close(fig) # Don't Show Plot


if __name__ == '__main__':
    pass