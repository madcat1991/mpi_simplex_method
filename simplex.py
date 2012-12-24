# coding: utf-8
u"""
    Решатель задач линейного программирования симплекс-методом.

    Использование:

        -- генерим задачу генератором из каталога problem_generator
        gena 50 100 100.txt
        -- решаем ее
        python simplex.py 100.txt

"""
import sys
import math
import time
from copy import deepcopy


# выводить описание задачи перед решением? (на больших задачах разумно отключить)
DUMP_TASK = True
# печатать на каждом шаге, какой столбец и строка выбираются?
PRINT_EVERY_STEP = True

def help():
    print __doc__

def with_sign(n):
    if n < 0:
        return '%6.3f' % n
    return '+%6.3f' % n

def print_tableau(tableau):
    print '['
    for row in tableau:
        print repr(row)
    print ']'

def get_n_from_tableau(tableau):
    m = len(tableau) - 1
    n = len(tableau[0]) - 1 - m
    return n

class Task(object):
    def __init__(self):
        self.n = 0  # number of variables
        self.m = 0  # number of constraints
        self.c = []
        self.A = []
        self.b = []

    def vars(self):
        for i in xrange(self.n):
            yield 'x%d' % i

    def eval(self, x):
        return sum(map(lambda c, x: c*x, self.c, x))

    def dump(self):
        print 'maximize'
        print ' '.join(map(lambda c, x: '%s%s' % (with_sign(c), x), self.c, self.vars()))
        print 'subject to'
        for j in xrange(self.m):
            print ' '.join(map(lambda c, x: '%s%s' % (with_sign(c), x), self.A[j], self.vars())), '<=', with_sign(self.b[j])

    def to_tableau(self):
        tableau = []

        row = []
        row.append(0.0)
        row.extend(map(lambda x: -x, self.c))
        row.extend([0.0] * self.m)
        tableau.append(row)

        for i in xrange(self.m):
            row = []
            row.append(self.b[i])
            row.extend(self.A[i])
            tail = [0.0] * self.m
            tail[i] = 1.0
            row.extend(tail)
            tableau.append(row)
        return tableau

    def dump_tableau(self):
        tableau = self.to_tableau()
        print_tableau(tableau)


def pivotOn(tableu, row, col):
    j = 0
    pivot = tableu[row][col]
    for x in tableu[row]:
        # updating rowth row
        tableu[row][j] = tableu[row][j] / pivot
        j += 1
    i = 0
    # for each row
    for xi in tableu:
        if i != row: # ...except the one just updated
            ratio = xi[col]
            j = 0
            # update all columns of this row
            for xij in xi:
                #xij -= ratio * tableu[row][j]
                #tableu[i][j] = xij
                tableu[i][j]  -= ratio * tableu[row][j]
                j += 1
        i += 1
    return tableu


def solve_simplex(tableu):
    THETA_INFINITE = -1
    opt = False
    unbounded = False
    n = len(tableu[0])
    m = len(tableu) - 1
    
    step_no = 0
    while ((not opt) and (not unbounded)):
        step_no += 1
        min = 0.0
        pivotCol = j = 0
        ########### COLUMN SELECTION
        # Requires: entire row 0
        while(j < (n-m)):
            # searching 0th row, selecting pivotCol
            cj = tableu[0][j]
            if (cj < min) and (j > 0):
                 min = cj
                 pivotCol = j
            j += 1   
        if min == 0.0:
            print 'step', step_no, 'optimal'
            opt = True
            continue

        ########## ROW SELECTION
        # Requires: column 0, selected column
        pivotRow = i = 0
        minTheta = THETA_INFINITE
        for xi in tableu:
            # searching 0th and pivotCol'th columns, selecting pivotRow
            if (i > 0):
                xij = xi[pivotCol]
                if xij > 0:
                    theta = (xi[0] / xij)
                    if (theta < minTheta) or (minTheta == THETA_INFINITE):
                        minTheta = theta
                        pivotRow = i
            i += 1
        if minTheta == THETA_INFINITE:
            unbounded = True
            continue
        
        #now we pivot on pivotRow and pivotCol
        print 'step', step_no, 'pivot x', pivotCol-1, 'row', pivotRow
        tableu = pivotOn(tableu, pivotRow, pivotCol)
        #print_tableau(tableu)
    print 'opt = %s' % opt
    print 'unbounded = %s' % unbounded
    return tableu


#============================================ MODEL OF PARALLEL IMPL =====================

class Worker(object):
    def __init__(self, tableau, x_lo, x_hi):
        self.x_lo = x_lo
        self.x_hi = x_hi
        self.m = len(tableau) - 1
        self.n = len(tableau[0]) - 1 - self.m
        self.fill_tab_from_tableau(tableau)

        for xi in xrange(0, self.n + self.m):
            if not (x_lo <= xi <= x_hi):
                self.clear_column(xi)

    def fill_tab_from_tableau(self, tableau):
        self.tab = [None] * (self.n + self.m + 1)  # list of columns
        for xi in [-1] + range(self.x_lo, self.x_hi + 1):
            self.tab[xi + 1] = [None] * (self.m + 1)
            for j in xrange(self.m + 1):
                self.tab[xi + 1][j] = (tableau[j][xi + 1])

    def get_column_for_x(self, xi):
        assert self.x_lo <= xi <= self.x_hi or xi == -1
        col = self.tab[xi + 1]
        return deepcopy(col)

    def receive_column(self, xi, col):
        self.tab[xi + 1] = deepcopy(col)

    def clear_column(self, xi):
        self.tab[xi + 1] = None

    def patch_in_tableau(self, tableau):
        for i in xrange(0, self.n + self.m + 1):
            if i != 0 and not (self.x_lo <= i-1 <= self.x_hi):
                continue
            for j in xrange(0, self.m + 1):
                tableau[j][i] = self.tab[i][j]

    def select_col_and_row(self):
        """ returns (xi, tab[0][xi+1], row) """
        selected_xi = -1
        max_xi = min(self.n - 1, self.x_hi)
        min_v = 0.0
        for xi in xrange(self.x_lo, max_xi + 1):
            if self.tab[xi + 1][0] < min_v:
                selected_xi = xi
                min_v = self.tab[xi + 1][0]
        if selected_xi == -1:
            return (-1, 0, 0, 0)

        selected_row = -1
        ratio = 0
        for j in xrange(1, self.m + 1):
            xij = self.tab[selected_xi + 1][j]
            if xij > 0:
                ratio = self.tab[0][j] / xij
                if (selected_row == -1) or ratio < min_ratio:
                    selected_row = j
                    min_ratio = ratio
        return (
            selected_xi, self.tab[selected_xi + 1][0], selected_row
        )

    def pivot(self, p_xi, p_row):
        # requires column 0, column p_xi
        pivot_value = self.tab[p_xi + 1][p_row]

        for xi in [-1] + range(self.x_lo, self.x_hi + 1):
            self.tab[xi + 1][p_row] /= pivot_value

        for row in xrange(0, self.m + 1):
            if row == p_row:
                continue
            ratio = self.tab[p_xi + 1][row]
            for xi in [-1] + range(self.x_lo, self.x_hi + 1):
                self.tab[xi + 1][row] -= ratio * self.tab[xi + 1][p_row]

        #print self.get_column_for_x(-1)

        if not (self.x_lo <= p_xi <= self.x_hi):
            self.clear_column(p_xi)


def split_x_between_n(x, n):
    parts = []
    n_per_part = int(math.ceil(x / float(n)))
    for i in xrange(n):
        i0 = n_per_part * i
        i1 = n_per_part * (i + 1) - 1
        parts.append((i0, min(x-1, i1)))
    return parts


def parallel_simplex(tableau, n_workers):
    n = get_n_from_tableau(tableau)
    m = len(tableau) - 1
    parts = split_x_between_n(n + m, n_workers)
    workers = dict([
        (i, Worker(tableau, *part))
        for i, part in enumerate(parts)
    ])

    def choose_best_vote(votes):
        best_worker = -1
        min_vote = (-1, 0, 0, 0)
        for wi, v in votes:
            if v[0] != -1 and v[2] != -1:
                if best_worker == -1 or v[1] < min_vote[1]:
                    best_worker = wi
                    min_vote = v
        return (best_worker, min_vote[0], min_vote[2])

    step_no = 0
    while True:
        step_no += 1
        col_row_votes = [(i, w.select_col_and_row()) for i, w in workers.iteritems()]
        best_wi, p_xi, p_row = choose_best_vote(col_row_votes)
        if best_wi == -1:
            print 'step', step_no, 'optimal'
            break

        if PRINT_EVERY_STEP:
            print 'step', step_no, 'pivot x', p_xi, 'row', p_row
        p_col = workers[best_wi].get_column_for_x(p_xi)
        for wi in workers:
            if wi == best_wi: continue
            workers[wi].receive_column(p_xi, p_col)

        for wi in workers:
            workers[wi].pivot(p_xi, p_row)

        """
        blank_tableau = [
            [0.0] * (n + m + 1)
            for j in xrange(len(tableau))
        ]

        for wi in workers:
            workers[wi].patch_in_tableau(blank_tableau)
        print_tableau(blank_tableau)
        """

    blank_tableau = [
        [0.0] * (n + m + 1)
        for j in xrange(len(tableau))
    ]

    for wi in workers:
        workers[wi].patch_in_tableau(blank_tableau)
    return blank_tableau


def get_solution_from_tableau(tableau):
    n = len(tableau[0]) - 1
    m = len(tableau) - 1
    x = []
    for xi in xrange(n):
        i = xi + 1
        if tableau[0][i] == 0:
            _x = 0
            for bi in xrange(m):
                j = bi + 1
                _x += tableau[j][0] * tableau[j][i]
        else:
            _x = 0
        x.append(_x)
    return x


def read_task(f):
    fint = lambda: int(f.readline().strip())
    ffloat = lambda: float(f.readline().strip())
    t = Task()
    t.n = fint()
    t.m = fint()
    for i in xrange(t.n):
        t.c.append(ffloat())
    for j in xrange(t.m):
        t.A.append([])
        for i in xrange(t.n):
            t.A[-1].append(ffloat())
        t.b.append(ffloat())
    return t


def show_input(filename):
    f = open(filename)
    t = read_task(f)
    t.dump()
    if 'tab' in sys.argv:
        tableau = t.to_tableau()
        print_tableau(tableau)


def test_parser(filename):
    f = open(filename)
    t = read_task(f)
    if DUMP_TASK:
        t.dump()
    #t.dump_tableau()

    tableau = t.to_tableau()
    time1 = time.time()
    tableau = solve_simplex(tableau)
    time2 = time.time()
    #print_tableau(tableau)
    print ''
    print ''
    print '================SOLUTION============='
    print ''
    print ''
    xs = get_solution_from_tableau(tableau)[:t.n]
    print xs
    print 'F =', t.eval(xs)
    print 'got solution in %.2f seconds' % (time2 - time1)


def test_parallel(filename):
    f = open(filename)
    t = read_task(f)
    tableau = t.to_tableau()
    print 'orig============'
    #print_tableau(tableau)

    ref_t = solve_simplex(deepcopy(tableau))
    print 'ref============='
    #print_tableau(ref_t)
    ref_xs = get_solution_from_tableau(ref_t)[:t.n]
    print ref_xs

    par_t = parallel_simplex(deepcopy(tableau), 3)
    print 'par============='
    #print_tableau(par_t)
    par_xs = get_solution_from_tableau(par_t)[:t.n]
    print par_xs

    print ref_xs == par_xs




if __name__ == '__main__':
    if len(sys.argv) == 1:
        help()
        exit()
    filename = sys.argv[-1]
    if not filename.endswith('.txt'):
        filename = 'input.txt'
    if 'show' in sys.argv:
        show_input(filename)
    elif 'par' in sys.argv:
        test_parallel(filename)
    else:
        test_parser(filename)

