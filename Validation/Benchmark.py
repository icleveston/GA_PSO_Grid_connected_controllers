from deap import benchmarks
from math import sin

class Benchmark:

    def __init__(self, nvar=2):

        self.fitness = None
        self.name = None
        self.limitInf = None
        self.limitSup = None
        self.step = None
        self.nvar = nvar

    def getRosenbrockInfo(self):

        self.fitness = self._rosenbrock
        self.name = "Rosenbrock"
        self.limitInf = [-2]*self.nvar
        self.limitSup = [2]*self.nvar
        self.step = 0.1

        return self

    def getSchafferInfo(self):

        self.fitness = self._schaffer
        self.name = "Schaffer"
        self.limitInf = [-100]*self.nvar
        self.limitSup = [100]*self.nvar
        self.step = 0.25

        return self

    def getHimmelblauInfo(self):

        self.fitness = self._himmelblau
        self.name = "Himmelblau"
        self.limitInf = [-6]*self.nvar
        self.limitSup = [6]*self.nvar
        self.step = 0.1

        return self

    def getSchwefelInfo(self):

        self.fitness = self._schwefel
        self.name = "Schwefel"
        self.limitInf = [-500]*self.nvar
        self.limitSup = [500]*self.nvar
        self.step = 10

        return self

    def getRastriginInfo(self):

        self.fitness = self._rastrigin
        self.name = "Rastrigin"
        self.limitInf = [-5]*self.nvar
        self.limitSup = [5]*self.nvar
        self.step = 0.1

        return self

    def getBohachevskyInfo(self):

        self.fitness = self._bohachevsky
        self.name = "Bohachevsky"
        self.limitInf = [-100]*self.nvar
        self.limitSup = [100]*self.nvar
        self.step = 0.5

        return self

    def getAckleyInfo(self):

        self.fitness = self._ackley
        self.name = "Ackley"
        self.limitInf = [-30]*self.nvar
        self.limitSup = [30]*self.nvar
        self.step = 0.5

        return self

    def getHauptInfo(self):

        self.fitness = self._haupt
        self.name = "Haupt"
        self.limitInf = 0
        self.limitSup = 10
        self.step = 0.05

        return self


    def _haupt(self, individual):
        # y = x1*sin(4*x1) + 1.1*x2*sin(2*x2)
        return individual[1] * sin(4 * individual[1]) + 1.1 * individual[0] * sin(2 * individual[0]),

    def _schwefel(self, sol):
        return benchmarks.schwefel(sol)

    def _himmelblau(self, sol):
        return benchmarks.himmelblau(sol)

    def _rastrigin(self, sol):
        return benchmarks.rastrigin(sol)

    def _rosenbrock(self, sol):
        return benchmarks.rosenbrock(sol)

    def _schaffer(self, sol):
        return benchmarks.schaffer(sol)

    def _ackley(self, sol):
        return benchmarks.ackley(sol)

    def _bohachevsky(self, sol):
        return benchmarks.bohachevsky(sol)

