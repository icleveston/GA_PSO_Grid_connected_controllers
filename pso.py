#!/usr/bin/env python
# -*- coding: utf-8 -*-

import operator
from deap import creator, base, tools, benchmarks
import Validation
from pyDOE import *
import random
import copy


class PSO:

    def __init__(self, fitnessFunction, limInf, limSup, path, weights, phi1=0.4, phi2=0.6, smin=-0.3, smax=0.3,
                 populationSize=100):

        self.fitnessFunction = fitnessFunction
        self.limInf = limInf
        self.limSup = limSup
        self.phi1 = phi1
        self.phi2 = phi2
        self.smin = smin
        self.smax = smax
        self.populationSize = populationSize
        self.path = path
        self.weights = weights

        creator.create("FitnessMulti", base.Fitness, weights=self.weights)
        creator.create("Particle", list, fitness=creator.FitnessMulti, speed=list, smin=None, smax=None, best=None)

        self.toolbox = base.Toolbox()

    def _generate(self, size, pmin, pmax, smin, smax):

        # Create the particle
        part = creator.Particle(random.uniform(pmin[x], pmax[x]) for x in range(size))

        # Generate the particle speed
        part.speed = [random.uniform(smin, smax) for _ in range(size)]

        # Set the particle speed
        part.smin = smin
        part.smax = smax

        return part

    def _updateParticle(self, part, best, phi1, phi2, w):

        u1 = (random.uniform(0, phi1) for _ in range(len(part)))
        u2 = (random.uniform(0, phi2) for _ in range(len(part)))
        v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
        v_u2 = map(operator.mul, u2, map(operator.sub, best, part))

        part.speed = list(map(operator.add, part.speed, map(operator.add, v_u1, v_u2)))

        for i, speed in enumerate(part.speed):
            if speed < part.smin:
                part.speed[i] = part.smin
            elif speed > part.smax:
                part.speed[i] = part.smax

        part[:] = list(map(operator.add, part, part.speed))

        # Random particle's position if it exceeds the limit
        if part[0] > 10 or part[0] < 0 or part[1] > 10 or part[1] < 0:
            part[:] = [random.uniform(0, 10) for _ in range(2)]

    def run(self, nGenerations=1000, saveEpoch=5):

        self.toolbox.register("particle", self._generate, size=len(self.limSup), pmin=self.limInf, pmax=self.limSup,
                              smin=self.smin, smax=self.smax)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.particle)
        self.toolbox.register("update", self._updateParticle, phi1=self.phi1, phi2=self.phi2)
        self.toolbox.register("evaluate", self.fitnessFunction)

        pop = self.toolbox.population(n=self.populationSize)

        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("min", np.min, axis=0)
        stats.register("mean", np.mean, axis=0)
        stats.register("max", np.max, axis=0)

        logbook = tools.Logbook()
        logbook.header = ['gen', 'evals'] + (stats.fields if stats else [])

        best = None
        results = []
        count = 0
        w = 1

        res = {'pop': copy.deepcopy(pop), 'best': False}

        results.append(res)

        for g in range(nGenerations+1):

            for part in pop:

                part.fitness.values = self.toolbox.evaluate(part)

                if not part.best or part.best.fitness < part.fitness:
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values

                if not best or best.fitness < part.fitness:
                    best = creator.Particle(part)
                    best.fitness.values = part.fitness.values

            for part in pop:
                self.toolbox.update(part, best, w=w)

            # Gather all the fitnesses in one list and print the stats
            logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
            print(logbook.stream + "\t\t" + str(best.fitness.values))

            if count == saveEpoch:

                res = {'pop': copy.deepcopy(pop), 'best': copy.deepcopy(best)}

                results.append(res)
                count = 0

            count = count + 1

            w = w - w / (nGenerations)

        return pop, logbook, best, results
