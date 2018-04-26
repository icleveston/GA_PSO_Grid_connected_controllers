#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing as mp
import operator
from deap import creator, base, tools
from pyDOE import *
import random
import pickle
from datetime import datetime
import copy


class PSO:

    def __init__(self, fitnessFunction, limInf, limSup, path, weights, phi1=0.4, phi2=0.6, smin=-0.3, smax=0.3,
                 populationSize=100, multiprocessing=False):

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
        self.multiprocessing = multiprocessing

        creator.create("FitnessMulti", base.Fitness, weights=self.weights)
        creator.create("Particle", list, fitness=creator.FitnessMulti, speed=list, smin=None, smax=None, best=None)

        self.toolbox = base.Toolbox()

        if self.multiprocessing:
            self.pool = mp.Pool()
            self.toolbox.register("map", self.pool.map)

    def _generate(self, size, pmin, pmax, smin, smax):

        # Create the particle
        part = creator.Particle(random.uniform(0, 500) for x in range(size))

        # Generate the particle speed
        part.speed = [random.uniform(smin, smax) for _ in range(size)]

        # Set the particle speed
        part.smin = smin
        part.smax = smax

        return part

    def _updateParticle(self, part, best, phi1, phi2, w):

        u1 = (w*random.uniform(0, phi1) for _ in range(len(part)))
        u2 = (w*random.uniform(0, phi2) for _ in range(len(part)))
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
        if part[0] > self.limSup[0] or part[0] < self.limInf[0] or part[1] > self.limSup[1] or part[1] < self.limInf[1]:
            part[:] = [random.uniform(self.limInf[x], self.limSup[x]) for x in range(len(self.limSup))]

    def run(self, nGenerations=1000, saveEpoch=5, method='modified', verbose=True):

        # Start time
        start_time = datetime.now()

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
        evaluationTotal = 0
        count = 0
        counter = 0
        w = 1
        g = 0
        minRaio = 0

        while (method == 'modified' and g < 5000) or (g < nGenerations and method != 'modified'):

            popFit = map(lambda x: list(x), pop)

            fitnesses = self.toolbox.map(self.toolbox.evaluate, popFit)

            for ind, fit in zip(pop, fitnesses):
                ind.fitness.values = fit

            for part in pop:

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

            if verbose:
                print(logbook.stream + "\t\t" + str(best.fitness.values))

            if count == saveEpoch:

                results.append(best.fitness.values[0])
                count = 0

            if method == 'modified':

                # Verifica estagnacao
                if round(best.fitness.values[0], 3) <= minRaio:
                    counter = counter + 1
                else:
                    counter = 0

                # Atualiza as referencias
                minRaio = round(best.fitness.values[0], 3)

                # Se estourou o limite, termina execucao
                if counter >= nGenerations:
                    break

            count = count + 1

            w = w - w / (nGenerations)
            g = g + 1

        # End time
        end_time = datetime.now()

        elapsed_time = end_time - start_time

        # Save parameters
        file = open(self.path + '/parameters.txt', "w")
        file.write("Method: PSO\n")
        file.write("N Population: " + str(self.populationSize) + "\n")
        file.write("N Generation: " + str(nGenerations) + "\n")
        file.write("phi1: " + str(self.phi1) + "\n")
        file.write("phi2: " + str(self.phi2) + "\n")
        file.write("Elapsed Time: " + str(elapsed_time) + "\n")
        file.close()

        # Save generation
        output = open(self.path + "/statsPopulation.pkl", 'wb')
        pickle.dump(results, output)
        output.close()

        # Select the historic
        genMin, genMean, genMax = logbook.select("min", 'mean', 'max')

        # Save historic
        output = open(self.path + "/historic.pkl", 'wb')
        pickle.dump([genMin, genMean, genMax], output)
        output.close()

        info = {
            'bestInd': list(best),
            'bestVal': self.fitnessFunction(best),
            'evalTotal': evaluationTotal,
            'elapsedTime': str(elapsed_time),
            'nGeneration': len(genMin)
        }

        # Save additional information
        output = open(self.path + "/info.pkl", 'wb')
        pickle.dump(info, output)
        output.close()

        return best
