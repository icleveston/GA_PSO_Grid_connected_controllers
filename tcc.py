#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Simulation import Simulation
from Validation.Validation import Validation

if __name__ == '__main__':

    #Create the validation object
    validation = Validation()

    #validation.plotFitness3D()
    #validation.plotFitnessContour()
    validation.runExperiment()
    #validation.validateNSGA2()
    #validation.validateMOPSO()
    #validation.compaeMO()
    # validation.summary('Validation/Results/04-18-2018-09-31-49.122773 - Experiment Schwefel', plot=False)

    #simulation = Simulation()

    # Execute the NSGA2 algorithm
    #pathNSGA2 = 'Simulation/04-11-2018-17-56-12 - NSGA2' #simulation.executeNSGA2()

    # Generate the report
    #simulation.getReport(pathNSGA2)

    # Execute the MOPSO algorithm
    #pathMOPSO = 'Simulation/04-13-2018-00-15-22 - MOPSO' #simulation.executeMOPSO()

    # Generate the report
    #simulation.getReport(pathMOPSO)

    # Generate comparison between the two algorithms
    #simulation.compare(pathNSGA2, pathMOPSO)

    # if action == 1:
    #
    #     #The new simulation path
    #     path = 'Simulation/' + time.strftime("%d-%m-%Y-%H-%M-%S")
    #
    #     # Create the simulation directory
    #     os.mkdir(path)
    #
    #     ga = Ga(evaluate, limInf=limInf, limSup=limSup, populationSize=500, path=path, weights=(-1, -1, -1))
    #
    #     hof = ga.run(nGenerations=20, crossOver=0.8, mutation=0.2, method='modified', saveGeneration=20)
    #     #hof = ga.run(nGenerations=100, crossOver=0.8, mutation=0.2, method='other', saveGeneration=20)
    #
    #     results = []
    #
    #     for ind in hof:
    #         res = {}
    #         res['ind'] = list(ind)
    #         res['val'] = evaluate(ind)
    #
    #         results.append(res)
    #
    #     output = open(path + '/results.pkl', 'wb')
    #     pickle.dump(results, output)
    #     output.close()
    #
    #     print("Finished")
    #