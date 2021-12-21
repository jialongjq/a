Per compilar basics, greedy, local_search i metaheuristic:

	make

Per executar basics:
	./basics -i [path de l'arxiu input (graf)] -r [path de l'arxiu del conjunt dominador]

Per executar greedy:

	./greedy -i [path de l'arxiu input]

Per executar cerca local (hill climbing first/best improvement):

	./local_search -i [path de l'arxiu input] -best	<- Best improvement
	./local_search -i [path de l'arxiu input] -first	<- First improvement

Per executar metaheuristica (genetic algorithm):

	./metaheuristic -i [path de l'arxiu input] -n_apps [numero d'iteracions] -t [temps limit en segons] -n_population [mida de la poblacio]
	
Per compilar i executar CPLEX:
NOTA: s'ha d'especificar al Makefile del directori ILP_CPLEX la ruta en la qual esta instalada el software de CPLEX.

	cd ILP_CPLEX
	make
	./cplex -i [path de l'arxiu input] -t [temps limit en segons]