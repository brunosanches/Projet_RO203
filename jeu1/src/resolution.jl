# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX

include("generation.jl")
include("io.jl")

TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function cplexSolve(limits::Array{Int64,2})

    # Create the model
    m = Model(with_optimizer(CPLEX.Optimizer))

    n = size(limits,2)

    #Defining variables
    @variable(m, x[1:n, 1:n, 1:n], Bin) # Variable 1 si tour i,j = k 0 sinon
    @variable(m, v[1:n, 1:n, 1:4], Bin) # Variable 1 si tour i,j est visible en visant de k
    @variable(m, g_l[1:n, 1:n, 1:n], Bin) # 1 si tour est plus grande qu'autre dans la meme ligne
    @variable(m, g_c[1:n, 1:n, 1:n], Bin) # 1 si tour est plus grande qu'autre dans la meme colonne
    @variable(m, tour[1:n, 1:n] >= 0, Int) # valeur d'hauteur de la tour

    #Defining contraints
    @constraint(m, [i in 1:n, j in 1:n], sum(x[i,j,k] for k in 1:n)==1) # 1 valeur par tour
    @constraint(m, [i in 1:n, k in 1:n], sum(x[i,j,k] for j in 1:n)==1) # 1 valeur par ligne
    @constraint(m, [j in 1:n, k in 1:n], sum(x[i,j,k] for i in 1:n)==1) # 1 valeur par colonne

    ## Constraintes de visualisation (combien de tours je peux voir)
    @constraint(m, [k in [1,3], j in 1:n], sum(v[i,j,k] for i in 1:n)==limits[k,j])
    @constraint(m, [k in [2,4], i in 1:n], sum(v[i,j,k] for j in 1:n)==limits[k,i])

    ## Premiere tour toujours visible
    @constraint(m, [j in 1:n], v[1, j, 1] == 1)
    @constraint(m, [i in 1:n], v[i, n, 2] == 1)
    @constraint(m, [j in 1:n], v[n, j, 3] == 1)
    @constraint(m, [i in 1:n], v[i, 1, 4] == 1)

    ## Si tour a n de hauteur, visible après toutes les perspectives
    @constraint(m, [i in 1:n, j in 1:n, k in 1:4], v[i,j,k] >= x[i,j,n])

    ## Prendre le valeur de la tour
    @constraint(m, [i in 1:n, j in 1:n], sum(k*x[i, j, k] for k in 1:n) == tour[i,j])

    ## Impose les bonnes valeurs pour les variables g_l et g_c
    @constraint(m, [i in 1:n, j in 1:n, j2 in 1:n], tour[i, j2] - tour[i,j] <= (1-g_l[i,j,j2])*n)
    @constraint(m, [i in 1:n, j in 1:n, j2 in 1:n], tour[i, j2] - tour[i,j] >= -g_l[i,j,j2]*n)

    @constraint(m, [i in 1:n, j in 1:n, i2 in 1:n], tour[i2, j] - tour[i,j] <= (1-g_c[i,j,i2])*n)
    @constraint(m, [i in 1:n, j in 1:n, i2 in 1:n], tour[i2, j] - tour[i,j] >= -g_c[i,j,i2]*n)

    # Perspective 1
    @constraint(m, [i in 2:n, j in 1:n], sum(g_c[i, j, i2] for i2 in 1:(i-1)) >= (i-1)*v[i,j,1])
    @constraint(m, [i in 2:n, j in 1:n], sum(g_c[i, j, i2] for i2 in 1:(i-1)) <= i-2 + v[i,j,1])

    # Perspective 2
    @constraint(m, [i in 1:n, j in 1:(n-1)], sum(g_l[i, j, j2] for j2 in (j+1):n) >= (n-j)*v[i,j,2])
    @constraint(m, [i in 1:n, j in 1:(n-1)], sum(g_l[i, j, j2] for j2 in (j+1):n) <= (n-j)-1 + v[i,j,2])

    # Perspective 3
    @constraint(m, [i in 1:(n-1), j in 1:n], sum(g_c[i, j, i2] for i2 in (i+1):n) >= (n-i)*v[i,j,3])
    @constraint(m, [i in 1:(n-1), j in 1:n], sum(g_c[i, j, i2] for i2 in (i+1):n) <= (n-i)-1 + v[i,j,3])

    # Perspective 4
    @constraint(m, [i in 1:n, j in 2:n], sum(g_l[i, j, j2] for j2 in 1:(j-1)) >= (j-1)*v[i,j,4])
    @constraint(m, [i in 1:n, j in 2:n], sum(g_l[i, j, j2] for j2 in 1:(j-1)) <= j-2 + v[i,j,4])

    @objective(m, Max, tour[1,1])
    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(m)


    if JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT
        vtour = Array{Int64}(floor.(JuMP.value.(tour)))
        displaySolution(limits, vtour)
    end
    # Return:
    # 1 - true if an optimum is found
    # 2 - the resolution time
    return JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, vtour
    
end

"""
Heuristically solve an instance
"""
function heuristicSolve(limits::Array{Int64,2})

    # TODO
    println("In file resolution.jl, in method heuristicSolve(), TODO: fix input and output, define the model")
    
end 

"""
Solve all the instances contained in "../data" through CPLEX and heuristics

The results are written in "../res/cplex" and "../res/heuristic"

Remark: If an instance has previously been solved (either by cplex or the heuristic) it will not be solved again
"""
function solveDataSet()

    dataFolder = "../data/"
    resFolder = "../res/"

    # Array which contains the name of the resolution methods
    resolutionMethod = ["cplex"]
    #resolutionMethod = ["cplex", "heuristique"]

    # Array which contains the result folder of each resolution method
    resolutionFolder = resFolder .* resolutionMethod

    # Create each result folder if it does not exist
    for folder in resolutionFolder
        if !isdir(folder)
            mkdir(folder)
        end
    end
            
    global isOptimal = false
    global solveTime = -1

    # For each instance
    # (for each file in folder dataFolder which ends by ".txt")
    for file in filter(x->occursin(".txt", x), readdir(dataFolder))  
        
        println("-- Resolution of ", file)
        instance = readInputFile(dataFolder * file)
        
        # For each resolution method
        for methodId in 1:size(resolutionMethod, 1)
            
            outputFile = resolutionFolder[methodId] * "/" * file

            # If the instance has not already been solved by this method
            if !isfile(outputFile)
                
                fout = open(outputFile, "w")  

                resolutionTime = -1
                isOptimal = false
                
                # If the method is cplex
                if resolutionMethod[methodId] == "cplex"
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime, solution = cplexSolve(instance)
                    
                    # If a solution is found, write it
                    if isOptimal
                        writeSolution(fout, solution)
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Start a chronometer 
                    startingTime = time()
                    
                    # While the grid is not solved and less than 100 seconds are elapsed
                    while !isOptimal && resolutionTime < 100
                        
                        # Solve it and get the results
                        isOptimal, resolutionTime, solution = heuristicSolve(instance)

                        # Stop the chronometer
                        resolutionTime = time() - startingTime
                        
                    end

                    # Write the solution (if any)
                    if isOptimal
                        writeSolution(fout, solution)
                    end 
                end

                println(fout, "solveTime = ", resolutionTime) 
                println(fout, "isOptimal = ", isOptimal)
                
                # TODO Necessary ?
                #println("In file resolution.jl, in method solveDataSet(), TODO: write the solution in fout") 
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end
