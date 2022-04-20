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

    # TODO
    println("In file resolution.jl, in method cplexSolve(), TODO: fix input and output, define the model")
    n = size(limits,2)
    half_n = Int64(floor(n/2))
    #Defining variables
    @variable(m, x[1:n, 1:n, 1:n], Bin) #Variable 1 si case i,j = k 0 sinon
    @variable(m, v[1:n, 1:n, 1:4], Bin) #Variable 1 si case i,j est visible en visant de k
    @variable(m, case[1:n, 1:n] >= 0, Int)

    #Defining contraints
    @constraint(m, [i in 1:n, j in 1:n], sum(x[i,j,k] for k in 1:n)==1) #1 valeur par case
    @constraint(m, [i in 1:n, k in 1:n], sum(x[i,j,k] for j in 1:n)==1) #1 valeur par ligne
    @constraint(m, [j in 1:n, k in 1:n], sum(x[i,j,k] for i in 1:n)==1) #1 valeur par colonne

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

    ## Prendre le valeur de la case
    @constraint(m, [i in 1:n, j in 1:n], sum(k*x[i, j, k] for k in 1:n) == case[i,j])

    ## Define les non visibles
    ## Perspective 1
    #constraint(m, [i in 2:n, j in 1:n], sum(x[i,j,k] for k in 1:(i-1)) + v[i,j,1] <= 1)
    ## Perspective 2
    #@constraint(m, [i in 1:n, j in 1:(n-1)], sum(x[i,j,k] for k in 1:(n-j)) + v[i,j,2] <= 1)
    ## Perspective 3
    #@constraint(m, [i in 1:(n-1), j in 1:n], sum(x[i,j,k] for k in 1:(n-i)) + v[i,j,3] <= 1)
    ## Perspective 4
    #@constraint(m, [i in 1:n, j in 2:n], sum(x[i,j,k] for k in 1:(j-1)) + v[i,j,4] <= 1)

    ## Perspective 1
    #@constraint(m, [i in 2:n, j in 1:n, k in i:n], x[i,j,k] +
     #   sum(x[i2,j,k2] for i2 in 1:(i-1) for k2 in (k+1):n) + v[i,j,1] <= 2)

    ## Perspective 2
    #@constraint(m, [i in 1:n, j in 1:(n-1), k in (n-j+1):n], x[i,j,k] +
    #    sum(x[i,j2,k2] for j2 in (j+1):n for k2 in (k+1):n) + v[i,j,2] <= 2)

    ## Perspective 3
    #@constraint(m, [i in 1:(n-1), j in 1:n, k in (n-i+1):n], x[i,j,k] +
    #    sum(x[i2,j,k2] for i2 in (i+1):n for k2 in (k+1):n) + v[i,j,3] <= 2)

    ## Perspective 4
    #@constraint(m, [i in 1:n, j in 2:n, k in j:n], x[i,j,k] +
        #sum(x[i,j2,k2] for j2 in 1:(j-1) for k2 in (k+1):n) + v[i,j,4] <= 2)

    # Perspective 1
    @constraint(m, [i in 2:n, i2 in 1:(i-1), j in 1:n], case[i2, j] - case[i, j] <= (1 - v[i,j,1])*n)

    #Perspective 2
    @constraint(m, [i in 1:n, j in 1:(n-1), j2 in (j+1):n], case[i, j2] - case[i, j] <= (1 - v[i,j,2])*n)

    #Perspective 3
    @constraint(m, [i in 1:(n-1), i2 in (i+1):n, j in 1:n], case[i2, j] - case[i, j] <= (1 - v[i,j,3])*n)

    #Perspective 4
    @constraint(m, [i in 1:n, j in 2:n, j2 in 1:(j-1)], case[i, j2] - case[i,j] <= (1 - v[i,j,4])*n)

    @objective(m, Max, case[1,1])
    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(m)


    if JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT
        vCase = Array{Int64}(JuMP.value.(case))
        displaySolution(limits, vCase)
        vV = Array{Int64}(JuMP.value.(v))
        print([vV[i, :, 1] for i in 1:n])
    end
    # Return:
    # 1 - true if an optimum is found
    # 2 - the resolution time
    return JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start
    
end

"""
Heuristically solve an instance
"""
function heuristicSolve()

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
        readInputFile(dataFolder * file)

        # TODO
        println("In file resolution.jl, in method solveDataSet(), TODO: read value returned by readInputFile()")
        
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
                    
                    # TODO 
                    println("In file resolution.jl, in method solveDataSet(), TODO: fix cplexSolve() arguments and returned values")
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime = cplexSolve()
                    
                    # If a solution is found, write it
                    if isOptimal
                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write cplex solution in fout") 
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Start a chronometer 
                    startingTime = time()
                    
                    # While the grid is not solved and less than 100 seconds are elapsed
                    while !isOptimal && resolutionTime < 100
                        
                        # TODO 
                        println("In file resolution.jl, in method solveDataSet(), TODO: fix heuristicSolve() arguments and returned values")
                        
                        # Solve it and get the results
                        isOptimal, resolutionTime = heuristicSolve()

                        # Stop the chronometer
                        resolutionTime = time() - startingTime
                        
                    end

                    # Write the solution (if any)
                    if isOptimal

                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write the heuristic solution in fout")
                        
                    end 
                end

                println(fout, "solveTime = ", resolutionTime) 
                println(fout, "isOptimal = ", isOptimal)
                
                # TODO
                println("In file resolution.jl, in method solveDataSet(), TODO: write the solution in fout") 
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end
