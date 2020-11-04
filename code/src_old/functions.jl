"""
    getData()

Reads the CSV data file and returns the dataframe, hosts and viruses.
"""
function getData()
    # reading the data file
    df = CSV.read("./code/data/virionette.csv");
    # making a sorted list of unique hosts and viruses
    hosts = sort(unique(df.host_species));
    viruses = sort(unique(df.virus_genus));

    return df, hosts, viruses
end

"""
    buildInteractionMatrix()

Prepares the interaction matrix with initial data and returns the interaction
matrix, the hosts, viruses, and the dataframe.
"""
function buildInteractionMatrix()
    # getting the data
    (df, hosts, viruses) = getData()
    # filling the matrix with ones where an interaction occurs
    interaction_matrix = zeros(Float64, (length(viruses), length(hosts)))
    for interaction in eachrow(df)
        host_idx = findfirst(hosts .== interaction.host_species)
        virus_idx = findfirst(viruses .== interaction.virus_genus)
        interaction_matrix[virus_idx, host_idx] = 1.0
    end

    return interaction_matrix, hosts, viruses, df
end

"""
    crossValidation(targetedValue, initialValueMatrix, α, rank,
        suspected_newData, unlikely_newData)

Computes the imputation for a given matrix if the targeted value is 0.
Computes the leave one out validation for a given matrix if the targeted value
is 1.
"""
function crossValidation(targetedValue, initialValueMatrix, α, rank, suspected_newData, unlikely_newData)
    # printing the parameters of the current test
    println("Targeted Value: $(targetedValue) ")
    println("Rank: $(rank)")
    println("Alpha: $(α)")
    # building the interaction matrix
    (interaction_matrix, hosts, viruses, df) = buildInteractionMatrix()
    # targeting the positions to impute
    positions_to_impute = findall(interaction_matrix .== targetedValue)
    # copying the interaction matrix aimed to be modified
    output_matrix = copy(interaction_matrix)
    # doing the imputation for every targeted value
    @showprogress for position in positions_to_impute
        output_matrix[position] = imputation(interaction_matrix,
            position, initialValueMatrix, rank)
    end
    # visualizing the interactions
    display(plot(generateHeatmap("Initial Matrix \n(targeted value: $(targetedValue), rank: $(rank))", interaction_matrix),
        generateHeatmap("Output matrix", output_matrix,)))

    # Calculating the variation between the initial and final matrices
    #println("Variation: $(calculateVariation(interaction_matrix, output_matrix))%")

    # getting the top 10 interactions
    (maxValues, delta) = getTopInteractions(10, interaction_matrix, output_matrix, hosts, viruses, df)
    # generating the table of results
    generateResultsTable(maxValues, α, rank, suspected_newData, unlikely_newData, delta, hosts)
end

"""
    imputation(matrix, position, initialValueMatrix, rank; tolerance=1e-2, maxiter=50)

This function returns the imputed value at a given position (expressed in
CartesianCoordinates) in the matrix, seeded from a value in initialeValueMatrix. The imputation
is done by iterating an SVD at a given rank, and stops when the iteration
difference is smaller than the absolute tolerance, or after maxiter steps
have been done. Returns the final value for the current imputed position.
"""
function imputation(matrix, position, initialValueMatrix, rank; tolerance = 1e-2, maxiter = 50)
    # making a temporary copy of the initial matrix
    tempMatrix = copy(matrix)
    # updating the current position to impute from the initial matrix to
    # the initial value determined by the filter
    tempMatrix[position] = initialValueMatrix[position]
    # setting the Δ and initial iteration values
    Δ = 1.0
    iter = 1
    # iterating the SVD
    while Δ > tolerance
        # increasing the iteration value
        iter += 1
        # approximating the matrix
        approx_matrix = lowrank(tempMatrix, rank)
        # the change in value is the absolute difference between the
        # next and current iteration
        Δ = abs(approx_matrix[position] - tempMatrix[position])
        #  updating the temporary matrix at the current imputed position
        tempMatrix[position] = approx_matrix[position]
        # we stop if there are more than a set number of iterations
        iter ≥ maxiter && break
    end

    return tempMatrix[position]
end

"""
    lowrank(matrix, svd_rank)

Returns the low rank appproximation of the matrix, after a SVD.
"""
function lowrank(matrix, svd_rank)
    # making sure that the rank of the matrix is greater than the chosen rank
    @assert svd_rank ≤ rank(matrix)
    # factorizing the matrix
    factorization = svd(matrix)
    # fixing to zero all the singular values found after the chosen rank
    factorization.S[(svd_rank + 1):end] .= zero(eltype(factorization.S))
    # Note that the SVD method returns Vt, which is the transposed
    # right singular matrix

    return factorization.U * Diagonal(factorization.S) * factorization.Vt
end

"""
    generateHeatMap(title, matrix)

Returns a heatmap based on the given matrix.
"""
function generateHeatmap(title, matrix)

    return heatmap(matrix, xlabel = "Hosts", ylabel = "Viruses", c = :Greys,
        leg = true, frame = :box, title = title)
end

"""
    calculateVariation(initial_matrix, output_matrix)

Calculates the variation between two matrices.
"""
function calculateVariation(initial_matrix, output_matrix)
    # calculating the variation between two matrices
    delta = sum(abs.(initial_matrix .- output_matrix))

    return round((delta/length(initial_matrix)) * 100, digits = 2);
end

"""
    getTopInteractions(top, initial_matrix, output_matrix, hosts, viruses, df)

Returns the top imputed interactions with the highest probabilities of occurrence.

"""
function getTopInteractions(top, initial_matrix, output_matrix, hosts, viruses, df)
    # Initializing the minimum value
    minValue = findmin(output_matrix)
    # Creating a buffer to stock the maximum probability values and their indexes
    maxValues = CircularBuffer{Tuple{Float64,Int64,Int64}}(top)
    # Initializing the maximum value array to the minimum value
    push!(maxValues, (minValue[1],minValue[2][1],minValue[2][2]));
    # Iterating through the matrix
    for r in 1:size(output_matrix,1)
        for c in 1:size(output_matrix,2)
            # making sure that the present interaction was initially missing, and is now a max
            # selecting only the betacoronaviruses and bat hosts
            if initial_matrix[r,c] == 0 && occursin("Betacoronavirus", viruses[r]) &&  first(df[(df.host_species .== hosts[c]),:]).host_order == "Chiroptera" && output_matrix[r,c] > findmin(maxValues)[1][1]
                # sorting the maximums
                sort!(maxValues)
                # adding the new value at the end of tthe array and overwriting the smallest one
                push!(maxValues, (output_matrix[r,c],r,c))
            end
        end
    end
    # sorting the array in decreasing order
    sort!(maxValues, rev=true)
    # printing the top10 interactions
    println("The top $(top) predicted interactions are:")
    for t in 1:top
        println("Virus: ", viruses[maxValues[t][2]], "  Host: ", hosts[maxValues[t][3]])
        #println("Virus: ", viruses[maxValues[t][2]], "  Host: ", hosts[maxValues[t][3]], " With: ", round(maxValues[t][1]*100, digits=1),"%")
    end
    # calculating the difference between the intial and finale value of the interaction with the highest score
    delta = maxValues[1][1] - initial_matrix[maxValues[1][2],maxValues[1][3]];
    println("The highest scoring interaction is:")
    println("Virus: ", viruses[maxValues[1][2]], "  Host: ", hosts[maxValues[1][3]], " With a Δ of: ", delta)

    # returning the array containing the top10
    return maxValues, delta
end


"""
    calculateInitialeValues(matrix, α)
Calculates the initial values to be attributed by applying a linear filter to the matrix,
using weighted parameters α.

"""
function calculateInitialeValues(matrix, α)
    # converting ones and zeros in boolean
    matrix_bool = convert(Array{Bool}, matrix.== 1)
    # converting matrix in bipartite network
    B = BipartiteNetwork(matrix_bool)
    # applying the linear filter
    F = linearfilter(B, α = α)

    #returning the adjacency matrix
    return(F.A)
end


"""
    generateResultsTable(top10, α, rank, suspected_newData, unlikely_newData, delta, hosts)

Generates an automated results table containing the combinations of alphas
used for the linear filter,the ranks used for the simulation,
and the number of species predicted in the top10 that are also included in the new dataset.
"""

function generateResultsTable(top10, α, rank, suspected_newData, unlikely_newData, delta, hosts)
    # creating a count to keep trace of the suspected and unlikely hosts
    suspected_count = 0;
    unlikely_count = 0;
    # updating the counts for suspected and unlikely hosts
    for i in 1:length(top10)
        if findfirst(suspected -> suspected == hosts[top10[i][3]], suspected_newData) != nothing
            suspected_count += 1;
        end
        if findfirst(unlikely -> unlikely == hosts[top10[i][3]], unlikely_newData) != nothing
            unlikely_count += 1;
        end
    end
    # collecting results
    results = DataFrame(Rank = [rank], Alpha = [α], Suspected = [suspected_count], Unlikely = [unlikely_count], Delta = [delta])
    # wrinting the results in a CSV file
    CSV.write("./figures/Results.csv", results, append = true, delim =';')
end
