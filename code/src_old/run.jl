include("./required.jl")
include("./functions.jl")

# Testing the imputation

# Creating the CSV file for the results
header_list = DataFrame(Rank = [], Alpha = [], Suspected = [], Unlikely = [], Delta = []);
CSV.write("./figures/Results.csv", header_list, delim =';');

# Setting the parameters used for the simulation
α_list = [[0.0, 0.0, 0.0, 1.0], [0.0, 1.0, 1.0, 0.0], [0.0, 1.0, 1.0, 1.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]];
rank_list = [1, 2, 3, 4, 5];

# Creating the list of suspected and unlikely hosts found in the new data
suspected_newData = ["Artibeus jamaicensis","Desmodus rotundus", "Hipposideros larvatus","Hipposideros pomona", "Megaerops kusnotoi", "Myonycteris angolensis", "Myotis pequinius", "Nanonycteris veldkampii", "Nycteris macrotis", "Pipistrellus deserti", "Plecotus auritus", "Pteropus lylei", "Scotophilus heathii", "Scotophilus kuhlii"];
unlikely_newData = ["Carollia sowelli", "Hipposideros gigas", "Hipposideros lekaguli", "Macroglossus minimus", "Myotis horsfieldii", "Pipistrellus coromandra", "Tadarida teniotis"];

# Interating the simulations with all parameters
for i in 1:length(α_list)[1]
    for j in 1:length(rank_list)
        (matrix, hosts, viruses, df) = buildInteractionMatrix();
        crossValidation(0.0, calculateInitialeValues(matrix, α_list[i]), α_list[i], rank_list[j], suspected_newData, unlikely_newData);
    end
end
