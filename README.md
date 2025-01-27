# sim_tree_mut
# **Estimating the Somatic Mutation Rate in Long-Lived Trees: Phylogenomic Approaches and ABC Alternatives**


## **Simulating Somatic Mutation Accumulation**
In response to the computation inefficiencies of the original pipeline (Orr et al., 2020), we developed a streamlined ‘approximate phylogenomic method’ optimised for high-throughput simulations.
A Python script denoted ‘sim_code.py’ is used to execute the simulations and can be found within the open GitHub repository sim_tree_mut (https://github.com/andreag186/sim_tree_mut). This script was written in Python and is compatible with version 3.11.4. The simulation code comprises several functions which handle various stages of the pipeline, from parameter sampling to simulation of somatic mutations and application of the approximate phylogenomic method, and finally, the estimation of a somatic mutation rate. Here, each function’s role is explained in the order in which they are executed. 

Latin Hypercube Sampling (LHS)  of Parameters (latin_hypercube_sampling)
LHS generates diverse parameter sets for each sample to explore parameter space comprehensively. LHS is implemented using the qmc.LatinHypercube class from scipy. The variables which form a parameter set are listed below.
mu_0: Mutation rate parameter, selected from 6 treatment levels specified:
[4, 6, 8, 40, 60, 80] * 1-10
GenSize: Genome size, fixed for each sample as proportional to mu_0 
tree_topologies: A list of samples, through a dictionary of predefined tree topologies and input tree parameters, of 20 treatment levels specified by a unique code.*
StD: Elongation parameter, 0 = Structured and 5=Stochastic (equal to stem cells)
biasVar: Branching parameter, 10 = Unbiased and 0.5 = Biassed. 
num_params : The number of parameters used in LHS (= 4, with GenSize constant). 
num_samples: The number of samples specified when the function is called
	
	* tree_topologies_dict is the dictionary that defines each input tree's number of branches, age, and topology. See the example entry of tree ‘bS4’ below:
 “bS4”: {
	“numBranch”:4,  -The number of terminal branches
	“age”: 123, -The root to tip ‘age’, not total age(= 300) , equivalent to height in cm
	“s10”: 10,  - The age of the trunk/root prior to the first split, 10 years for each tree (1m)
	“b11”: 75, “bb11”: 38, “b12”: 38 -Branches RIGHT of split
	“s40”: 75, “b41”:38, “s41”: 38   -Branches LEFT of split
}
Here,  “bS4” indicates this topology is balanced (b) with short terminal branches (S) and 4 terminal branches. For branches RIGHT of the split, internal branches are denoted by b11-bxx, with the last ‘internal branch’ treated as a terminal branch. Terminal branches are denoted by bb11-bbxx, with the number corresponding to the internal branch on which they stem. For branches LEFT of the split, internal branches are denoted by s40-sxx, and terminal branches by b41-bxx. This nomenclature structure is preserved from the original Tomimoto and Satake (2023) source code. See Figure 1, Appendix. 


Result Compilation and Output to CSV (write_samples_to_csv)
This function writes all simulation results to a CSV file using csv.writer and pandas for handling tabular data. The sample number, input parameter values, and all output values are recorded. 


‘Decode’ Tree Topology Configuration (create_tree_list_and_dict):
This function ‘decodes’ the input trees defined in the tree_topologies_dict, extracting all information needed for future functions. This function allows ANY input tree with any number of branches and topology to be compatible with the simulation, given that the tree is correctly structured in tree_topologies_dict. This function returns the following:
tree_list: An output list of all branch names and corresponding age (in years)
tree_dict: An output dictionary classifying branches into right or left and terminal and internal branches for future applications. 
numBranch and age are also outputted here for subsequent calculations.

Average Meristem Time Calculation (calculate_ave_meristem_time)
This function calculates the average time (in years) spent by meristem cells in each branch, using root-to-tip distances of left and right branches. This value is used later to estimate mutation accumulation over time (in mut_dist_func )
tree_dict is inputted specifying the root, right_internal, right_terminal, left_internal and left_terminal values for each input tree. 
ave_meristem_time: The average meristem time, calculated as the sum of root-to-tip distances divided by numBranch


Main Simulation Loop (simulate_somatic_mutations)
This function iterates through each specified replicate (NumTime), generating mutations in unique genome positions across branches based on input parameters (e.g. model). A combined matrix (br_brmutmatrix) is created (using the package numpy) containing all the mutations accumulated in the tree and their unique site patterns, specifying in which branch or branches each mutation occurs. * 
NumStem (number of stem cells = 5), NumTime (number of simulation replications = 1), mu_0 (input mutation rate), GenSize (genome size = 500Mb), NumDiv (branching division constant = 7), nDiv (elongation division constant = 1), StD and biasVar (model parameters): The core input parameters which define mutation simulation. 
List10, List_br, and List_Left: The data structures storing mutations for the root (growth prior to the first split), right branches, and left branches, respectively.
	
	* The function simulate_somatic_mutations, calls upon functions 6-9 defined by Tomimoto and Satake (2023) in their source code to simulate somatic mutations by applying their mathematical models. We did not alter these functions, only their application.

Stem Cell Mutation Simulation (mutInStemCells)
This function initialises a list of stem cells and simulates mutations across consecutive divisions; in other words, it simulates elongation in the ‘stem’ of the tree prior to the first split of branches. copy.deepcopy ensures each stem cell's mutation history is preserved accurately, with random.random() introducing mutations based on the probability mu_0 (mutation rate).  The mutation history of each stem cell is saved in the object tCells and updated after each division cycle.

Branch Mutation Simulation (mutInBrStemCells)
Similar to the function mutInStemCells, this function models mutation accumulation in branch-specific cells (axillary meristems formed from the SAM). In other words, it simulates elongation along a branch. The stemCells list is passed as input. 
stemCells: A list of the initial set of stem cells with their mutation history, as sampled from the SAM
ccells and ccells2: Temporary lists for the daughter cells produced in each division cycle (during axillary meristem formation, aka. branching AND elongation along a branch)
tCells: Final list of all mutated cells for the given branch.

Sampling and Mutation Assignment (sample_mutations)
This function samples cells for the next branch generation using weights defined by weightList, which is based on spatial (x [0,2) ) and bias (  or biasVar) parameters. In other words, it simulates the branching process. numpy.random.choice selects cells for each division.
wList: A list of weights generated to reflect spatial positioning and bias
NumDiv: Number of cell divisions for each branching event ( rb = 7)
sample: List of sampled cells for the new meristem, used as input for further elongation simulations along that branch. 

Mutation Matrix Generation (makeMutMatrix)
This function generates a binary mutation matrix, indicating the presence or absence of across stem cells in each branch. The matrix provides a detailed mutation map.
List: An input list of meristems with mutation histories
mutList: A list of all unique mutations used to populate the matrix
np_per_mutMatrix: The final mutation matrix with binary entries, representing mutation presence in each cell of one branch. The final br_brmutMatrix output for the whole tree collates and filters np_per_mutMatrix objects for each branch.

Mutation Frequency and Distribution Analysis (mut_dist_func)
This function computes mutation distribution patterns, calculating each branch's total, shared and unique mutations. It summarises the mutation frequencies and patterns, capturing distribution across branches.
mutShapeTemplate: A template of all possible mutation combinations across branches (dependent on numBranch)
allMutDist: An array of mutation distribution patterns averaged across simulation runs (as defined by NumTime)
total_mutations, shared_mutations and unique_mutations: Dictionaries storing counts of specified mutations corresponding relevant branch number(s).

Distance Matrix Construction (gen_matrices)
This function constructs genetic and physical distance matrices to analyse each tree's mutational and spatial distances. The genetic matrix (genetic_distance_matrix) captures shared mutations between branches, while the physical matrix (physical_distance_matrix) represents distances (in cm) between branches based on root-to-tip calculations (root_to_tip_distances). 

Mutation Rate and Variant Calculations (calc_variants)
This function performs a linear regression using sklearn.linear_model.LinearRegression. This function estimates the relationship between genetic and physical distances within a tree, thereby estimating a somatic mutation rate using our approximate phylogenomic method. The function outputs the estimated mutation rate as well as the regression equation. 
genetic_distances and physical_distances: Arrays of distances extracted from the upper triangle of matrices, used as inputs for the regression model. 
regressor: A linear regression model fits the genetic and physical distances.
output_mut_rate: The estimated mutation rate per nucleotide, derived directly from the regressor. This value is later converted to per site per year by dividing by age (300 years). 

