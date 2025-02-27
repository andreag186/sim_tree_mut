# sim_tree_mut
# **Estimating the Somatic Mutation Rate in Long-Lived Trees: Phylogenomic Approaches and ABC Alternatives**


## **Abstract**
Long-lived trees accumulate mutations throughout their lifetimes, creating intra-individual genetic diversity that persists for centuries. These mutations shape evolutionary potential and species adaptation, yet accurately estimating their rates remains challenging due to the complexity of tree growth and mutation inheritance. In this thesis, I evaluate the phylogenomic method developed by Orr et al. (2020), which links genetic variation to tree topology to estimate mutation rates. However, this method assumes mutations follow the tree's branching structure, an assumption that may not always hold. Using simulations based on Tomimoto and Satake's (2023) models of somatic mutation accumulation, the method is assessed across various tree topologies, mutation rates, and meristem dynamics.

Unbalanced long-terminal topologies performed best under low mutation rates (<10⁻⁹ per site per year), as their extended terminal branches reduce mutation overlap and enhance genetic signal distinctiveness. In contrast, balanced long-terminal topologies systematically underpredicted mutation rates due to diluted mutation signals, while unbalanced short-terminal topologies overpredicted rates due to uneven mutation distributions. Additionally, I observed a substantial increase in shared mutations at higher mutation rates, with the coefficient of variation remaining nearly constant. This suggests that as mutation rates rise, branches share mutations more evenly, reducing phylogenetic distinctiveness. This effect was particularly pronounced in topologies with long shared internal nodes relative to terminal branches, further diminishing the phylogenomic method’s effectiveness in resolving tree topologies.

Given these findings, I recommend using the phylogenomic method for low mutation rates and topologies with minimal shared mutations and distinct genetic signals.

To address its limitations, I present a prototype approach based on Approximate Bayesian Computation (ABC), which does not rely on topological assumptions. The ABC simulation framework approximates the highest posterior density (HPD) for somatic mutation rates and meristem parameters. Validation showed 99.4% of true parameter values were captured within the 95% HPD interval. This prototype offers significant potential as a flexible, scalable alternative for estimating somatic mutation rates, particularly in complex tree topologies. My findings suggest that while the phylogenomic method can be effective under certain conditions, the ABC approach provides a promising direction for future research in tree mutation dynamics.


## **Research Objectives**
The phylogenomic method introduced by Orr et al. (2020) provides a promising framework for estimating somatic mutation rates in long-lived trees by aligning mutation accumulation with tree topology. However, as outlined in Section 1.2, critiques of this method suggest it may underestimate mutation rates by filtering out low-frequency mutations and making assumptions about stem cell replacement dynamics. Additionally, it has not been systematically tested across diverse tree architectures, meristem behaviours, and mutation rate conditions throughout plant taxa.

This thesis addresses these gaps by employing Tomimoto and Satake’s (2023) hierarchical modular growth model to simulate somatic mutation accumulation in a tree under controlled conditions, providing a framework to assess and refine somatic mutation rate estimation. The research is structured around two core objectives:

1. **Testing the Phylogenomic Method**: I aim to evaluate the accuracy and limitations of the methodology as described by Orr et al. (2020) by investigating two key factors through simulation.

   - **Mutation Rate Recovery**: The accuracy of the phylogenomic method in recovering inputted mutation rates across varying tree topologies and other biological parameters is assessed using simulated mutations.

   - **Topology Recovery**: A phylogenetic tree is reconstructed from simulated somatic mutations across a tree and compared to the true topology of the tree.

2. **Developing a Novel Approximate Bayesian Computation (ABC) Method**: I introduce a proof-of-concept, non-phylogenetic approach for estimating the somatic mutation rates of long-lived trees utilising Tomimoto and Satake’s (2023) simulation framework and an ABC-Reject methodology. This novel method does not assume that mutations strictly follow topology, thus avoiding over-filtering. Instead, the ABC-Reject framework approximates a plausible interval of somatic mutation rates which produce observed mutation frequencies across branches, as well as approximating parameters of interest of meristem behaviour.

## *See the attached pdf for the full Masters Thesis relating to this github repository*
