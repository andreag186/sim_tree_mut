# FULL NEWICK!!!!
#CALC RF DISTANCE FOR SIMULATION TREES

#PREP
library(ape)
library(phangorn)
library(tidyverse)
library(reticulate)

#PYTHON RETICULATE SETUP
use_python("C:/Users/amgst/AppData/Local/Programs/Python/Python38/python.exe", required = TRUE)  
py_config()

py_run_string("
import ast

def parse_unique_mutations(mutations_str):
    try:
        mutations_list = ast.literal_eval(mutations_str)
        return mutations_list
    except (ValueError, SyntaxError) as e:
        print(f'Error parsing Unique_Mutations: {mutations_str}')
        return []

def parse_shared_mutations(mutations_str):
    try:
        shared_mutations = ast.literal_eval(mutations_str)
        shared_mutations = {str(k): v for k, v in shared_mutations.items()}
        return shared_mutations
    except (ValueError, SyntaxError) as e:
        print(f'Error parsing Shared_Mutations: {mutations_str}')
        return {}
")

#UPLOAD DATA
df <- read.csv("new_F.csv")
#SAVE VERSION
filtered_df_saved <- df

#PARSING FIX R/PYTHON INCOMPATIBILITY
df$Unique_Mutations <- lapply(df$Unique_Mutations, function(x) {
  py$parse_unique_mutations(as.character(x))
})

df$Shared_Mutations <- lapply(df$Shared_Mutations, function(x) {
  py$parse_shared_mutations(as.character(x))
})
#CHECK ERROR
any_na_unique <- sapply(df$Unique_Mutations, function(x) any(is.na(x)))
any_na_shared <- sapply(df$Shared_Mutations, function(x) any(sapply(x, is.null)))


if (any(any_na_unique)) {
  warning("NAs introduced by coercion in Unique_Mutations")
}
if (any(any_na_shared)) {
  warning("NAs introduced by coercion in Shared_Mutations")
}


#MAKE PHYLO FROM MUT DIST
create_phylo_tree <- function(unique_mutations, shared_mutations, sample_number) {
  if (is.null(unique_mutations) || is.null(shared_mutations)) {
    return(NULL)
  }
  
  terminals <- length(unique_mutations)
  dist_matrix <- matrix(0, nrow = terminals, ncol = terminals)
  
  # dist matrix from shared mut
  for (shared_branches in names(shared_mutations)) {
    shared_value <- shared_mutations[[shared_branches]]
    branches <- as.numeric(unlist(strsplit(gsub("[()]", "", shared_branches), ",")))
    
    if (any(is.na(branches)) || is.na(shared_value)) {
      cat("Parsing issue with shared_mutations:", shared_branches, "->", shared_value, "\n")
      next
    }
    
    for (j in 1:(length(branches) - 1)) {
      for (k in (j + 1):length(branches)) {
        branch1 <- branches[j]
        branch2 <- branches[k]
        if (!is.na(branch1) && !is.na(branch2) && branch1 <= terminals && branch2 <= terminals) {
          dist_matrix[branch1, branch2] <- shared_value
          dist_matrix[branch2, branch1] <- shared_value
        }
      }
    }
  }
  
  # DIST MATRIX ADJUST
  for (i in 1:terminals) {
    for (j in 1:terminals) {
      if (i != j) {
        dist_matrix[i, j] <- unique_mutations[[i]] + unique_mutations[[j]]
        
        for (shared_branches in names(shared_mutations)) {
          branches <- as.numeric(unlist(strsplit(gsub("[()]", "", shared_branches), ",")))
          if (i %in% branches && !(j %in% branches)) {
            dist_matrix[i, j] <- dist_matrix[i, j] + shared_mutations[[shared_branches]]
          }
          if (j %in% branches && !(i %in% branches)) {
            dist_matrix[i, j] <- dist_matrix[i, j] + shared_mutations[[shared_branches]]
          }
        }
      }
    }
  }
  
  
  
  tree <- nj(as.dist(dist_matrix))
  tree$tip.label <- as.character(1:terminals)
  
  return(tree)
}


#RF DISTS DATA FRAME
rf_distances <- data.frame(
  Sample_Number = df$Sample_Number,
  tree_topology = df$tree_topology,
  RF_Distance = NA
)

#STANDARDIZE TIP LABELS
standardize_labels <- function(tree1, tree2) {
  common_labels <- intersect(tree1$tip.label, tree2$tip.label)
  tree1 <- drop.tip(tree1, setdiff(tree1$tip.label, common_labels))
  tree2 <- drop.tip(tree2, setdiff(tree2$tip.label, common_labels))
  return(list(tree1, tree2))
}

#STANDARDIZE NEWICK
standardize_newick <- function(newick, terminals) {
  tree <- read.tree(text = newick)
  tree$tip.label <- as.character(1:terminals)
  return(tree)
}


#FULL NEWICK DICT
newick_dict <- list(
  test = "((((Lt3:70,Lt2:64)Li2:8,Lt1:80)Li1:77,Lt0:84)Li0:36,(((Rt3:98,Rt2:93)Ri2:11,Rt1:103)Ri1:11,Rt0:130)Ri0:36)s10:1;",
  bS4 = "(((Rt1:38,Rt0:38)Ri0:75),(Lt1:38,Lt0:38)Li0:75);", 
  ubS4 = "(((Rt2:38,Rt1:38)Ri1:75,Rt0:38)Ri0:75,Lt0:38);", 
  bL4 = "(((Rt1:60,Rt0:60)Ri0:30),(Lt1:60,Lt0:60)Li0:30);", 
  ubL4 = "(((Rt2:60,Rt1:60)Ri1:30,Rt0:60)Ri0:30,Lt0:60);",
  bS6 = "((((Rt2:21,Rt1:21)Ri1:42,Rt0:21)Ri0:42),((Lt2:21,Lt1:21)Li1:42,Lt0:21)Li0:42);", 
  ubS6 = "(((((Rt4:21,Rt3:21)Ri3:42,Rt2:21)Ri2:42,Rt1:21)Ri1:42,Rt0:21)Ri0:42,Lt0:21);",
  bL6 = "((((Rt2:38,Rt1:38)Ri1:19,Rt0:38)Ri0:19),((Lt2:38,Lt1:38)Li1:19,Lt0:38)Li0:19);", 
  ubL6 = "(((((Rt4:38,Rt3:38)Ri3:19,Rt2:38)Ri2:19,Rt1:38)Ri1:19,Rt0:38)Ri0:19,Lt0:38);",
  bS8 = "(((((Rt3:15,Rt2:15)Ri2:30,Rt1:15)Ri1:30,Rt0:15)Ri0:30),(((Lt3:15,Lt2:15)Li2:30,Lt1:15)Li1:30,Lt0:15)Li0:30);",
  ubS8 = "(((((((Rt6:15,Rt5:15)Ri5:30,Rt4:15)Ri4:30,Rt3:15)Ri3:30,Rt2:15)Ri2:30,Rt1:15)Ri1:30,Rt0:15)Ri0:30,Lt0:15);",
  bL8 = "(((((Rt3:26,Rt2:26)Ri2:13,Rt1:26)Ri1:13,Rt0:26)Ri0:13),(((Lt3:26,Lt2:26)Li2:13,Lt1:26)Li1:13,Lt0:26)Li0:13);", 
  ubL8 = "(((((((Rt6:26,Rt5:26)Ri5:13,Rt4:26)Ri4:13,Rt3:26)Ri3:13,Rt2:26)Ri2:13,Rt1:26)Ri1:13,Rt0:26)Ri0:13,Lt0:26);", 
  bS10 = "((((((Rt4:12,Rt3:12)Ri3:24,Rt2:12)Ri2:24,Rt1:12)Ri1:24,Rt0:12)Ri0:24),((((Lt4:12,Lt3:12)Li3:24,Lt2:12)Li2:24,Lt1:12)Li1:24,Lt0:12)Li0:24);",
  ubS10 = "(Lt0:24,(Rt0:12,(Rt1:12,(Rt2:12,(Rt3:12,(Rt4:12,(Rt5:12,(Rt6:12,(Rt8:12,Rt7:12)Ri7:24)Ri6:24)Ri5:24)Ri4:24)Ri3:24)Ri2:24)Ri1:24)Ri0:24);",
  bL10 = "((((((Rt4:22,Rt3:22)Ri3:11,Rt2:22)Ri2:11,Rt1:22)Ri1:11,Rt0:22)Ri0:11),((((Lt4:22,Lt3:22)Li3:11,Lt2:22)Li2:11,Lt1:22)Li1:11,Lt0:22)Li0:11);", 
  ubL10 = "(Lt0:22,(Rt0:22,(Rt1:22,(Rt2:22,(Rt3:22,(Rt4:22,(Rt5:22,(Rt6:22,(Rt8:22,Rt7:22)Ri7:11)Ri6:11)Ri5:11)Ri4:11)Ri3:11)Ri2:11)Ri1:11)Ri0:11);",
  bS12 = "(((((((Rt5:9,Rt4:9)Ri4:18,Rt3:9)Ri3:18,Rt2:9)Ri2:18,Rt1:9)Ri1:18,Rt0:9)Ri0:18),(((((Lt5:9,Lt4:9)Li4:18,Lt3:9)Li3:18,Lt2:9)Li2:18,Lt1:9)Li1:18,Lt0:9)Li0:18);", 
  ubS12 = "(((((((((((Rt10:9,Rt9:9)Ri9:18,Rt8:9)Ri8:18,Rt7:9)Ri7:18,Rt6:9)Ri6:18,Rt5:9)Ri5:18,Rt4:9)Ri4:18,Rt3:9)Ri3:18,Rt2:9)Ri2:18,Rt1:9)Ri1:18,Rt0:9)Ri0:18,Lt0:9);",
  bL12 = "(((((((Rt5:16,Rt4:16)Ri4:8,Rt3:16)Ri3:8,Rt2:16)Ri2:8,Rt1:16)Ri1:8,Rt0:16)Ri0:8),(((((Lt5:16,Lt4:16)Li4:8,Lt3:16)Li3:8,Lt2:16)Li2:8,Lt1:16)Li1:8,Lt0:16)Li0:8);",
  ubL12 = "(((((((((((Rt10:16,Rt9:16)Ri9:8,Rt8:16)Ri8:8,Rt7:16)Ri7:8,Rt6:16)Ri6:8,Rt5:16)Ri5:8,Rt4:16)Ri4:8,Rt3:16)Ri3:8,Rt2:16)Ri2:8,Rt1:16)Ri1:8,Rt0:16)Ri0:8,Lt0:16);"
)

# LOOP AND CALC RF
for (i in seq_len(nrow(df))) {
  unique_mut <- df$Unique_Mutations[[i]]
  shared_mut <- df$Shared_Mutations[[i]]
  tree <- create_phylo_tree(unique_mut, shared_mut, i)
  
  if (!is.null(tree)) {
    newick_tree <- standardize_newick(newick_dict[[df$tree_topology[i]]], length(tree$tip.label))
    
    trees <- standardize_labels(tree, newick_tree)
    tree <- trees[[1]]
    newick_tree <- trees[[2]]
    
    if (i <= 10) {
      cat("Sample:", i, "\n")
      cat("Tree tip labels:", tree$tip.label, "\n")
      cat("Newick tree tip labels:", newick_tree$tip.label, "\n")
      cat("Shared Mutations for Sample", i, ":\n")
      print(shared_mut)
      
      plot(tree, main = paste("Generated Tree - Sample", i))
      plot(newick_tree, main = paste("Reference Tree - Sample", i))
    }
    
    rf_distances$RF_Distance[i] <- RF.dist(tree, newick_tree)
  } else {
    rf_distances$RF_Distance[i] <- NA
  }
}

print(head(rf_distances))
#check rows
print(nrow(rf_distances))

# merge rf and og df
merged_df <- filtered_df_saved %>%
  left_join(rf_distances, by = c("Sample_Number", "tree_topology"), relationship = "many-to-many")


# check
print(nrow(merged_df))

#save
write.csv(merged_df, "nj_full.csv", row.names = FALSE)

# ANOVA + PLOT
variables <- c("Model", "tree_topology", "biasVar", "StD")

anova_results <- lapply(variables, function(var) {
  formula <- as.formula(paste("RF_Distance ~", var))
  aov_model <- aov(formula, data = merged_df)
  summary(aov_model)
})

names(anova_results) <- variables

for (var in variables) {
  cat("ANOVA results for", var, ":\n")
  print(anova_results[[var]])
  cat("\n")
}

extract_anova_table <- function(aov_summary) {
  df <- as.data.frame(aov_summary[[1]])
  rownames_to_column(df, var = "Variable")
}

summary_table <- do.call(rbind, lapply(anova_results, extract_anova_table))
print(summary_table)

merged_df <- merged_df %>%
  mutate(Balance = ifelse(startsWith(tree_topology, "b"), "Balanced", "Unbalanced"))

ggplot(merged_df, aes(x = Balance, y = RF_Distance, fill = Balance)) +
  geom_violin() +
  labs(title = "RF Distances of Balanced and Unbalanced Trees",
       x = "Tree Type",
       y = "RF Distance") +
  theme_minimal() +
  theme(legend.position = "none")