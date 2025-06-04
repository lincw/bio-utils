# Load data (if not already loaded)
if (!exists("physical") || !exists("merged")) {
  if (exists("DATA_DIR")) {
    physical <- readRDS(file.path(DATA_DIR, 'results/linkcomm/physical.RDS'))
    merged <- readRDS(file.path(DATA_DIR, 'results/linkcomm/merged_community.RDS'))
  } else {
    stop("DATA_DIR not found. Please set DATA_DIR or load data manually.")
  }
}

# Method 1: Enhanced approach for finding identical and subset communities
find_identical_communities <- function(list1, list2) {
  # Create hash signatures for each community
  create_hash <- function(proteins) {
    paste(sort(proteins), collapse = "|")
  }
  
  # Hash all communities and convert to character sets
  hash1 <- sapply(list1, create_hash)
  hash2 <- sapply(list2, create_hash)
  
  # Convert all communities to sets for subset checking
  set1 <- lapply(list1, function(x) unique(as.character(x)))
  set2 <- lapply(list2, function(x) unique(as.character(x)))
  
  identical_matches <- list()
  
  for(i in seq_along(set1)) {
    set_i <- set1[[i]]
    size_i <- length(set_i)
    
    for(j in seq_along(set2)) {
      set_j <- set2[[j]]
      size_j <- length(set_j)
      
      # Check for perfect match
      if(hash1[i] == hash2[j]) {
        identical_matches[[length(identical_matches) + 1]] <- list(
          list1_id = i,
          list2_id = j,
          jaccard_index = 1.0,
          overlap_type = "perfect",
          list1_size = size_i,
          list2_size = size_j,
          overlap_count = size_i,
          proteins = set_i
        )
      } 
      # Check if set_i is a subset of set_j
      else if(all(set_i %in% set_j)) {
        intersection <- length(intersect(set_i, set_j))
        union_size <- length(union(set_i, set_j))
        jaccard <- intersection / union_size
        
        identical_matches[[length(identical_matches) + 1]] <- list(
          list1_id = i,
          list2_id = j,
          jaccard_index = jaccard,
          overlap_type = "subset",
          list1_size = size_i,
          list2_size = size_j,
          overlap_count = intersection,
          proteins = set_i
        )
      }
      # Check if set_j is a subset of set_i
      else if(all(set_j %in% set_i)) {
        intersection <- length(intersect(set_i, set_j))
        union_size <- length(union(set_i, set_j))
        jaccard <- intersection / union_size
        
        identical_matches[[length(identical_matches) + 1]] <- list(
          list1_id = i,
          list2_id = j,
          jaccard_index = jaccard,
          overlap_type = "superset",
          list1_size = size_i,
          list2_size = size_j,
          overlap_count = intersection,
          proteins = set_j
        )
      }
    }
  }
  
  return(identical_matches)
}

# Method 2: Pre-filtering approach for similar communities
find_similar_communities_fast <- function(list1, list2, 
                                        size_tolerance = 0.5, 
                                        min_jaccard = 0.3) {
  
  # Pre-compute sizes
  sizes1 <- sapply(list1, length)
  sizes2 <- sapply(list2, length)
  
  # Pre-filter by size similarity
  potential_matches <- list()
  comparison_count <- 0
  
  for(i in 1:length(list1)) {
    size1 <- sizes1[i]
    
    # Only compare with communities of similar size
    size_filter <- which(abs(sizes2 - size1) <= size1 * size_tolerance)
    
    for(j in size_filter) {
      comparison_count <- comparison_count + 1
      
      cluster1 <- list1[[i]]
      cluster2 <- list2[[j]]
      
      # Quick overlap check
      intersection_size <- length(intersect(cluster1, cluster2))
      
      # Skip if no overlap
      if(intersection_size == 0) next
      
      # Calculate Jaccard only if there's overlap
      union_size <- length(union(cluster1, cluster2))
      jaccard <- intersection_size / union_size
      
      if(jaccard >= min_jaccard) {
        potential_matches[[length(potential_matches) + 1]] <- list(
          list1_id = i,
          list2_id = j,
          jaccard_index = jaccard,
          overlap_type = "similar",
          overlap_count = intersection_size,
          list1_size = size1,
          list2_size = sizes2[j]
        )
      }
    }
  }
  
  # Sort by Jaccard index
  if(length(potential_matches) > 0) {
    jaccard_scores <- sapply(potential_matches, function(x) x$jaccard_index)
    potential_matches <- potential_matches[order(jaccard_scores, decreasing = TRUE)]
  }
  
  cat("Total comparisons performed:", comparison_count, "\n")
  
  return(potential_matches)
}

# Method 3: Two-step approach combining both methods
find_best_matches_optimized <- function(list1, list2, 
                                      size_tolerance = 0.3, 
                                      min_jaccard = 0.3) {
  
  cat("Step 1: Finding identical communities...\n")
  identical <- find_identical_communities(list1, list2)
  cat("Found", length(identical), "identical communities\n\n")
  
  cat("Step 2: Finding similar communities...\n")
  similar <- find_similar_communities_fast(list1, list2, 
                                         size_tolerance = size_tolerance,
                                         min_jaccard = min_jaccard)
  
  # Remove identical matches from similar matches
  identical_pairs <- sapply(identical, function(x) paste(x$list1_id, x$list2_id, sep = "_"))
  
  if(length(similar) > 0) {
    similar_pairs <- sapply(similar, function(x) paste(x$list1_id, x$list2_id, sep = "_"))
    keep_similar <- !similar_pairs %in% identical_pairs
    similar <- similar[keep_similar]
  }
  
  cat("Found", length(similar), "similar communities\n")
  
  return(list(
    identical = identical,
    similar = similar,
    all_matches = c(identical, similar)
  ))
}

# Create summary for optimized results
create_optimized_summary <- function(optimized_results) {
  all_matches <- optimized_results$all_matches
  
  if(length(all_matches) == 0) {
    return(data.frame(message = "No matches found"))
  }
  
  summary_df <- data.frame(
    match_type = c(rep("Identical", length(optimized_results$identical)),
                   rep("Similar", length(optimized_results$similar))),
    list1_community = sapply(all_matches, function(x) x$list1_id),
    list2_community = sapply(all_matches, function(x) x$list2_id),
    jaccard_index = sapply(all_matches, function(x) x$jaccard_index),
    overlap_count = sapply(all_matches, function(x) {
      if(is.null(x$overlap_count)) x$size else x$overlap_count
    }),
    list1_size = sapply(all_matches, function(x) {
      if(is.null(x$list1_size)) x$size else x$list1_size
    }),
    list2_size = sapply(all_matches, function(x) {
      if(is.null(x$list2_size)) x$size else x$list2_size
    })
  )
  
  return(summary_df[order(summary_df$jaccard_index, decreasing = TRUE), ])
}

# Example usage:
if (!exists("physical") || !exists("merged")) {
  cat("Loading sample data...\n")
  # This is just an example - replace with your actual data loading
  # physical <- readRDS("path/to/physical.RDS")
  # merged <- readRDS("path/to/merged_community.RDS")
}

# Run the comparison
results <- find_best_matches_optimized(physical, merged, 
                                     size_tolerance = 0.3, 
                                     min_jaccard = 0.5)

# Create and display summary
summary_table <- create_optimized_summary(results)
print(head(summary_table, 20))  # Show top 20 matches

# Save results if needed
# write.csv(summary_table, "module_comparison_results.csv", row.names = FALSE)