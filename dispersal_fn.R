#Adult dispersal
dispersal <- function(patch_area, number_patches, S) {
  
  
  # Create the common pool matrix
  D_common_pool <- matrix(rep(patch_area, number_patches), ncol = number_patches, byrow = TRUE) # common pool dispersal such that the movement is proportional to patch relative area
  
  
  require(ramify) # need this for matrix functions (e.g., triu)
  
  if(S > 0){
    tmp1 <- D_common_pool
    tmp1_triu <- triu(tmp1, 1) # set values in the upper triangular part of the matrix (in the 2x2 matrix this mean the cross patch movement)
    tmp1_tril <- tril(tmp1,-1) # set values in the lower triangular part of the matrix (in the 2x2 matrix this mean the cross patch movement)
    tmp1_offdiag <- tmp1_triu + tmp1_tril # Original off diagonal, or cross patch movement (used below)
  } else {
    return(D_common_pool)
  }
  
  for(D_row in 1:length(diag(tmp1))){ # Repeat for each patch, taking the length of the diagonal is equivalent to the total number of patches
    if(tmp1[D_row,D_row] > 0){ # If that site has positive habitat area  
      if((sum(tmp1[D_row, ]) - tmp1[D_row, D_row]) > 0){ # and the cross recruit sites have positive habitat area, this should always be true for our case
        # calculate the difference between larval pool recruitment and 100% self-recruitment
        diff_self_pool <- 1 - tmp1[D_row, D_row] # the 1 being 100% self-recruitment
        # add enhanced site-fidelity to common pool
        enhanced_site_fidelity <- tmp1[D_row, D_row] + (diff_self_pool * S) # tmp1[D_row, D_row] indicates self recruitment to a patch 
        # Standardize the proportion of the off diagonals such that they sum to 1 (see offdiag above)
        tmp1_offdiag_row_standardized <- tmp1_offdiag[D_row, ] / sum(tmp1_offdiag[D_row, ])      
        # Subtract from the common pool cross-recruit probability the standardized diff*proportion enhancement
        reduced_cross_recruit <- tmp1_offdiag[D_row, ] - (tmp1_offdiag_row_standardized * (diff_self_pool*S))
        # Insert the new values
        tmp1[D_row, ] <- reduced_cross_recruit 
        tmp1[D_row, D_row] <- enhanced_site_fidelity
        
      }
    } 
  }
  
  return(tmp1) 
}

##JUVENILE DISPERSAL
dispersal_juvenile <- function(patch_area, number_patches, S_recruit) {
  if(S_recruit > 0){
    tmp2 <- D_common_pool
    tmp2_triu <- triu(tmp2, 1) # set values in the upper triangular part of the matrix (in the 2x2 matrix this mean the cross patch movement)
    tmp2_tril <- tril(tmp2,-1) # set values in the lower triangular part of the matrix (in the 2x2 matrix this mean the cross patch movement)
    tmp2_offdiag <- tmp2_triu + tmp2_tril # Original off diagonal, or cross patch movement (used below)
  } else {
    return(D_common_pool)
  }
  
  for(D_row in 1:length(diag(tmp2))){ # Repeat for each patch, taking the length of the diagonal is equivalent to the total number of patches
    if(tmp2[D_row,D_row] > 0){ # If that site has positive habitat area  
      if((sum(tmp2[D_row, ]) - tmp2[D_row, D_row]) > 0){ # and the cross recruit sites have positive habitat area, this should always be true for our case
        # calculate the difference between larval pool recruitment and 100% self-recruitment
        diff_self_pool <- 1 - tmp2[D_row, D_row] # the 1 being 100% self-recruitment
        # add enhanced site-fidelity to common pool
        enhanced_site_fidelity <- tmp2[D_row, D_row] + (diff_self_pool * S_recruit) # tmp1[D_row, D_row] indicates self recruitment to a patch 
        # Standardize the proportion of the off diagonals such that they sum to 1 (see offdiag above)
        tmp2_offdiag_row_standardized <- tmp2_offdiag[D_row, ] / sum(tmp2_offdiag[D_row, ])      
        # Subtract from the common pool cross-recruit probability the standardized diff*proportion enhancement
        reduced_cross_recruit <- tmp2_offdiag[D_row, ] - (tmp2_offdiag_row_standardized * (diff_self_pool*S_recruit))
        # Insert the new values
        tmp2[D_row, ] <- reduced_cross_recruit 
        tmp2[D_row, D_row] <- enhanced_site_fidelity
        
      }
    } 
  }
  
  return(tmp2) 
}