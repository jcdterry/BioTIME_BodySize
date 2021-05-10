

### this script containts the source functions to calculate tau


ExtractSpeciesSlope <- function( sp, df){
  df %>%
    filter( TidyBTName == sp)%>%
    lm(  RelRank ~ TimeStd,  data =  .  ) %>% 
    tidy() %>% 
    filter( term == 'TimeStd') %>%
    mutate( TidyBTName = sp) %>%
    select( Slope = estimate, TidyBTName)%>%
    return()
}


FindTrait_Rank_Cor <- function(StudyCell, all_df, trait, response = 'Abundance'){
  cat(StudyCell)
  cat('-')
  
  x<- filter(all_df, STUDY_ID_CELL ==StudyCell)
  
  
  x$TraitToUse <- x[,trait]
  
  ## Extract traits  
  TraitDf <-   distinct(x,  TidyBTName , TraitToUse) 
  
  
  ## Extract start and end years
  StartYear = min(x$YEAR)
  EndYear = max(x$YEAR)
  
  # Calculate relative ranks, by using abundance or biomass as specified
  
  if( response ==  'Abundance'){
    x %>%
      group_by(YEAR) %>%
      summarise(RelRank =  rank(CrossCellAbundance ,  ties.method = "average")/n(), 
                TidyBTName =TidyBTName , .groups= 'keep') -> YYY
  }else{
    x %>%
      group_by(YEAR) %>%
      summarise(RelRank =  rank(CrossCellBiomass ,  ties.method = "average")/n(), 
                TidyBTName =TidyBTName , .groups= 'keep') -> YYY
  }
  
  YYY %>% 
    arrange( YEAR, RelRank) %>%
    spread( YEAR, RelRank, fill = 0) %>%       #  fill in zeros, across all years not observed 
    gather( 'Year',  'RelRank', -TidyBTName  )%>%
    mutate( Year = as.numeric( Year),   # This sets 'time' to be on a 0-1 axis, just to make it a little easier to interpret model fits and intercepts
            TimeStd =(Year-StartYear)/ (EndYear- StartYear)   )-> RelativeRankTimeSeries
  
  ## For each species, calculate its relative rank change through time:
  
  SpeciesSlopes<- map_df( unique(RelativeRankTimeSeries$TidyBTName),  
                          RelativeRankTimeSeries,
                          .f=ExtractSpeciesSlope ) %>%
    mutate( Slope = ifelse(abs(Slope) < 0.00001  , 0 ,Slope    ) )  %>%   ## set any very low values (i.e. flat lines) just to 0
    right_join( TraitDf, by = "TidyBTName")
  
  
  try({  ### known errors: if all traits are the same, then can't get a good correlation out. (most frequent with qualitative size)
    
    OUT <- tibble( trait= trait, 
                   STUDY_ID_CELL =StudyCell, 
                   kendall_cor= cor(SpeciesSlopes$TraitToUse, SpeciesSlopes$Slope, 
                                    method = 'kendall',
                                    use = "na.or.complete")[1])
    return( OUT)
  })
  
  Fail_OUT<- tibble( trait = NA, STUDY_ID_CELL =StudyCell, kendall_cor = NA)
  return( Fail_OUT)
  
}


FindTrait_Rank_Cor_RANDOMISE <- function(StudyCell, all_df, trait, response = 'Abundance', N_REPS=50){
  cat(StudyCell)
  cat('-')
  
  x<- filter(all_df, STUDY_ID_CELL ==StudyCell)
  
  
  x$TraitToUse <- x[,trait]
  
  
  ## Extract start and end years
  StartYear = min(x$YEAR)
  EndYear = max(x$YEAR)
  
  if( response ==  'Abundance'){
    x %>%
      group_by(YEAR) %>%
      summarise(RelRank =  rank(CrossCellAbundance ,  ties.method = "average")/n(), 
                TidyBTName =TidyBTName , .groups= 'keep') -> YYY
  }else{
    x %>%
      group_by(YEAR) %>%
      summarise(RelRank =  rank(CrossCellBiomass ,  ties.method = "average")/n(), 
                TidyBTName =TidyBTName , .groups= 'keep') -> YYY
  }
  
  YYY %>% 
    arrange( YEAR, RelRank) %>%
    spread( YEAR, RelRank, fill = 0) %>%       #  fill in zeros, across all years not observed 
    gather( 'Year',  'RelRank', -TidyBTName  )%>%
    mutate( Year = as.numeric( Year),   # This sets 'time' to be on a 0-1 axis, just to make it a little easier to interpret model fits and intercepts
            TimeStd =(Year-StartYear)/ (EndYear- StartYear)   )-> RelativeRankTimeSeries
  
  ## For each species, calculate its relative rank change through time:
  
  SpeciesSlopes<- map_df( unique(RelativeRankTimeSeries$TidyBTName),  
                          RelativeRankTimeSeries,
                          .f=ExtractSpeciesSlope ) %>%
    mutate( Slope = ifelse(abs(Slope) < 0.00001  , 0 ,Slope    ) )    ## set any very low values (i.e. flat lines) just to 0
  
  ## Randomize Traits and test for correlations, 50 times each 
  TraitDf <-   distinct(x,  TidyBTName , TraitToUse) 
  
  N_reps  = N_REPS
  
  Out<- data.frame(Rep = 1:N_reps,
                   tau_rand = NA, 
                   trait = trait,
                   STUDY_ID_CELL =StudyCell)
  
  for( rep in 1:N_reps){
    try({  ### known errors if all traits are the same, then can't get a good correlation out. (most frequent with qualitative size)
      SpeciesSlopes%>% 
        right_join( TraitDf, by = "TidyBTName") %>%
        mutate (TraitDf,  Trait_Shuffle = sample_frac( TraitToUse)) -> xx
      Out$tau_rand[rep]= cor(xx$Trait_Shuffle, xx$Slope, 
                             method = 'kendall',
                             use = "na.or.complete")[1]
    })
  }
  return( Out)
}

