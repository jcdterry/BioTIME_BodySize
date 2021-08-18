#### All together function#

## this calculates at an individual species level, the full range of values for later use


# 1) the original relative rank based beta
# 2a) the slope through the chopped, scaled data, following dornelas
# 2b) the std error and significance of the fit ^^
# 3a) the slope through the raw data, divided by the mean, following Axel's method
# 3b) the std error and significance of the fit ^^

# also:

# 1) the species trait value + relative trait rank in that assemblage
# 2) the number of observations of the species
# 3) the mean realtive rank of the species through time
# 4) the average absolute relative rank change between each time step. 

## Various tau values can be calculated further down the line, as can the randomisations 

# assemblage level stats:

# study_cell_ID
# trait to be used 


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


Cut_Stdise_BothWays_Fit <- function( sp, df){
  sp_df <- filter(df, TidyBTName == sp)
  
  N_times_observed <- nrow( filter(sp_df, sl_VALUE>0))
  
  ### Cut off extra zeros at either end (i.e. still leave one), to avoid pull towards zero for extinctions or colonisations
  VV <- sp_df$sl_VALUE
  
  
  if(all(VV==0)){  ## very rarely this is occurs, which then buggers up the time series length counting process (e.g. Amphipoda in 176_353543)
    sp.out <- data.frame(TidyBTName = sp,
                         D19_slope = 0,
                         D19_pvalue = 0, # a placeholder values
                         D19_StdErr = 0 ,
                         Elas_slope = 0,
                         Elas_pvalue = 0,
                         Elas_StdErr = 0,
                         N_times_observed = N_times_observed, 
                         N_used = 0)
    
    return( sp.out)
  }
  
  ValuesToUse <- (min( which ( VV != 0 ))-1) : (max( which( VV != 0 ))+1) 
  ## ^^ this could include 0 indicies or be too long hence ->
  ValuesToUse <-ValuesToUse[ ValuesToUse>0 & ValuesToUse <= length(VV)  ]
  
  VV_trimmed <-   VV[ValuesToUse ]
  
  Time <- sp_df$Time[ ValuesToUse]
  
  #####
  # D19 standardisation
  VV_trimmed_D19 <- scale(sqrt(VV_trimmed))
  
  if(is.nan(VV_trimmed_D19[1] )){  ## this  occurs when the counts are completely fixed through time, so there is no way to scale the SD
    ###This is also likely to cause problems with the fitting, so best to jump out early
    sp.out <- data.frame(TidyBTName = sp,
                         D19_slope = 0,
                         D19_pvalue = 0, # a placeholder values
                         D19_StdErr = 0 ,
                         Elas_slope = 0,
                         Elas_pvalue = 0,
                         Elas_StdErr = 0,
                         N_times_observed = N_times_observed, 
                         N_used = length(VV_trimmed))
    
  }else{
    lm(VV_trimmed_D19~Time  ) -> fit_D19
    
    ##############
    # Divide by mean standardisation
    
    VV_trimmed_Elas <- VV_trimmed / mean(VV_trimmed)
    lm(VV_trimmed_Elas~Time  ) -> fit_Elas
    
    sp.out <- data.frame(TidyBTName = sp,
                         D19_slope = coefficients(fit_D19)[2],
                         D19_pvalue = summary(fit_D19)$coefficients[2,4],
                         D19_StdErr = summary(fit_D19)$coefficients[2,2],
                         
                         Elas_slope = coefficients(fit_Elas)[2],
                         Elas_pvalue = summary(fit_Elas)$coefficients[2,4],
                         Elas_StdErr = summary(fit_Elas)$coefficients[2,2],
                         
                         N_times_observed = N_times_observed, 
                         N_used = length(VV_trimmed))
  }
  
  return(sp.out)
}



Find_all_fits_splLvl <- function(StudyCell, all_df, trait, response = 'Abundance'){
  cat(StudyCell)
  cat('-')
  
  x<- filter(all_df, STUDY_ID_CELL ==StudyCell)
  
  x$trait_value <-  unname(unlist(select(x,    !!!trait)))
  
  TraitDf <-   distinct(x,  TidyBTName , trait_value)  ## Extract traits  
  
  
  ## Extract start and end years
  StartYear = min(x$YEAR)
  EndYear = max(x$YEAR)
  
  #############
  # 1. Calculate relative ranks, by using abundance or biomass as specified
  ################
  
  if( response ==  'Abundance'){
    x %>%
      group_by(YEAR) %>%
      summarise(RelRank =  rank(CrossCellAbundance ,  ties.method = "average")/n(), 
                TidyBTName =TidyBTName , .groups= 'keep') -> RanksDF
  }else{
    x %>%
      group_by(YEAR) %>%
      summarise(RelRank =  rank(CrossCellBiomass ,  ties.method = "average")/n(), 
                TidyBTName =TidyBTName , .groups= 'keep') -> RanksDF
  }
  
  RanksDF %>% 
    arrange( YEAR, RelRank) %>%
    spread( YEAR, RelRank, fill = 0) %>%       #  fill in zeros, across all years not observed 
    gather( 'Year',  'RelRank', -TidyBTName  )%>%
    mutate( Year = as.numeric( Year),   # This sets 'time' to be on a 0-1 axis, just to make it a little easier to interpret model fits and intercepts
            TimeStd =(Year-StartYear)/ (EndYear- StartYear)   )-> RelativeRankTimeSeries
  
  ## For each species, calculate its relative rank change through time:
  
  RankSlopes<- map_df( unique(RelativeRankTimeSeries$TidyBTName),  
                       RelativeRankTimeSeries,
                       .f=ExtractSpeciesSlope ) %>%
    mutate( Slope = ifelse(abs(Slope) < 0.00001  , 0 ,Slope    ) ) 
  
  
  
  MeanRankShift <- RelativeRankTimeSeries %>%
    group_by(TidyBTName)%>%
    arrange( TimeStd) %>%
    mutate( Rank_prev_time = lag(RelRank),
            RankChange = RelRank - Rank_prev_time  ) %>%
    summarise(MeanRank = mean(RelRank ),
              MeanAbsRankChange = mean( abs(RankChange), na.rm = TRUE),
              .groups= 'keep')
  
  RankSlopes <- left_join(RankSlopes,MeanRankShift, by = "TidyBTName") 
  
  #########
  ### 2 slopes through quantitative data. 
  
  # identify abundance or biomass to use, then fill in zeros
  
  
  if( response ==  'Abundance'){
    mutate(x, sl_VALUE = CrossCellAbundance ) -> timeSeries
  }else{
    mutate(x, sl_VALUE = CrossCellBiomass ) -> timeSeries
  }
  
  timeSeries %>%
    select( YEAR, sl_VALUE,TidyBTName )%>% 
    arrange( YEAR, sl_VALUE) %>%
    spread( YEAR, sl_VALUE, fill = 0) %>%       #  fill in zeros, across all sample years not observed 
    gather( 'Year',  'sl_VALUE', -TidyBTName  )%>%
    mutate( Time = as.numeric(Year)) -> TimeSeries
  
  FittedLMs_SpeciesResponses<- map_df(unique(TimeSeries$TidyBTName),  
                                      .f=Cut_Stdise_BothWays_Fit, 
                                      df=TimeSeries) 
  
  ######
  ### rate of rank change
  
  ##################
  #### 4.  Calculate relative trait rank, and combine all outputs
  #################
  
  FittedLMs_SpeciesResponses %>%
    left_join( RankSlopes, by = "TidyBTName")%>%
    left_join( TraitDf, by = "TidyBTName")%>%
    mutate( RelTraitRank =  rank(trait_value, 
                                 ties.method = "average", 
                                 na.last = 'keep')/sum(!is.na( trait_value)),
            trait= trait,STUDY_ID_CELL =StudyCell)-> OUT
  
  return(OUT)
}
