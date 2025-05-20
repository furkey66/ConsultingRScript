
library(haven)         # reading in SPSS files
library(skimr)         # overview data
library(tidyverse)     # data wrangling, manipulation, piping
library(naniar)        # missing data analysis
library(mice)          # multiple imputation
library(magrittr)      # combines piping and assignment ( %<>% )
library(caret)         # useful for data-preparation 
library(pROC)          # AUC values and ROC plots
library(datawizard)    # standardization of variables
library(glmnet)        # GLM with LASSO and cross-validation

#------------------------------------------------------------------------------#
####------------------------------Data Wrangling----------------------------####
#------------------------------------------------------------------------------#

# Read data
DF_full <- read_sav ( "20250313_Df_labelled.sav" )

DF_complete <- DF_full %>%
  # Step 1: remove observations of which the values are improbable/impossible.
  filter ( PPI != 33 | is.na(PPI), 
           BMI < 200 | is.na(BMI), 
           BMI > 0 | is.na(BMI),
           wijzigingVasoactief != 2 | is.na(wijzigingVasoactief),
           ZP != 300) %>%
  # Step 2: make sex dichotomous and transform values of -98 and -99 into NA. 
  mutate ( sex = sex - 1,
           across ( everything ( ), ~ ifelse ( . %in% c ( -98, -99 ), NA, . ) ) ) %>% 
  # Step 3: give proper names for the ID and Visit variables.
  rename ( ID = ZP,
           Visit = Index1 ) %>% 
  mutate ( VisitDate = as.Date ( VisitDate, format = "%Y-%m-%d" ) ) %>% 
  # Step 4: remove observations with missing values on outcome variables
  filter ( !is.na ( wijzigingIS ) ) %>% 
  # Step 5: correct visit numbers
  mutate ( Visit = case_when ( row_number ( ) %in% 979:989 ~ Visit - 9,         # ID 167
                               TRUE ~ Visit ),
           Visit = case_when ( row_number ( ) %in% 2949:2951 ~ Visit - 4,       # ID 790
                               TRUE ~ Visit ),
           Visit = case_when ( row_number ( ) %in% 3127:3128 ~ Visit - 1,       # ID 936
                               TRUE ~ Visit ),
           Visit = case_when ( row_number ( ) %in% 3142 ~ Visit - 2,            # ID 954
                               TRUE ~ Visit ) )

# Recode lab parameters such that values are 0 when they fall inside reference
# range, and not 0 where the value represents the deviation or distance from
# this reference with its corresponding upper and lower bound. 
DF_ranged <- DF_complete %>%
  mutate (
    BSE = ifelse ( BSE <= 20 & BSE >= 0, 0, 
                   ifelse ( BSE < 0, 
                            0 - BSE, 
                            BSE - 20 ) ),
    Hemoglobine = ifelse ( Hemoglobine <= 10 & Hemoglobine >= 7.5, 0,
                           ifelse ( Hemoglobine < 7.5, 
                                    7.5 - Hemoglobine, 
                                    Hemoglobine - 10 ) ),
    Hematocriet = ifelse ( Hematocriet <= 0.45 & Hematocriet >= 0.35, 0,
                           ifelse ( Hematocriet < 0.35, 
                                    0.35 - Hematocriet, 
                                    Hematocriet - 0.45 ) ),
    Erytrocyten = ifelse ( Erytrocyten <= 5 & Erytrocyten >= 4, 0,
                           ifelse ( Erytrocyten < 4, 
                                    0 - Erytrocyten, 
                                    Erytrocyten - 0 ) ),
    Erytroblasten = ifelse ( Erytroblasten < 0.1, 0, 
                             Erytroblasten - 0.1 ),
    MCV = ifelse ( MCV <= 100 & MCV >= 80, 0,
                   ifelse ( MCV < 80, 
                            80 - MCV, 
                            MCV - 100 ) ),
    MCH = ifelse ( MCH <= 2.1 & MCH >= 1.7, 0,
                   ifelse ( MCH < 1.7, 
                            1.7 - MCH, 
                            MCH - 2.1 ) ),
    MCHC = ifelse ( MCHC <= 16 & MCHC >= 12, 0,
                    ifelse ( MCHC < 12, 
                             12 - MCHC, 
                             MCHC - 16 ) ),
    Trombocyten = ifelse ( Trombocyten <= 400 & Trombocyten >= 150, 0,
                           ifelse ( Trombocyten < 150, 
                                    150 - Trombocyten, 
                                    Trombocyten - 400 ) ),
    Leukocyten = ifelse ( Leukocyten <= 10 & Leukocyten >= 4, 0,
                          ifelse ( Leukocyten < 4, 
                                   4 - Leukocyten, 
                                   Leukocyten - 10 ) ),
    Natrium = ifelse ( Natrium <= 145 & Natrium >= 135, 0,
                       ifelse ( Natrium < 135, 
                                135 - Natrium, 
                                Natrium - 145 ) ),
    Kalium = ifelse ( Kalium <= 5.1 & Kalium >= 3.5, 0,
                      ifelse ( Kalium < 3.5, 
                               3.5 - Kalium, 
                               Kalium - 5.1 ) ),
    Calcium = ifelse ( Calcium <= 2.55 & Calcium >= 2.15, 0,
                       ifelse ( Calcium < 2.15, 
                                2.15 - Calcium, 
                                2.55 - Calcium ) ),
    Calcium_gecorrigeerd = ifelse ( Calcium_gecorrigeerd <= 2.55 & Calcium_gecorrigeerd >= 2.15, 0,
                                    ifelse ( Calcium_gecorrigeerd < 2.15, 
                                             2.15 - Calcium_gecorrigeerd, 
                                             2.55 - Calcium_gecorrigeerd ) ),
    Fosfaat = ifelse ( Fosfaat <= 1.5 & Fosfaat >= 0.9, 0,
                       ifelse ( Fosfaat < 0.9, 
                                0.9 - Fosfaat, 
                                Fosfaat - 1.5 ) ),
    Kreatinine = ifelse ( Kreatinine <= 90 & Kreatinine >= 49, 0,
                          ifelse ( Kreatinine < 49, 
                                   49 - Kreatinine, 
                                   Kreatinine - 90 ) ),
    eGFR = ifelse ( eGFR > 60, 0, 60 - eGFR ),
    Urinezuur = ifelse ( Urinezuur <= 0.34 & Urinezuur >= 0.14, 0,
                         ifelse ( Urinezuur < 0.14, 
                                  0.14 - Urinezuur, 
                                  Urinezuur - 0.34 ) ),
    Albumine = ifelse ( Albumine <= 48 & Albumine >= 34, 0,
                        ifelse ( Albumine < 34, 
                                 34 - Albumine, 
                                 Albumine - 48 ) ),
    LDH = ifelse ( LDH < 247, 0, LDH - 247 ),
    ASAT = ifelse ( ASAT < 31, 0, ASAT - 31 ),
    ALAT = ifelse ( ALAT < 34, 0, ALAT - 34 ),
    AlkalischFosfatase = ifelse ( AlkalischFosfatase < 98, 0, AlkalischFosfatase - 98 ),
    GammaGT = ifelse ( GammaGT <= 40 & GammaGT >= 5, 0, GammaGT - 40 ),
    CRP = ifelse ( CRP < 5, 0, CRP - 5 ),
    CK = ifelse ( CK < 145, 0, CK - 145 ),
    TropoT = ifelse ( TropoT < 14, 0, TropoT - 14 ),
    NTProBNP = ifelse ( NTProBNP < 247, 0, NTProBNP - 247 ),
    Cholesterol = ifelse ( Cholesterol < 5, 0, Cholesterol - 5 ),
    HDL = ifelse ( HDL > 1, 0, 1 - HDL ),
    LDL = ifelse ( LDL < 3, 0, LDL - 3 ),
    CholHDLCratio = ifelse ( CholHDLCratio < 5, 0, CholHDLCratio - 5 ),
    Ijzer = ifelse ( Ijzer <= 25 & Ijzer >= 10, 0, 
                     ifelse ( Ijzer < 10, 
                              10 - Ijzer, 
                              Ijzer - 25 ) ),
    Triglyceriden = ifelse ( Triglyceriden < 2, 0, Triglyceriden - 2 ),
    Ferritine = ifelse ( Ferritine <= 150 & Ferritine >= 10, 0,
                         ifelse ( Ferritine < 10, 
                                  10 - Ferritine, 
                                  Ferritine - 150 ) ),
    Transferrine = ifelse ( Transferrine <= 3.6 & Transferrine >= 2, 0,
                            ifelse ( Transferrine < 2, 
                                     2 - Transferrine, 
                                     Transferrine - 3.6 ) ),
    Ijzerverzadiging = ifelse ( Ijzerverzadiging <= 45 & Ijzerverzadiging >= 16, 0,
                                ifelse ( Ijzerverzadiging < 16, 
                                         16 - Ijzerverzadiging, 
                                         Ijzerverzadiging - 45 ) ),
    TSH = ifelse ( TSH <= 4.8 & TSH >= 0.3, 0,
                   ifelse ( TSH < 0.3, 
                            0.3 - TSH, 
                            TSH - 4.8 ) ),
    VrijT4 = ifelse ( VrijT4 <= 25 & VrijT4 >= 12, 0,
                      ifelse ( VrijT4 < 12, 
                               12 - VrijT4, 
                               VrijT4 - 25 ) ) ) %>%
  select ( -Erytroblasten )                                                     # Removing erytroblasten

# Specifying groups of variables
lab_params        <- names ( DF_complete [, c ( 6:9, 11:44 ) ] )                # All lab parameters ( - erytroblasten )
med_info_params   <- names ( DF_complete [ , 45:53 ] )                          # Medication changes
med_intake_params <- names ( DF_complete [ , 54:68 ] )                          # Medication use
med_Y_params      <- c ( "wijzigingIS",                                         # Outcome variables
                         "wijzigingVasoactief", 
                         "wijzingSupportmed" )
cov_params        <- names ( DF_complete [ , c ( 4:5, 69:79 ) ] )               # Covariates = Disease Manifestations + Sex + YearofBirth + BMI


# Relevant variables for analysis only
DF_relevant <- DF_ranged %>% 
  select ( ID, Visit, lab_params, med_Y_params, cov_params )

################################################################################
####--------------------------Model Validation------------------------------####
################################################################################

# Set seed for reproducibility 
set.seed ( 66 )

# Make a train-test split based on observations
train_ID <- sample ( seq_len ( nrow ( DF_relevant ) ), 
                     size = ceiling ( 0.7 * nrow ( DF_relevant ) ) )
DF_train <- DF_relevant [ train_ID, ]
DF_test  <- DF_relevant [ -train_ID, ]

all_predictors <- c ( lab_params, cov_params, "Visit" )                         # Shouldn't we also include visit number?

#------------------------------Imputation Data---------------------------------#

DF_train_imp <- mice::complete ( mice ( DF_train, m = 1, maxit = 5,             # Training data
                                        print = FALSE ), 1 )
DF_test_imp <- mice::complete ( mice ( DF_test, m = 1, maxit = 5,               # Test data
                                       print = FALSE ), 1 )
DF_full_imp <- rbind ( DF_train_imp, DF_test_imp )                              # Full data

#----------------------Define predictors and outcomes--------------------------#

# Training Data
x_train <- model.matrix ( as.formula ( paste ( "~ ", paste ( all_predictors, collapse = " + " ) ) ), 
                          data = DF_train_imp) [ , -1 ]
x_train_lab <- model.matrix ( as.formula ( paste ( "~ ", paste ( all_predictors, collapse = " + " ) ) ), 
                              data = DF_train_imp) [ , 2 : 39  ]
y_train_IS <- DF_train_imp$wijzigingIS
y_train_Vaso <- DF_train_imp$wijzigingVasoactief
y_train_Support <- DF_train_imp$wijzingSupportmed

# Test Data
x_test <- model.matrix ( as.formula ( paste ( "~ ", paste ( all_predictors, collapse = " + " ) ) ), 
                         data = DF_test_imp) [ , -1 ]
x_test_lab <- model.matrix ( as.formula ( paste ( "~ ", paste ( all_predictors, collapse = " + " ) ) ), 
                             data = DF_test_imp) [ , 2 : 39 ]
y_test_IS <- DF_test$wijzigingIS                                                # Does not matter if you use 'DF_test' or 'DF_test_imp' as the outcome variables are not imputed
y_test_Vaso <- DF_test$wijzigingVasoactief
y_test_Support <- DF_test$wijzingSupportmed

# Full Data
x_full <- model.matrix ( as.formula ( paste ( "~ ", paste ( all_predictors, collapse = " + " ) ) ), 
                         data = DF_full_imp) [ , -1 ]
x_full_lab <- model.matrix ( as.formula ( paste ( "~ ", paste ( all_predictors, collapse = " + " ) ) ), 
                             data = DF_full_imp) [ , 2 : 39 ]
y_full_IS <- DF_full_imp$wijzigingIS
y_full_Vaso <- DF_full_imp$wijzigingVasoactief
y_full_Support <- DF_full_imp$wijzigingVasoactief

# Define which variables should be penalised 
penalty_vars <- ifelse ( colnames ( x_train ) %in% lab_params, 1, 0 )

#------------------------------------------------------------------------------#
#----------------------------Immunosuppressants--------------------------------#
#------------------------------------------------------------------------------#

                      ### COVARIATES AND LAB PARAMETERS ###
# Define model 
model_IS <- glmnet ( x_train, y_train_IS, alpha = 1, family = "binomial", 
                     type.measure = "auc", penalty.factor = penalty_vars )

# Determine best lambda with cross-validation
train_cv_IS <- cv.glmnet ( x_train, y_train_IS, alpha = 1, family = "binomial", 
                           type.measure = "auc", penalty.factor = penalty_vars )
plot ( train_cv_IS )

best_lambda_IS <- train_cv_IS$lambda.1se                                        # Best lambda based on the 1SE rule

# Compute AUC for test set
pred_IS <- predict ( model_IS, s = best_lambda_IS, newx = x_test, type = "response" )
roc_IS <- roc ( y_test_IS, as.numeric ( pred_IS ) )
auc ( roc_IS )
plot ( roc_IS, print.auc = TRUE ) 

# Determine coefficients for full data
model_final_IS <- glmnet ( x_full, y_full_IS, alpha = 1,
                           family = "binomial", type.measure = "auc",
                           penalty.factor = penalty_vars )
coef_IS <- predict ( model_final_IS, type = "coefficients", s = best_lambda_IS )
coef_IS

# Compute AUC for full data
pred_full_IS <- predict ( model_final_IS, s = best_lambda_IS, newx = x_full, type = "response" )
roc_full_IS <- roc ( y_full_IS, as.numeric ( pred_full_IS ) )
auc ( roc_full_IS )

                            ### LAB PARAMETERS ONLY ###

# Define model 
model_lab_IS <- glmnet ( x_train_lab, y_train_IS, alpha = 1, family = "binomial", 
                         type.measure = "auc" )

# Determine best lambda with cross-validation
train_cv_lab_IS <- cv.glmnet ( x_train_lab, y_train_IS, alpha = 1, family = "binomial", 
                               type.measure = "auc" )
plot ( train_cv_lab_IS )

best_lambda_lab_IS <- train_cv_lab_IS$lambda.1se                                 # Best lambda based on the 1SE rule

# Compute AUC for test set
pred_lab_IS <- predict ( model_lab_IS, s = best_lambda_lab_IS, newx = x_test_lab, type = "response" )
roc_lab_IS <- roc ( y_test_IS, as.numeric ( pred_lab_IS ) )
auc ( roc_lab_IS )
plot ( roc_lab_IS, print.auc = TRUE ) 

# Determine coefficients for full data
model_final_lab_IS <- glmnet ( x_full_lab, y_full_IS, alpha = 1,
                               family = "binomial", type.measure = "auc", intercept = FALSE)
coef_lab_IS <- predict ( model_final_lab_IS, type = "coefficients", 
                         s = best_lambda_lab_IS )
coef_lab_IS

# Compute AUC for full data
pred_full_lab_IS <- predict ( model_final_lab_IS, s = best_lambda_lab_IS, 
                              newx = x_full_lab, type = "response" )
roc_full_lab_IS <- roc ( y_full_IS, as.numeric ( pred_full_lab_IS ) )
auc ( roc_full_lab_IS )

#----------------------------------Table IS------------------------------------#
                        ### COVARIATES AND LAB PARAMETERS ###
best_lambda_IS                                                                  # Optimal lambda based on 1SE rule
train_cv_IS$cvm [ which ( train_cv_IS$lambda == best_lambda_IS ) ]              # AUC training set
auc ( roc_IS )                                                                  # AUC test set
auc ( roc_full_IS )                                                             # AUC full data
sum ( coef_IS != 0 ) - 1                                                        # Number of predictors

                            ### LAB PARAMETERS ONLY ###
best_lambda_lab_IS                                                              # Optimal lambda based on 1SE rule
train_cv_lab_IS$cvm [ which ( train_cv_lab_IS$lambda == best_lambda_lab_IS ) ]  # AUC training set
auc ( roc_lab_IS )                                                              # AUC test set
auc ( roc_full_lab_IS )                                                         # AUC full data
sum ( coef_lab_IS != 0 ) - 1                                                    # Number of predictors



#------------------------------------------------------------------------------#
#-----------------------------Vasoactive Agents--------------------------------#
#------------------------------------------------------------------------------#

                      ### COVARIATES AND LAB PARAMETERS ###
# Define model 
model_Vaso <- glmnet ( x_train, y_train_Vaso, alpha = 1, family = "binomial", 
                       type.measure = "auc", penalty.factor = penalty_vars )

# Determine best lambda with cross-validation
train_cv_Vaso <- cv.glmnet ( x_train, y_train_Vaso, alpha = 1, family = "binomial", 
                             type.measure = "auc", penalty.factor = penalty_vars )
plot ( train_cv_Vaso )

best_lambda_Vaso <- train_cv_Vaso $ lambda.1se                                  # Best lambda based on the 1SE rule

# Compute AUC for test set
pred_Vaso <- predict ( model_Vaso, s = best_lambda_Vaso, newx = x_test, type = "response" )
roc_Vaso <- roc ( y_test_Vaso, as.numeric ( pred_Vaso ) )
auc ( roc_Vaso )

# Determine coefficients for full data
model_final_Vaso <- glmnet ( x_full, y_full_Vaso, alpha = 1,
                             family = "binomial", type.measure = "auc",
                             penalty.factor = penalty_vars )
coef_Vaso <- predict ( model_final_Vaso, type = "coefficients", s = best_lambda_Vaso )
coef_Vaso

# Compute AUC for full data
pred_full_Vaso <- predict ( model_final_Vaso, s = best_lambda_Vaso, newx = x_full, type = "response" )
roc_full_Vaso <- roc ( y_full_Vaso, as.numeric ( pred_full_Vaso ) )
auc ( roc_full_Vaso )

                          ### LAB PARAMETERS ONLY ###

# Define model 
model_lab_Vaso <- glmnet ( x_train_lab, y_train_Vaso, alpha = 1, family = "binomial", 
                           type.measure = "auc" )

# Determine best lambda with cross-validation
train_cv_lab_Vaso <- cv.glmnet ( x_train_lab, y_train_Vaso, alpha = 1, family = "binomial", 
                                 type.measure = "auc" )
plot ( train_cv_lab_Vaso )

best_lambda_lab_Vaso <- train_cv_lab_Vaso$lambda.1se                            # Best lambda based on the 1SE rule

# Compute AUC for test set
pred_lab_Vaso <- predict ( model_lab_Vaso, s = best_lambda_lab_Vaso, newx = x_test_lab, type = "response" )
roc_lab_Vaso <- roc ( y_test_Vaso, as.numeric ( pred_lab_Vaso ) )
auc ( roc_lab_Vaso )
plot ( roc_lab_Vaso, print.auc = TRUE ) 

# Determine coefficients for full data
model_final_lab_Vaso <- glmnet ( x_full_lab, y_full_Vaso, alpha = 1,
                                 family = "binomial", type.measure = "auc" )
coef_lab_Vaso <- predict ( model_final_lab_Vaso, type = "coefficients", 
                           s = best_lambda_lab_Vaso )
coef_lab_Vaso

# Compute AUC for full data
pred_full_lab_Vaso <- predict ( model_final_lab_Vaso, s = best_lambda_lab_Vaso, 
                                newx = x_full_lab, type = "response" )
roc_full_lab_Vaso <- roc ( y_full_IS, as.numeric ( pred_full_lab_Vaso ) )
auc ( roc_full_lab_Vaso )

#----------------------------------Table Vaso----------------------------------#

                        ### COVARIATES AND LAB PARAMETERS ###
best_lambda_Vaso                                                                # Optimal lambda based on 1SE rule
train_cv_Vaso$cvm [ which ( train_cv_Vaso$lambda == best_lambda_Vaso ) ]        # AUC training set
auc ( roc_Vaso )                                                                # AUC test set
auc ( roc_full_Vaso )                                                           # AUC full data
sum ( coef_Vaso != 0 ) - 1                                                      # Number of predictors

                            ### LAB PARAMETERS ONLY ###
best_lambda_lab_Vaso                                                            # Optimal lambda based on 1SE rule
train_cv_lab_Vaso$cvm [ which ( train_cv_lab_Vaso$lambda == best_lambda_lab_Vaso ) ]  # AUC training set
auc ( roc_lab_Vaso )                                                            # AUC test set
auc ( roc_full_lab_Vaso )                                                       # AUC full data
sum ( coef_lab_Vaso != 0 ) - 1                                                  # Number of predictors

#------------------------------------------------------------------------------#
#----------------------------Support Medication--------------------------------#
#------------------------------------------------------------------------------#

                      ### COVARIATES AND LAB PARAMETERS ###
# Define model 
model_Support <- glmnet ( x_train, y_train_Support, alpha = 1, 
                          family = "binomial", type.measure = "auc", 
                          penalty.factor = penalty_vars )

# Determine best lambda with cross-validation
train_cv_Support <- cv.glmnet ( x_train, y_train_Support, alpha = 1, 
                                family = "binomial", type.measure = "auc", 
                                penalty.factor = penalty_vars )
plot ( train_cv_Support )

best_lambda_Support <- train_cv_Support $ lambda.1se                            # Best lambda based on the 1SE rule

# Compute AUC for test set
pred_Support <- predict ( model_Support, s = best_lambda_Support, newx = x_test, 
                          type = "response" )
roc_Support <- roc ( y_test_Support, as.numeric ( pred_Support ) )
auc ( roc_Support )

# Determine coefficients for full data
model_final_Support <- glmnet ( x_full, y_full_Support, alpha = 1,
                                family = "binomial", type.measure = "auc",
                                penalty.factor = penalty_vars )
coef_Support <- predict ( model_final_Support, type = "coefficients", s = best_lambda_Support )
coef_Support

# Compute AUC for full data
pred_full_Support <- predict ( model_final_Support, s = best_lambda_Support, newx = x_full, type = "response" )
roc_full_Support <- roc ( y_full_Support, as.numeric ( pred_full_Support ) )
auc ( roc_full_Support )

                              ### LAB PARAMETERS ONLY ###
# Define model 
model_lab_Support <- glmnet ( x_train_lab, y_train_Support, alpha = 1, family = "binomial", 
                              type.measure = "auc" )

# Determine best lambda with cross-validation
train_cv_lab_Support <- cv.glmnet ( x_train_lab, y_train_Support, alpha = 1, family = "binomial", 
                                    type.measure = "auc" )
plot ( train_cv_lab_Support )

best_lambda_lab_Support <- train_cv_lab_Support$lambda.1se                      # Best lambda based on the 1SE rule

# Compute AUC for test set
pred_lab_Support <- predict ( model_lab_Support, s = best_lambda_lab_Support, newx = x_test_lab, type = "response" )
roc_lab_Support <- roc ( y_test_Support, as.numeric ( pred_lab_Support ) )
auc ( roc_lab_Support )
plot ( roc_lab_Support, print.auc = TRUE ) 

# Determine coefficients for full data
model_final_lab_Support <- glmnet ( x_full_lab, y_full_Support, alpha = 1,
                                    family = "binomial", type.measure = "auc" )
coef_lab_Support <- predict ( model_final_lab_Support, type = "coefficients", 
                              s = best_lambda_lab_Support )
coef_lab_Support

# Compute AUC for full data
pred_full_lab_Support <- predict ( model_final_lab_Support, s = best_lambda_lab_Support, 
                                   newx = x_full_lab, type = "response" )
roc_full_lab_Support <- roc ( y_full_Support, as.numeric ( pred_full_lab_Support ) )
auc ( roc_full_lab_Support )

#----------------------------------Table Support-------------------------------#

                        ### COVARIATES AND LAB PARAMETERS ###
best_lambda_Support                                                             # Optimal lambda based on 1SE rule
train_cv_Support$cvm [ which ( train_cv_Support$lambda == best_lambda_Support ) ]# AUC training set
auc ( roc_Support )                                                             # AUC test set
auc ( roc_full_Support )                                                        # AUC full data
sum ( coef_Support != 0 ) - 1                                                   # Number of predictors

                            ### LAB PARAMETERS ONLY ###
best_lambda_lab_Support                                                         # Optimal lambda based on 1SE rule
train_cv_lab_Support$cvm [ which ( train_cv_lab_Support$lambda == best_lambda_lab_Support ) ]  # AUC training set
auc ( roc_lab_Support )                                                         # AUC test set
auc ( roc_full_lab_Support )                                                    # AUC full data
sum ( coef_lab_Support != 0 ) - 1                                               # Number of predictors


