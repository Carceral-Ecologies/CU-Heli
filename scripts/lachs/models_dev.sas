/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine

RWJF - Noise helicopted and public health outcomes

LACHS Methodology
http://publichealth.lacounty.gov/ha/docs/2022LACHS/Methodology/2023%20LACHS%20Methodology%20Report%20Final.pdf

Current Data are 2018

Weight Variables:
POP_WGT - Represents the sampling weights, their sum is close to 8million
which is reasonable value for adult population of LA county in 2018.

POP_SAMWGT - Represents the normalized population weight, it has a mean of 1, 
values range from 0.1 to 9.53
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Map a library to a directory name ----------------------------------------*/
/* Set a variable to a the root data directory */
%let root = D:/rwjf_noise;

/* Import data ------------------------------------------------------------- */

proc import datafile = "&root./data/prepped_data/merged_lachs_2018.csv"
    out = data
    dbms = csv
    replace;
run;

/* Set a library name to view raw data */
libname lacph "D:\rwjf_noise\data\LACPH";

/* ///////////////////////////////////////////////////////////////////////// */
/* ////////////////////////PHQ Unadjusted Models /////////////////////////// */
/* ///////////////////////////////////////////////////////////////////////// */

/* --------------------------- BEGIN WORD OUTPUT --------------------------- */
ods word file = "&root./results/lachs/models.docx";

proc freq data=data;
    tables PHQ_num;
run;


/* //////////////////////////////////////////////////////////////////////// */
/* Unadjusted GLM Model */
/* //////////////////////////////////////////////////////////////////////// */
proc glimmix data=data method=laplace;
    model PHQ_num(event='1') = pw_exp_db_cen
        / dist=binomial link=logit solution;
run;

/* //////////////////////////////////////////////////////////////////////// */
/* Unadjusted GLMM Model */
/* //////////////////////////////////////////////////////////////////////// */
proc glimmix data=data method=laplace;
    class GEOID;
    
    model PHQ_num(event='1') = pw_exp_db_cen
        / dist=binomial link=logit solution;
    
    random intercept / subject=GEOID;
run;

/* //////////////////////////////////////////////////////////////////////// */
/* Weighted Unadjusted GLMM Model WEIGHT argument (ignoring clustering)*/
/* //////////////////////////////////////////////////////////////////////// */

/* The coefficients change when adding the weight and the SEs grow
substantially. The change in the coefficient suggests a mis-specification */
proc glimmix data=data method=laplace;
    class GEOID;

    model PHQ_num(event='1') = pw_exp_db_cen
        / dist=binomial link=logit solution;

    random intercept / subject=GEOID;

    weight POP_WGT;
run;

/* //////////////////////////////////////////////////////////////////////// */
/* Weighted Unadjusted GLMM Model OBS WEIGHT argument (Carter)  */
/* Uses the normalized weights, otherwise the model fails  */
/* obsweights are only supported for Method = QUAD */
/* //////////////////////////////////////////////////////////////////////// */

ods exclude ClassLevels;
proc glimmix data=data method=quad empirical;

    class GEOID;

    model PHQ_num(event='1') = pw_exp_db_cen
        / dist=binomial
          link=logit
          obsweight=POP_SAMWGT
          solution
          oddsratio;

    random intercept / subject=GEOID;

    covtest 0;

run;



/* //////////////////////////////////////////////////////////////////////// */
/* Unweighted Adjusted GLM */
/* //////////////////////////////////////////////////////////////////////// */

title "UnWeighted Adjusted GLM";
proc glimmix data=data;
    class GEOID AGEGROUP_c GENDERO RACESOP_c FPL_FIN;

    model PHQ_num(event='1') = pw_exp_db_cen AGEGROUP_c GENDERO RACESOP_c FPL_FIN 
        / dist=binomial link=logit solution;
run;

/* //////////////////////////////////////////////////////////////////////// */
/* Weighted Adjusted GLM */
/* //////////////////////////////////////////////////////////////////////// */

title "Population Weighted Adjusted GLM";
proc glimmix data=data;
    class GEOID AGEGROUP_c GENDERO RACESOP_c FPL_FIN;

    model PHQ_num(event='1') = pw_exp_db_cen AGEGROUP_c GENDERO RACESOP_c FPL_FIN 
        / dist=binomial link=logit solution;

    weight POP_WGT;
run;

/* //////////////////////////////////////////////////////////////////////// */
/* Normalized Weighted Adjusted GLM */
/* //////////////////////////////////////////////////////////////////////// */

/* Shows that the PE is similar as weighted adjusted but the SE is larger, more like the unweighted adjusted model */
title "Normalized Weighted Adjusted GLM";
proc glimmix data=data;
    class GEOID AGEGROUP_c GENDERO RACESOP_c FPL_FIN;

    model PHQ_num(event='1') = pw_exp_db_cen AGEGROUP_c GENDERO RACESOP_c FPL_FIN 
        / dist=binomial link=logit solution;

    weight POP_SAMWGT;
run;


/* //////////////////////////////////////////////////////////////////////// */
/* Unweighted Adjusted GLMM */
/* //////////////////////////////////////////////////////////////////////// */
title "Unweighted Adjusted GLMM";
proc glimmix data=data;
    class GEOID AGEGROUP_c GENDERO RACESOP_c FPL_FIN;

    model PHQ_num(event='1') = pw_exp_db_cen AGEGROUP_c GENDERO RACESOP_c FPL_FIN 
        / dist=binomial link=logit solution;

    random intercept / subject=GEOID;
run;

/* //////////////////////////////////////////////////////////////////////// */
/* Weighted Adjusted GLMM */
/* //////////////////////////////////////////////////////////////////////// */
title "Weighted Adjusted GLMM";
proc glimmix data=data;
    class GEOID AGEGROUP_c GENDERO RACESOP_c FPL_FIN;

    model PHQ_num(event='1') = pw_exp_db_cen AGEGROUP_c GENDERO RACESOP_c FPL_FIN 
        / dist=binomial link=logit solution;

    random intercept / subject=GEOID;

    /* This is analagous to lme4 weight argument (prior/precision weights) */
    weight POP_WGT;
run;


/* //////////////////////////////////////////////////////////////////////// */
/* Weighted adjusted GLMM Model OBS WEIGHT argument  */
/* Uses the normalized weights, otherwise the model fails  */
/* obsweights are only supported for Method = QUAD */
/* //////////////////////////////////////////////////////////////////////// */

ods exclude ClassLevels;
proc glimmix data=data method=quad;

    class GEOID AGEGROUP_c GENDERO RACESOP_c FPL_FIN;

    model PHQ_num(event='1') = pw_exp_db_cen AGEGROUP_c GENDERO RACESOP_c FPL_FIN
        / dist=binomial
          link=logit
          obsweight=POP_SAMWGT
          solution
          oddsratio;

    random intercept / subject=GEOID;

    covtest 0;

run;



ods word close;
/* ---------------------------- END WORD OUTPUT ---------------------------- */

