/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine

RWJF - Noise helicopted and public health outcomes

This script fits adjusted and unadjusted models to the LACHS count dependent 
variables:
  HSMH - Number of poor mental health days in the past 30-days
  HSPH - Number of poor physical health days in the past 30-days
  HSALD - Number of activity limited days in the past 30-days, due to poor 
  physical or mental health


LACHS Methodology
http://publichealth.lacounty.gov/ha/docs/2022LACHS/Methodology/2023%20LACHS%20Methodology%20Report%20Final.pdf


Weight Variables:
POP_WGT - Represents the population sampling weights.

POP_SAMWGT - Represents the normalized population weight. Sums to the number of
individual responses in the sample.
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Map a library to a directory name ----------------------------------------*/
/* Set a variable to a the root data directory */
%let root = D:/rwjf_noise;

/* Loop through the years */
%let years = 
  2018
  2023
;

/* Test MACRO */
/* 
%let year=2023;
%let dvar=HSPH;
%let dvar=HSMH;
%let dvar=HSALD;
*/


/* /////////////////////////// Unadjusted models /////////////////////////// */
%macro fit_unadjusted_models;
%do i = 1 %to %sysfunc(countw(&years));
  %let year = %scan(&years, &i);

  /* Import the data file with the corresponding year */
  proc import datafile = "&root./data/prepped_data/merged_lachs_&year..csv"
    out = data
    dbms = csv
    replace;
  run;

  /* Loop through the variables */
  %let dvars = 
    HSPH
    HSMH
    HSALD
  ;


  %do j = 1 %to %sysfunc(countw(&dvars));
    %let dvar = %scan(&dvars, &j);

    ods word file = "&root./results/lachs/count_&dvar._&year._glimmix_unadj.docx";

    title "Weighted Unadjusted GLMM - &dvar &year";
    
    proc means data=data mean var std;
      var &dvar;
    run;

    proc sgplot data=data;
      where not missing(&dvar);
      histogram &dvar;
    run;

    ods output OneWayFreqs=freq_out;
    proc freq data=data;
      where not missing(GEOID)
        and not missing(pw_exp_db_cen)
        and not missing(&dvar);
      tables &dvar;
    run;

    /* Set the output for the unadjusted models */
    ods output
        Estimates= _pe
        NObs = nobs
    ;

    ods exclude ClassLevels;
    proc glimmix data=data method=quad empirical;

        class GEOID;

        model &dvar = pw_exp_db_cen
            / dist=negbin
              link=log
              obsweight=POP_SAMWGT
              solution
              oddsratio;

        random intercept / subject=GEOID;

        estimate "pw_exp_db_cen" pw_exp_db_cen 1 / exp cl;

        covtest 0;

    run;
    ods word close;

    /* sub_table_1 */
    proc sql;
      create table sub_tab_1 as
      select
          "&year." as Year,
          "&dvar." as DV,
          N
      from nobs
      where Label = "Number of Observations Used";
    quit;

    /* sub_table_2 */
    proc sql outobs=1;
      create table sub_tab_2 as
      select
          "&year." as Year,
          "&dvar." as DV,
          ExpEstimate as IRR,
          ExpLower as Lower,
          ExpUpper as Upper,
          Probt as P
      from _PE;
    quit;

    /* Combine the sub tables */
    proc sql;
        create table unadjusted_row_out as
        select a.*,
               b.IRR,
               b.Lower,
               b.Upper,
               b.P
        from sub_tab_1 as a

        left join sub_tab_2 as b
            on a.DV = b.DV
            and a.Year = b.Year
        ;
    quit;

  /* Append to a master table  */
    proc append base=all_unadjusted data= unadjusted_row_out force;
    run;

  %end;

%end;

%mend;

%fit_unadjusted_models;

/* Output to .csv files */
proc export data=ALL_PE
    outfile="&root./results/lachs/sas/unadjusted_count_dv_glimmix_PEs.csv"
    dbms=csv
    replace;
run;



/* //////////////////////////// ADJUSTED MODELS //////////////////////////// */
/* Loop through the years */
%let root = D:/rwjf_noise;

%let years = 
  2018
  2023
;


/* 
%let year = 2018; 
%let dvar = HSPH;
*/
%macro fit_models;
%do i = 1 %to %sysfunc(countw(&years));
  %let year = %scan(&years, &i);

  /* Import the data file with the corresponding year */
  proc import datafile = "&root./data/prepped_data/merged_lachs_&year..csv"
    out = data
    dbms = csv
    replace;
  run;

  /* Loop through the variables */
  %let dvars = 
    HSPH
    HSMH
    HSALD
  ;

  %do j = 1 %to %sysfunc(countw(&dvars));
    %let dvar = %scan(&dvars, &j);
    
    ods word file = "&root./results/lachs/count_&dvar._&year._glimmix_adj.docx";

    title "Weighted Adjusted GLMM - &dvar &year";

    /* Set the output for the unadjusted models */
    ods output clear;

    ods output
        Estimates= _pe_adj
        NObs = nobs
    ;

    ods exclude ClassLevels;
    proc glimmix data=data method = quad empirical;

        class GEOID AGEGROUP_c GENDER RACESOP_c FPL_FIN GEO_SUPD;
        model &dvar = pw_exp_db_cen AGEGROUP_c GENDER RACESOP_c FPL_FIN GEO_SUPD
            / dist=negbin
              link=log
              obsweight=POP_SAMWGT
              solution
              oddsratio;

        random intercept / subject=GEOID; 

        /* To output the IRR */
        estimate "pw_exp_db_cen" pw_exp_db_cen 1 / exp cl;

        covtest 0;

    run;
    ods word close;

    /* sub_table_1 */
    proc sql;
      create table sub_tab_1 as
      select
          "&year." as Year,
          "&dvar." as DV,
          N
      from nobs
      where Label = "Number of Observations Used";
    quit;

    /* sub_table_2 */
    proc sql outobs=1;
      create table sub_tab_2 as
      select
          "&year." as Year,
          "&dvar." as DV,
          ExpEstimate as IRR,
          ExpLower as Lower,
          ExpUpper as Upper,
          Probt as P
      from _PE_ADJ;
    quit;

    /* Combine the sub tables */
    proc sql;
        create table adjusted_row_out as
        select a.*,
               b.IRR,
               b.Lower,
               b.Upper,
               b.P
        from sub_tab_1 as a

        left join sub_tab_2 as b
            on a.DV = b.DV
            and a.Year = b.Year
        ;
    quit;

    /* Append to a master table  */
    proc append base=all_adjusted data= adjusted_row_out force;
    run;

  %end;

%end;

%mend;

%fit_adjusted_models;

/* Output to .csv files */
proc export data=ALL_PE_ADJ
    outfile="&root./results/lachs/sas/adjusted_count_dv_glimmix.csv"
    dbms=csv
    replace;
run;