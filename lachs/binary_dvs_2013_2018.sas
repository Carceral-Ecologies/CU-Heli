/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Carlos Rodriguez, PhD. CU Anschutz Dept. of Family Medicine

RWJF - Noise helicopted and public health outcomes

This script fits adjusted and unadjusted models to the LACHS binary dependent 
variables:
  PHQ_num - 0,1 binary version of whether or not someone meets criteria for 
    psychological distress
  DEPNOW_num - 0,1, binary version of whether or not someone self reports 
    having depression
  ASTHMA_num - 0,1, binary version of whether or not someone has asthma

LACHS Methodology
http://publichealth.lacounty.gov/ha/docs/2022LACHS/Methodology/2023%20LACHS%20Methodology%20Report%20Final.pdf


Weight Variables:
POP_WGT - Represents the population sampling weights.

POP_SAMWGT - Represents the normalized population weight. Sums to the number of
individual responses in the sample.
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
proc format;
  picture pctfmt (round) other="009.9%";
run;


/* Map a library to a directory name ----------------------------------------*/
/* Set a variable to a the root data directory */
%let root = D:/rwjf_noise;

/* Loop through the years */
%let years = 
  2018
  2023
;

/* 
%let year = 2018; 
%let dvar = PHQ_num;
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
    PHQ_num
    DEPNOW_num
    ASTHMA_num
  ;


  %do j = 1 %to %sysfunc(countw(&dvars));
    %let dvar = %scan(&dvars, &j);

    ods word file = "&root./results/lachs/bin_&dvar._&year._glimmix_unadj.docx";
    title "Weighted Unadjusted GLMM - &dvar &year";

    ods output OneWayFreqs=freq_out;
    proc freq data=data;
      where not missing(GEOID)
        and not missing(pw_exp_db_cen)
        and not missing(&dvar);
      tables &dvar;
    run;
    ods output close;

    /* sub_table_1 */
    proc sql;
      create table sub_tab_1 as
      select
          "&year." as Year,
          "&dvar." as DV,
          CumFrequency as N,
          Frequency,
          Percent
      from freq_out
      where &dvar = 1;
    quit;

    /* Set the output for the unadjusted models */
    /* _pe is for the p value */
    ods output
        Tests3 = _pe
        OddsRatios = _or
    ;

    ods exclude ClassLevels;
    proc glimmix data=data method=quad empirical;

        class GEOID GEO_SUPD;

        model &dvar(event='1') = pw_exp_db_cen
            / dist=binomial
              link=logit
              obsweight=POP_SAMWGT
              solution
              oddsratio;

        random intercept / subject=GEOID;
        covtest 0;

    run;
    ods word close;


    /* Add the identifying dependent variable and month */
    /* Sub table 2 */
    proc sql;
      create table sub_tab_2 as
      select
          "&year." as Year,
          "&dvar." as DV,
          Estimate,
          Lower,
          Upper
      from _or;
    quit;

    /* Sub table 3 */
    proc sql;
      create table sub_tab_3 as
      select
        "&year." as Year,
        "&dvar." as DV,
        ProbF as P
      from _pe;
    quit;

    /* Combine the sub tables */
    proc sql;
        create table unadjusted_row_out as
        select a.*,
               b.Estimate,
               b.Lower,
               b.Upper,
               c.P
        from sub_tab_1 as a

        left join sub_tab_2 as b
            on a.DV = b.DV
            and a.Year = b.Year
            
        left join sub_tab_3 as c
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
proc export data=ALL_UNADJUSTED
    outfile="&root./results/lachs/sas/unadjusted_bin_dv_glimmix_.csv"
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

/* Test MACRO */
/* 
%let year=2023;
%let dvar=PHQ_num;
*/

%macro fit_adjusted_models;
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
    PHQ_num
    DEPNOW_num
    ASTHMA_num
  ;

  %do j = 1 %to %sysfunc(countw(&dvars));
    %let dvar = %scan(&dvars, &j);

    ods word file = "&root./results/lachs/bin_&dvar._&year._glimmix_adj.docx";

    title "Weighted Adjusted GLMM - &dvar &year";

    ods output OneWayFreqs=freq_out;
    proc freq data=data;
      where not missing(GEOID)
        and not missing(pw_exp_db_cen)
        and not missing(AGEGROUP_c)
        and not missing(GENDER)
        and not missing(RACESOP_c)
        and not missing(FPL_FIN)
        and not missing(GEO_SUPD);
      tables &dvar;
    run;

    /* sub_table_1 */
    proc sql;
      create table sub_tab_1 as
      select
          "&year." as Year,
          "&dvar." as DV,
          CumFrequency as N,
          Frequency,
          Percent
      from freq_out
      where &dvar = 1;
    quit;

    /* Set the output for the unadjusted models */
    ods output clear;

    ods output
        Tests3 = _pe_adj
        OddsRatios = _or_adj
    ;

    ods exclude ClassLevels;
    proc glimmix data=data method = quad empirical;

        class GEOID AGEGROUP_c GENDER RACESOP_c FPL_FIN GEO_SUPD;
        model &dvar(event='1') = pw_exp_db_cen AGEGROUP_c GENDER RACESOP_c FPL_FIN GEO_SUPD
            / dist=binomial
              link=logit
              obsweight=POP_SAMWGT
              solution
              oddsratio;

        random intercept / subject=GEOID; 
        
        covtest 0;

    run;
    ods word close;


    /* Add the identifying dependent variable and month and filter to the noise
      predictor*/
    /* Sub table 2 */
    proc sql outobs=1;
      create table sub_tab_2 as
      select
          "&year." as Year,
          "&dvar." as DV,
          Estimate,
          Lower,
          Upper
      from _or_adj;
    quit;


    /* Sub table 3 */
    proc sql;
      create table sub_tab_3 as
      select
        "&year." as Year,
        "&dvar." as DV,
        ProbF as P
      from _pe_adj
      where effect = "pw_exp_db_cen";
    quit;

    /* Combine the sub tables */
    proc sql;
        create table adjusted_row_out as
        select a.*,
               b.Estimate,
               b.Lower,
               b.Upper,
               c.P
        from sub_tab_1 as a

        left join sub_tab_2 as b
            on a.DV = b.DV
            and a.Year = b.Year
            
        left join sub_tab_3 as c
            on a.DV = b.DV
            and a.Year = b.Year
        ;
    quit;

    /* Append to a master table  */
    proc append base=all__adjusted data= adjusted_row_out force;
    run;

  %end;

%end;

%mend;

%fit_adjusted_models;

/* Output to .csv files */
proc export data=ALL_ADJUSTED
    outfile="&root./results/lachs/sas/adjusted_bin_dv_glimmix.csv"
    dbms=csv
    replace;
run;