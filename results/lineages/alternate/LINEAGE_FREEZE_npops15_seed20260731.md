# Pre-effect transmission-lineage freeze

- Alignment: frozen 11,550-position core-SNP alignment (not redistributed)
- Genomes: 989
- Parsimony-informative core SNPs: 11550
- hierBAPS maximum populations: 15
- Random seed: 20260731
- The clustering was defined without using case counts or post-clustering growth estimates.
- Level-1 populations are candidate model lineages; level-2 populations are sensitivity sublineages.
- A primary lineage requires at least 20 genomes and cross-country or cross-period coverage.

## Level-1 populations

   model_lineage_id n_genomes n_focal n_countries
             <char>     <int>   <int>       <int>
1:            L1_02       750     574          26
2:            L1_01       190     174          13
3:            L1_03        46      26           9
4:            L1_04         3       0           3
                                                                                                 countries
                                                                                                    <char>
1: ARG;AUS;AUT;BEL;BRA;CAN;CHN;DEU;DNK;ESP;FIN;FRA;GTM;HTI;IND;IRL;IRN;ISR;ITA;JPN;MEX;NLD;NOR;SWE;TUN;USA
2:                                                     AUS;BEL;CHN;DNK;ESP;FRA;IND;IRN;JPN;NOR;POL;SWE;USA
3:                                                                     AUT;CHN;DNK;IND;IRN;JPN;KEN;NLD;USA
4:                                                                                             DEU;IRN;USA
   n_periods                         periods min_year max_year
       <int>                          <char>    <int>    <int>
1:         3 pandemic;prepandemic;resurgence     1994     2025
2:         3 pandemic;prepandemic;resurgence     1962     2024
3:         2          prepandemic;resurgence     1939     2024
4:         1                     prepandemic     1952     2015
   primary_model_eligible
                   <lgcl>
1:                   TRUE
2:                   TRUE
3:                   TRUE
4:                  FALSE
                                                            exclusion_reason
                                                                      <char>
1:                                                                          
2:                                                                          
3:                                                                          
4: Fewer than 20 genomes or insufficient cross-country/cross-period coverage
