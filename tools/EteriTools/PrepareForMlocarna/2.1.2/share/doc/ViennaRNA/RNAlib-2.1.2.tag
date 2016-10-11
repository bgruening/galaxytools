<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>main</name>
    <title>ViennaRNA Package core - RNAlib</title>
    <filename>main</filename>
    <docanchor file="main">mp_intro</docanchor>
  </compound>
  <compound kind="file">
    <name>2Dfold.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>2Dfold_8h</filename>
    <member kind="function">
      <type>TwoDfold_vars *</type>
      <name>get_TwoDfold_variables</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>g490d7a00dc1a6188bb401438fd975dbe</anchor>
      <arglist>(const char *seq, const char *structure1, const char *structure2, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroy_TwoDfold_variables</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>g9837b781df6e7274d848c137fb438bb5</anchor>
      <arglist>(TwoDfold_vars *our_variables)</arglist>
    </member>
    <member kind="function">
      <type>TwoDfold_solution *</type>
      <name>TwoDfoldList</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>g131c9c01e5203a1d0f8ab14ab9eb8315</anchor>
      <arglist>(TwoDfold_vars *vars, int distance1, int distance2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDfold_backtrack_f5</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>gf453ac72055554d5057ddbb166a67999</anchor>
      <arglist>(unsigned int j, int k, int l, TwoDfold_vars *vars)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>2Dpfold.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>2Dpfold_8h</filename>
    <member kind="function">
      <type>TwoDpfold_vars *</type>
      <name>get_TwoDpfold_variables</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>g674dad118a59dc002516fe02b6613912</anchor>
      <arglist>(const char *seq, const char *structure1, char *structure2, int circ)</arglist>
    </member>
    <member kind="function">
      <type>TwoDpfold_vars *</type>
      <name>get_TwoDpfold_variables_from_MFE</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>g1a04119ebae7e52b9f6e05be864cdb91</anchor>
      <arglist>(TwoDfold_vars *mfe_vars)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroy_TwoDpfold_variables</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>ge6a9242efcf43798acad404334bccad4</anchor>
      <arglist>(TwoDpfold_vars *vars)</arglist>
    </member>
    <member kind="function">
      <type>TwoDpfold_solution *</type>
      <name>TwoDpfoldList</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>g4b9ecb276d9ae3c1da9e5ebc04c61be0</anchor>
      <arglist>(TwoDpfold_vars *vars, int maxDistance1, int maxDistance2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDpfold_pbacktrack</name>
      <anchorfile>group__kl__neighborhood__stochbt.html</anchorfile>
      <anchor>g0b1553be92cca053c3bb5aa1c2b2fd23</anchor>
      <arglist>(TwoDpfold_vars *vars, int d1, int d2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDpfold_pbacktrack5</name>
      <anchorfile>group__kl__neighborhood__stochbt.html</anchorfile>
      <anchor>g319e7822ddbe0a5a3c9a6568c2e2a33a</anchor>
      <arglist>(TwoDpfold_vars *vars, int d1, int d2, unsigned int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>alifold.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>alifold_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>update_alifold_params</name>
      <anchorfile>alifold_8h.html</anchorfile>
      <anchor>92a7f6d293182abbaaa349849127e36b</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alifold</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>g86c43a79871b43bbee379b7b8f39d564</anchor>
      <arglist>(const char **strings, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>circalifold</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>gc245f44507e3324a7cf04879a5a7f8c3</anchor>
      <arglist>(const char **strings, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_alifold_arrays</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>g2f8550298ee654e25e44c30bb242bbd1</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_mpi</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gcffae3b2d3a3f38b830e89488f2242c0</anchor>
      <arglist>(char *Alseq[], int n_seq, int length, int *mini)</arglist>
    </member>
    <member kind="function">
      <type>float **</type>
      <name>readribosum</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g50ef3f5b786c98b20169c7bb6dded8d3</anchor>
      <arglist>(char *name)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_alistruct</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g13133e191923ba173dc40c9358ed861f</anchor>
      <arglist>(const char **sequences, const char *structure, int n_seq, float *energy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>encode_ali_sequence</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g12ba60c4e7278aa0cead337071aebc49</anchor>
      <arglist>(const char *sequence, short *S, short *s5, short *s3, char *ss, unsigned short *as, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>alloc_sequence_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g393d56ec257e2ec4bcb26b27dd9df7c1</anchor>
      <arglist>(const char **sequences, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_sequence_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g160c03a96eb91c6dd35aabed49f215e8</anchor>
      <arglist>(unsigned int n_seq, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_fold_par</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>g4159d3d8bb0fbb5b0be41fd54ad69030</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl, pf_paramT *parameters, int calculate_bppm, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_fold</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gca4fc5a8337db6695f0fb9fd4f5f0530</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_circ_fold</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gd6ee46b5bd800138407994fd1efc2197</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_ali_bppm</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gc59b5cbb1810e9331bc91edfd388365a</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>alipbacktrack</name>
      <anchorfile>group__consensus__stochbt.html</anchorfile>
      <anchor>g59ce07041a74911ca72ea5efa0bb3478</anchor>
      <arglist>(double *prob)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_alipf_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gd6f7cd96aeaa838e32ada56b350a22f3</anchor>
      <arglist>(short ***S_p, short ***S5_p, short ***S3_p, unsigned short ***a2s_p, char ***Ss_p, double **qb_p, double **qm_p, double **q1k_p, double **qln_p, short **pscore)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>cv_fact</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gd6fac190cfd3af07934dfbba7dfdf636</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>nc_fact</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g74ff4bb41b72ed967db7dce00fa66edc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>cofold.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>cofold_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>cofold</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>g6de084332decd07c52ede3d20af76d95</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>cofold_par</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>g83d04c7fb08ac05039575183545e31df</anchor>
      <arglist>(const char *string, char *structure, paramT *parameters, int is_constrained)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_co_arrays</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>gf31159249b8aaadfb95aa3997d422541</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_cofold_params</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>g95d130f8c6f091db0d90816857001e21</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>export_cofold_arrays_gq</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>ga3206b80db2f82bbff8296f0673af991</anchor>
      <arglist>(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **fc_p, int **ggg_p, int **indx_p, char **ptype_p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>export_cofold_arrays</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>gcdfda1fefb12716b1bfb2a7cacd0e5dc</anchor>
      <arglist>(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **fc_p, int **indx_p, char **ptype_p)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>zukersubopt</name>
      <anchorfile>group__subopt__zuker.html</anchorfile>
      <anchor>gec07305955df39b1d579611c0039ddb9</anchor>
      <arglist>(const char *string)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>zukersubopt_par</name>
      <anchorfile>group__subopt__zuker.html</anchorfile>
      <anchor>g38aa4496326af32ea9ec137a9bf71ffe</anchor>
      <arglist>(const char *string, paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_monomere_mfes</name>
      <anchorfile>cofold_8h.html</anchorfile>
      <anchor>fdb567aa963260db6415b2eb2f6bba31</anchor>
      <arglist>(float *e1, float *e2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize_cofold</name>
      <anchorfile>cofold_8h.html</anchorfile>
      <anchor>2b4f882e015d5bf57d54c9f89f1e7170</anchor>
      <arglist>(int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>convert_epars.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>convert__epars_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_ALL</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g80ca4b03c4e33f8224bd74cb65ea8f21</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g48f634f8a30d7a76e9f671216a8e0a41</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_STACK</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gcc7cee41b96c907b1661acce4a1bb20c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g758f179bbca896409b22df663fa72ae7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga7bc7f81c5188c4e9d7d72fefec8417a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT_1N</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gb3defd4125dac42a477b74ef967b416d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT_23</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g4812fbe4be82f8d3a34b08d5a809d868</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_MULTI</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g8c9ffd43ac2cf9dee4b5923c2367d5c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_EXT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g2bc3edd21d014b0891a567eaa6aff758</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DANGLE5</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gd3988d62d6ebe9e5d220aea78512c3e8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DANGLE3</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g6c3eab5261e2cd3fe6a30e51406bf9d2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_11</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g54502f72c1a16d08acd233d06242cc0f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_21</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga0b6db305805e42ca08488f8436e9536</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_22</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gfe4f9deeb34ccc2400203aed4f494685</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_BULGE</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g22004edd37d0efff575f28d0f79f363d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g8415ba389c7c6d0fe9621c395cfb9fe0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_ML</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g5356c521b724fad6b5272d22cb14adb9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MISC</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g08304cb51f65c09dae2fc4eeb5e6277e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_SPECIAL_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ge42ee1ec21bd62f632bb554faf336667</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_VANILLA</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g6a704a05e53ee38f37564ed36122134e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_NINIO</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga299c4d0a5ea47e7301a304d2d6aa9ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DUMP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g1aabc3c11923e5253d6f001cf19a8f95</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convert_parameter_file</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g3855fb2de2943aed520512eab2f1c7c3</anchor>
      <arglist>(const char *iname, const char *oname, unsigned int options)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>data_structures.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>data__structures_8h</filename>
    <class kind="struct">plist</class>
    <class kind="struct">cpair</class>
    <class kind="struct">COORDINATE</class>
    <class kind="struct">sect</class>
    <class kind="struct">bondT</class>
    <class kind="struct">bondTEn</class>
    <class kind="struct">model_detailsT</class>
    <class kind="struct">paramT</class>
    <class kind="struct">pf_paramT</class>
    <class kind="struct">PAIR</class>
    <class kind="struct">INTERVAL</class>
    <class kind="struct">SOLUTION</class>
    <class kind="struct">cofoldF</class>
    <class kind="struct">ConcEnt</class>
    <class kind="struct">pairpro</class>
    <class kind="struct">pair_info</class>
    <class kind="struct">move</class>
    <class kind="struct">intermediate</class>
    <class kind="struct">path</class>
    <class kind="struct">pu_contrib</class>
    <class kind="struct">interact</class>
    <class kind="struct">pu_out</class>
    <class kind="struct">constrain</class>
    <class kind="struct">duplexT</class>
    <class kind="struct">node</class>
    <class kind="struct">snoopT</class>
    <class kind="struct">dupVar</class>
    <class kind="struct">TwoDfold_solution</class>
    <class kind="struct">TwoDfold_vars</class>
    <class kind="struct">TwoDpfold_solution</class>
    <class kind="struct">TwoDpfold_vars</class>
    <member kind="define">
      <type>#define</type>
      <name>MAXALPHA</name>
      <anchorfile>data__structures_8h.html</anchorfile>
      <anchor>89a2a7432783b689b8913a3486cf04ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAXDOS</name>
      <anchorfile>data__structures_8h.html</anchorfile>
      <anchor>2280be174b42d3ae9f8350a0f2302bbe</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>dist_vars.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>dist__vars_8h</filename>
    <class kind="struct">Postorder_list</class>
    <class kind="struct">Tree</class>
    <class kind="struct">swString</class>
    <member kind="variable">
      <type>int</type>
      <name>edit_backtrack</name>
      <anchorfile>dist__vars_8h.html</anchorfile>
      <anchor>17ff12d422d6e292bf715cbcce0fdcbc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>aligned_line</name>
      <anchorfile>dist__vars_8h.html</anchorfile>
      <anchor>1a89dc61fb31085208cf7af48e6a32a2</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>cost_matrix</name>
      <anchorfile>dist__vars_8h.html</anchorfile>
      <anchor>bfa2482b8ad5a3bf84df152a5cd395a2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>duplex.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>duplex_8h</filename>
  </compound>
  <compound kind="file">
    <name>edit_cost.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>edit__cost_8h</filename>
  </compound>
  <compound kind="file">
    <name>energy_const.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>energy__const_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>GASCONST</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>d17105da6a9a514fdb8c560adea0f282</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>K0</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>161bf96c4308ee2d9a532be866ecbb84</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>INF</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>16860b9e6d3b0aa03435fe499b1597ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>FORBIDDEN</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>fcc7369ed3414bd0bb74494b01e421c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>BONUS</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>720489e16392267ce01036f1182f666a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>NBPAIRS</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>11aeaf479a69853f8dacb0ed779d416e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>TURN</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>ffd7252901b04ee30e14034f9fa58ec7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAXLOOP</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>3fa4041bfda5ee538c126a1ec2c471e6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>findpath.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>findpath_8h</filename>
    <member kind="function">
      <type>int</type>
      <name>find_saddle</name>
      <anchorfile>findpath_8h.html</anchorfile>
      <anchor>78576db358b3b9e26accc1257378625b</anchor>
      <arglist>(const char *seq, const char *struc1, const char *struc2, int max)</arglist>
    </member>
    <member kind="function">
      <type>path_t *</type>
      <name>get_path</name>
      <anchorfile>findpath_8h.html</anchorfile>
      <anchor>be361edc74af6080856eec0f52c7b420</anchor>
      <arglist>(const char *seq, const char *s1, const char *s2, int maxkeep)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_path</name>
      <anchorfile>findpath_8h.html</anchorfile>
      <anchor>2de5fbb89658040d45f10db89d896531</anchor>
      <arglist>(path_t *path)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>fold.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>fold_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>fold_par</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gf3d4bf3ba3c37140d643203b67112495</anchor>
      <arglist>(const char *sequence, char *structure, paramT *parameters, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>fold</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gd11bb36c69e7e75f0bbae81a6642a015</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>circfold</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>g85be769330de2b504a1c130fe55af166</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_structure</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g3532b3f697f33a33cf02c774b6f49c15</anchor>
      <arglist>(const char *string, const char *structure, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_struct_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>gfe8e6f53f3641396e718bb45269e09e6</anchor>
      <arglist>(const char *string, const char *structure, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_structure</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g164a4b025da054f7f4445eba3652fd4a</anchor>
      <arglist>(const char *string, const char *structure, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_struct_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g46602e8377fad2621b873654bfd976f1</anchor>
      <arglist>(const char *string, const char *structure, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_structure_pt</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>ge478bf8b55469a82b580ceb9d198849d</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_struct_pt_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g2d5a0f1e2797a05d09d1867a8c61554a</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_arrays</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gb823f35d01c03e9c10823c0fd002c8f6</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parenthesis_structure</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>f8c83f08bcda93b6e6c7de80e5c5e57d</anchor>
      <arglist>(char *structure, bondT *bp, int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parenthesis_zuker</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>ef15155a540163009457eb39147b8382</anchor>
      <arglist>(char *structure, bondT *bp, int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_fold_params</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>g44d079d97f158f294c61acf250f22582</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_move</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>c9330596465ca14f0522dfe0fd56dbf6</anchor>
      <arglist>(const char *string, const char *structure, int m1, int m2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_move_pt</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>aa3fa66452fc0aa837832b22e10e4767</anchor>
      <arglist>(short *pt, short *s, short *s1, int m1, int m2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>loop_energy</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>49152b530b2c0f63d93d24b0beb8a5c3</anchor>
      <arglist>(short *ptable, short *s, short *s1, int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>assign_plist_from_db</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>9527a44c5061809435182d58912437d0</anchor>
      <arglist>(plist **pl, const char *struc, float pr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>LoopEnergy</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>13c10a2632985fb581dc6b059ebb1baa</anchor>
      <arglist>(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>HairpinE</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>474304e5456d413c2741b681e30efd0c</anchor>
      <arglist>(int size, int type, int si1, int sj1, const char *string)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize_fold</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>2f4e4f07ce2159c55e83998eca1d5199</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_struct</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a752ff04d38c41bd3767849e7bb7e022</anchor>
      <arglist>(const char *string, const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_struct_pt</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>ffa0ea3c9377335cf52bd441ae254dd9</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_struct</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>d7d32e68fa5e1576abe4b0936d67f94f</anchor>
      <arglist>(const char *string, const char *structure)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>logML</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>f66b4b895e811bf11e6d9d9a261a1632</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>uniq_ML</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a7aa8c598788d963d1c7e16c6eb859ed</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>cut_point</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>d8b86ee9a485d93a0d0a52e683e8ca61</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>eos_debug</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g05f1af959998612960ddcf1436fa1837</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>fold_vars.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>fold__vars_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>set_model_details</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>1235f6733808c132ddba9f14d74f7fbf</anchor>
      <arglist>(model_detailsT *md)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>fold_constrained</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>85ce01be2d212e2043761bdb64fcce60</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noLonelyPairs</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>f65a8bb2c18daac48c2ec0f12ff34759</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>dangles</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>b389914dc21770ef9da67aac4936f234</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noGU</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>8b2312ab1f28ae79a65ecf8a6444e702</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>no_closingGU</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>c1cfeadb4b4d880f545bb87c6afd4378</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>tetra_loop</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>996b2aaf20a173b09b3ee371ab9a0bd8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>energy_set</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>c6cf88260813dd4f2f67672183ffe5b4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>circ</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>48f19216ecf3d2d651abc107171c80bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>csv</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>275ef10a691bcc254206c0cf4c8ebba9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>oldAliEn</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>77ecfb7cd618d7bf69bbd0e910dde871</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ribo</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>169a0471a50fb17454bd5fa262c7e519</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>RibosumFile</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>53536913ff8d5f2f1581d6e8576ccea9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>nonstandards</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>9e83de561900f59058bd523eb071cabd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>temperature</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>0b516316c5f7ab02f8d1a9fbc1ba2fbb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>james_rule</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>28c3777fb816f026fb32e51340de5ea3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>logML</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>f66b4b895e811bf11e6d9d9a261a1632</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>cut_point</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>d8b86ee9a485d93a0d0a52e683e8ca61</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bondT *</type>
      <name>base_pair</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>aec1cda581406cac521e8a55c8c1a4cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>pr</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>75a6fe753fbc832775d21e537ad1b21e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>iindx</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>e2e64f77fe25400c55cd38b8270db573</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>pf_scale</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>1f54ec895a6e5961f82e188bafbea432</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>do_backtrack</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>d18a5668c071dfcf0cb088899aa32106</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>backtrack_type</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>eb0d56969f1dcd94da7fca7d1bf9c1f6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>gquad</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>6f60ecae21447dde4e31a4872690c817</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>gquad.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>gquad_8h</filename>
    <member kind="function">
      <type>int *</type>
      <name>get_gquad_matrix</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>f7270ad770e537ea838f65776f666df7</anchor>
      <arglist>(short *S, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>parse_gquad</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>46cb26be977b21caf37f6d00504ff303</anchor>
      <arglist>(const char *struc, int *L, int l[3])</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>backtrack_GQuad_IntLoop</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>28444e8a5d771a46fd2167f352de8a3a</anchor>
      <arglist>(int c, int i, int j, int type, short *S, int *ggg, int *index, int *p, int *q, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>backtrack_GQuad_IntLoop_L</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>fda29497e04c4d31160db111d4f3b8e3</anchor>
      <arglist>(int c, int i, int j, int type, short *S, int **ggg, int maxdist, int *p, int *q, paramT *P)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>inverse.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>inverse_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>inverse_fold</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>geabaaebf882055ba2be4ee79bb9f7b68</anchor>
      <arglist>(char *start, const char *target)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>inverse_pf_fold</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>g5a6b5e7e690d82eb9f76fcdca4c9ba35</anchor>
      <arglist>(char *start, const char *target)</arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>symbolset</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>gac1b233b632dc5493d7db6faaabacea2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>final_cost</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>g324d08a838c433d88c61274296f2e430</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>give_up</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>g42a09fe340bc850eceecb1d250b82c4c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>inv_verbose</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>ge080aad5ddc389b259dc8ee473bc961e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Lfold.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>Lfold_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>Lfold</name>
      <anchorfile>group__local__mfe__fold.html</anchorfile>
      <anchor>g0e2866cb1e62e7dd232401f30150e5f2</anchor>
      <arglist>(const char *string, char *structure, int maxdist)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>Lfoldz</name>
      <anchorfile>group__local__mfe__fold.html</anchorfile>
      <anchor>g533fba4cfbbff061762d922fcffa55b4</anchor>
      <arglist>(const char *string, char *structure, int maxdist, int zsc, double min_z)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>aliLfold</name>
      <anchorfile>group__local__consensus__fold.html</anchorfile>
      <anchor>ged2aeb7e512657929eb777e791a735f7</anchor>
      <arglist>(const char **strings, char *structure, int maxdist)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>loop_energies.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>loop__energies_8h</filename>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>E_IntLoop</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>1465eff5af77b784f95efc54c373e84d</anchor>
      <arglist>(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>E_Hairpin</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>80003bd877fa7cb8536a3f35f9ecbefc</anchor>
      <arglist>(int size, int type, int si1, int sj1, const char *string, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>E_Stem</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>da87d5360451964827082265b2d8bef4</anchor>
      <arglist>(int type, int si1, int sj1, int extLoop, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE double</type>
      <name>exp_E_Stem</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>626bd3b6e7be3d16431e0841d498e4b9</anchor>
      <arglist>(int type, int si1, int sj1, int extLoop, pf_paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE double</type>
      <name>exp_E_Hairpin</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>fb36ead5c01af61dd4e9872774c835f6</anchor>
      <arglist>(int u, int type, short si1, short sj1, const char *string, pf_paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE double</type>
      <name>exp_E_IntLoop</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>0513d1a53231ebaa0e1c4a2a9dd5d600</anchor>
      <arglist>(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>LPfold.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>LPfold_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>update_pf_paramsLP</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g7f109477f6756d03dbc780102b4079b2</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>pfl_fold</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>gb89faedaac434b9d540881f9d83c6bf9</anchor>
      <arglist>(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>pfl_fold_par</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g3c654120b40e7e5024521009680b34e8</anchor>
      <arglist>(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putoutpU_prob</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g6c82279376586dfe303f89bf473b2cb5</anchor>
      <arglist>(double **pU, int length, int ulength, FILE *fp, int energies)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putoutpU_prob_bin</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g70499cf6be97458642c4d3941a37df6e</anchor>
      <arglist>(double **pU, int length, int ulength, FILE *fp, int energies)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_pf_foldLP</name>
      <anchorfile>LPfold_8h.html</anchorfile>
      <anchor>69efd43313fc5ccf582c6cbadc142ab2</anchor>
      <arglist>(int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>MEA.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>MEA_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>MEA</name>
      <anchorfile>MEA_8h.html</anchorfile>
      <anchor>d59407009f795f1afe3f64e50f3b9286</anchor>
      <arglist>(plist *p, char *structure, double gamma)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mm.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>mm_8h</filename>
  </compound>
  <compound kind="file">
    <name>naview.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>naview_8h</filename>
  </compound>
  <compound kind="file">
    <name>params.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>params_8h</filename>
    <member kind="function">
      <type>paramT *</type>
      <name>scale_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g7046ab25440b97f1bcf686368c74942e</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>paramT *</type>
      <name>get_scaled_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gb423e01f73b60aaf952a79fcf5200796</anchor>
      <arglist>(double temperature, model_detailsT md)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_scaled_pf_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g99791e2d86bdbe49978022054b44263d</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_boltzmann_factors</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gceda7aaccb2d4719f4821fe16b93c8d3</anchor>
      <arglist>(double temperature, double betaScale, model_detailsT md, double pf_scale)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_boltzmann_factor_copy</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g8977a1422372354759395a0eb603b8ea</anchor>
      <arglist>(pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_scaled_alipf_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g6f060963d3357405f724e10e447cb1e0</anchor>
      <arglist>(unsigned int n_seq)</arglist>
    </member>
    <member kind="function">
      <type>PUBLIC pf_paramT *</type>
      <name>get_boltzmann_factors_ali</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g3d89dad88ca4e64c7b137a6abc79526a</anchor>
      <arglist>(unsigned int n_seq, double temperature, double betaScale, model_detailsT md, double pf_scale)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>part_func.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>part__func_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>pf_fold_par</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gf5fd8cc57e35ddd713122c886f215e15</anchor>
      <arglist>(const char *sequence, char *structure, pf_paramT *parameters, int calculate_bppm, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>pf_fold</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g25c82afd90b779becea09444dbe39bcf</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>pf_circ_fold</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gc9027d62aebf1b8fa7c767fb221992be</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>pbacktrack</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>g2d831ab73b5b166053ae661bf2b0eeaa</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>pbacktrack_circ</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>g037f22eb5ea8de98ca02510bdaa9c893</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_pf_arrays</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gba5cac9edbc48e9d0e8c8fe83f4d0781</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_pf_params</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g72007173e035f9605597fd0b0f101984</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_pf_params_par</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g60d7c405d606c8b2756e0891b5b19d6e</anchor>
      <arglist>(int length, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_bppm</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g3a9f8c18a6043cfc8bbc30a8a37e3969</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>assign_plist_from_pr</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g4099fa99067719874f214e3d329fbeea</anchor>
      <arglist>(plist **pl, double *probs, int length, double cutoff)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_pf_arrays</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g537854c3a777887c94bb0bed2baa26b0</anchor>
      <arglist>(short **S_p, short **S1_p, char **ptype_p, double **qb_p, double **qm_p, double **q1k_p, double **qln_p)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_subseq_F</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>bb2b5ba650533b0d3fa9b8d6c34b454f</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>get_centroid_struct_pl</name>
      <anchorfile>group__centroid__fold.html</anchorfile>
      <anchor>g1a3ec2649a3c531d453181685f2936ac</anchor>
      <arglist>(int length, double *dist, plist *pl)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>get_centroid_struct_pr</name>
      <anchorfile>group__centroid__fold.html</anchorfile>
      <anchor>g19f026a960a09d083e0fe0445e4672bf</anchor>
      <arglist>(int length, double *dist, double *pr)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_distance</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ge6a7ca930ed3c671b4dbdd1bcd9717b0</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_distance_pr</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gb98e05fef617c0333942a089e93be018</anchor>
      <arglist>(int length, double *pr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>bppm_to_structure</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>8bc7fb2b347512e67f413e5327c8f1a2</anchor>
      <arglist>(char *structure, double *pr, unsigned int length)</arglist>
    </member>
    <member kind="function">
      <type>char</type>
      <name>bppm_symbol</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>0e7a82f456e4597fb8f6be3061985e6c</anchor>
      <arglist>(const float *x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_pf_fold</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>30e926483c83e17b65095c2d1217f8d8</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>centroid</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>03929e89d3848604c7a283ae6743eb31</anchor>
      <arglist>(int length, double *dist)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_dist</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>b83fc8dea4ec3add6b237b86e2ce7286</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>expLoopEnergy</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>d04dc5cd85ef6577bb54b796815f2f50</anchor>
      <arglist>(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>expHairpinEnergy</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>742a67171d96ce562e949f3ddb7690b5</anchor>
      <arglist>(int u, int type, short si1, short sj1, const char *string)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>st_back</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>g9e37de73e77ae3cd5542d188602e8b9c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>part_func_co.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>part__func__co_8h</filename>
    <member kind="function">
      <type>cofoldF</type>
      <name>co_pf_fold</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gea46dc8f7374b7b93d51cbae3e162691</anchor>
      <arglist>(char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>cofoldF</type>
      <name>co_pf_fold_par</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gd0d429a4951be49455dcb0d20dd58643</anchor>
      <arglist>(char *sequence, char *structure, pf_paramT *parameters, int calculate_bppm, int is_constrained)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_co_bppm</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g98fd6158a5edc4911278b50e0ef94cb1</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_co_pf_arrays</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g91678df4db58bb0fd63ab9b98990394d</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_co_pf_params</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>geba781005bdd5793d554cc9ca9c56ce6</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_co_pf_params_par</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g2571a9e830f3f37d8027396c6efa807c</anchor>
      <arglist>(int length, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_probabilities</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gfefe156be07bf6f984dc1e165b77d3ce</anchor>
      <arglist>(double FAB, double FEA, double FEB, struct plist *prAB, struct plist *prA, struct plist *prB, int Alength)</arglist>
    </member>
    <member kind="function">
      <type>ConcEnt *</type>
      <name>get_concentrations</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gd1ca7df9c139484587aefaee6d1a3ce8</anchor>
      <arglist>(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>get_plist</name>
      <anchorfile>part__func__co_8h.html</anchorfile>
      <anchor>9d5278afa65f1348efbc6453a435714f</anchor>
      <arglist>(struct plist *pl, int length, double cut_off)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_co_pf_fold</name>
      <anchorfile>part__func__co_8h.html</anchorfile>
      <anchor>7e048f282df0bc954ab9fa4313237543</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>mirnatog</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g879522a85196a68fc0296c3cc9b10bcc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>F_monomer</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g6229b16b4482944bebea9414a53e1d0c</anchor>
      <arglist>[2]</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>part_func_up.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>part__func__up_8h</filename>
    <member kind="function">
      <type>pu_contrib *</type>
      <name>pf_unstru</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>g89dd9589f7d1f6582d0f9c428b9166ab</anchor>
      <arglist>(char *sequence, int max_w)</arglist>
    </member>
    <member kind="function">
      <type>interact *</type>
      <name>pf_interact</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>g3592ae8158ff790597f99a174df14408</anchor>
      <arglist>(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int max_w, char *cstruc, int incr3, int incr5)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_interact</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>gb273eebc8db290c0634fd769ade1d965</anchor>
      <arglist>(interact *pin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_pu_contrib_struct</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>gdd95a8c66d7a48219bef7ecc031770cf</anchor>
      <arglist>(pu_contrib *pu)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>plot_layouts.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>plot__layouts_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_PLOT_TYPE_SIMPLE</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>419d8d4772b14f63d1f4764d3fcdf8cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_PLOT_TYPE_NAVIEW</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>f1bb5bb724fa94cee207888118ecffde</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_PLOT_TYPE_CIRCULAR</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>87394f3cbf57ea0154faed2bc8c82bc5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>simple_xy_coordinates</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>7ce9669657183f9bdafa8d669bc3b7f3</anchor>
      <arglist>(short *pair_table, float *X, float *Y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>simple_circplot_coordinates</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>7e57e12dc7dfe9494bb6f807da3c9ef1</anchor>
      <arglist>(short *pair_table, float *x, float *y)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>rna_plot_type</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>40f24c96f6f5ebf22db60a3c820e24b4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>profiledist.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>profiledist_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>profile_edit_distance</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>e03b78cc482af7812f750c7f1cb2a40c</anchor>
      <arglist>(const float *T1, const float *T2)</arglist>
    </member>
    <member kind="function">
      <type>float *</type>
      <name>Make_bp_profile_bppm</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>b19c3d52d27182862a0da09150a3f852</anchor>
      <arglist>(double *bppm, int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_bppm</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>fb926c3587cdde9c4923925eb9b71199</anchor>
      <arglist>(const float *T)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_profile</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>22d18ea73ab4faad46872bffa2fcb3de</anchor>
      <arglist>(float *T)</arglist>
    </member>
    <member kind="function">
      <type>float *</type>
      <name>Make_bp_profile</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>f88f1864eb3bd9a53419decfbbb5e56c</anchor>
      <arglist>(int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>PS_dot.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>PS__dot_8h</filename>
    <member kind="function">
      <type>int</type>
      <name>PS_rna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>f7ce1bb2a3480d201c052c8dcac1d848</anchor>
      <arglist>(char *string, char *structure, char *file)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PS_rna_plot_a</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>c72680a16fe448908d3a6463b465b4c8</anchor>
      <arglist>(char *string, char *structure, char *file, char *pre, char *post)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>gmlRNA</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>1a2cab4e1ebf31b24261a25b6741253b</anchor>
      <arglist>(char *string, char *structure, char *ssfile, char option)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ssv_rna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>acfde04e76d83eee9071a182a864f9d5</anchor>
      <arglist>(char *string, char *structure, char *ssfile)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>svg_rna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>be977068d1eb713dc3153a977245cadb</anchor>
      <arglist>(char *string, char *structure, char *ssfile)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>xrna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>8b339233acda2304031ebafbad51e277</anchor>
      <arglist>(char *string, char *structure, char *ssfile)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PS_dot_plot_list</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>e19d5a95620133d9d0b9e686e8b7b21a</anchor>
      <arglist>(char *seq, char *filename, plist *pl, plist *mf, char *comment)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>aliPS_color_aln</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>f4e3fe44dea15be4374f66350f44824c</anchor>
      <arglist>(const char *structure, const char *filename, const char *seqs[], const char *names[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PS_dot_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>e39c043bd1a57731f0f7752a66cbf1d3</anchor>
      <arglist>(char *string, char *file)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>read_epars.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>read__epars_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>read_parameter_file</name>
      <anchorfile>group__energy__parameters__rw.html</anchorfile>
      <anchor>gc450f1b828dd7af183d3c134f25b124d</anchor>
      <arglist>(const char fname[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_parameter_file</name>
      <anchorfile>group__energy__parameters__rw.html</anchorfile>
      <anchor>g0312c701ca5326156ab6ec28c8fe8364</anchor>
      <arglist>(const char fname[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>RNAstruct.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>RNAstruct_8h</filename>
    <member kind="function">
      <type>char *</type>
      <name>b2HIT</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>55d8fb0508073cdabdf9aa0cc57447a5</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>b2C</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>d7cef34da9999bdc3d435f3abd8cdef1</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>b2Shapiro</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>26e751a615f4d2ae1e4551c7f876323f</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>add_root</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>dd70f6943687be2699c4ded3919db74b</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>expand_Shapiro</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>e1833510fb3f6430a0d7644161731eb9</anchor>
      <arglist>(const char *coarse)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>expand_Full</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>9ffcc30130c70ae29d780a8c1e5f6b25</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>unexpand_Full</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>ee45e1b11a35750d75f58eba6b730632</anchor>
      <arglist>(const char *ffull)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>unweight</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>7e2fe5dfd6f5956f3f340edf815ad258</anchor>
      <arglist>(const char *wcoarse)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unexpand_aligned_F</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>ad14cf08f537c6d25a636ab51a092470</anchor>
      <arglist>(char *align[2])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parse_structure</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>b621c5c9d5eb55ffc6f8a5ea03f71fe2</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>loop_size</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>709ca8d061b6093a57da17e948c88f69</anchor>
      <arglist>[STRUC]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>helix_size</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>973e5f8d5ba5652deb88ade43464e6d3</anchor>
      <arglist>[STRUC]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>loop_degree</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>bce4bda4d865eb4f6873efbc2815ab34</anchor>
      <arglist>[STRUC]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>loops</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>f0e57656a58f6928ccc1c73bf07561bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>unpaired</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>f47d2f021af446e98d3f51b6faa7db01</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>pairs</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>e2fa6b8e55d7964f8f3fd4ff6ce45477</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>stringdist.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>stringdist_8h</filename>
    <member kind="function">
      <type>swString *</type>
      <name>Make_swString</name>
      <anchorfile>stringdist_8h.html</anchorfile>
      <anchor>03f514f5a70a7b15e7e89cc05514fb98</anchor>
      <arglist>(char *string)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>string_edit_distance</name>
      <anchorfile>stringdist_8h.html</anchorfile>
      <anchor>389803f56854a9247cd2ea2c61348a2d</anchor>
      <arglist>(swString *T1, swString *T2)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>subopt.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>subopt_8h</filename>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g31909e554f5083c4f59cf2b771fed238</anchor>
      <arglist>(char *seq, char *structure, int delta, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt_par</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g1f43921d8b66091434e458b8c556ac7f</anchor>
      <arglist>(char *seq, char *structure, paramT *parameters, int delta, int is_constrained, int is_circular, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt_circ</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g8851d3250b1abc95b8734795e459f5c5</anchor>
      <arglist>(char *seq, char *sequence, int delta, FILE *fp)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>subopt_sorted</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>gc6a9413cc2b66a5ed9f0a20756772ec1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>print_energy</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g4dd8be3fb1adaa95baf6a7449aec1f68</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>density_of_states</name>
      <anchorfile>group__dos.html</anchorfile>
      <anchor>g18bf4f164623970156f3ac9f2986f437</anchor>
      <arglist>[MAXDOS+1]</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>treedist.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>treedist_8h</filename>
    <member kind="function">
      <type>Tree *</type>
      <name>make_tree</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>f364f296fa9a725b009ffe85487c7f3f</anchor>
      <arglist>(char *struc)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>tree_edit_distance</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>18dc4c201254d05aff231de759a9f7c4</anchor>
      <arglist>(Tree *T1, Tree *T2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tree</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>80ac708cdd15da644346d5190bdc5d11</anchor>
      <arglist>(Tree *t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_tree</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>a6583700fc10496e45b48a2b246702a9</anchor>
      <arglist>(Tree *t)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>utils.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/</path>
    <filename>utils_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_ERROR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>e6c620ed43bd41a2a96c2e54cd95adeb</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_QUIT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>95be86637d235cd3476bf77d69cceff7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_MISC</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>fb3db7eb518a29c4ed59e811644e6197</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_FASTA_HEADER</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>64458d0477ccb293b4b5ceac29664825</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_SEQUENCE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>6a627616dddc5cb4b66b636b2610db26</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_CONSTRAINT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>031ece0f176323b711e92aed2914b3b4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NO_TRUNCATION</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>c4ce83926b02fb7e70d27f431d8422e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NO_REST</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>47b0bb814624673f033a596a85bfa0fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NO_SPAN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>c0d8625fc8d82096cbc71ec066b118c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NOSKIP_BLANK_LINES</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>021b4d981b488ca57f0e29859f402a67</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_BLANK_LINE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>b3d6d56347592ce5fd1530b154c3c1a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NOSKIP_COMMENTS</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>9adc80b99094860389cac919dc486b72</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_COMMENT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>32a57751523205c2b25e72ad385332c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_PIPE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>89f0aba598c99b34edcabb75265f574e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_DOT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>cf452a44749d65687e085b14a85397ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_X</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>f0b50c1eb1147c9883ae5f4a7e85d7c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_ANG_BRACK</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>e20845c118ca551f2b2b7ad03590dc5a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_RND_BRACK</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ff1a166c7fad52d1075d73ef1888fe6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_MULTILINE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aaed46fca3c1e8a44bad466b04bc0b3a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_NO_HEADER</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a5ef5ebdeb4793654a8535a0c9cf35aa</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_ALL</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>cf41822e812dc8215471e6de46126f57</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_G</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>4f68a5cacc65d2a456e710be2a975c3c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_OPTION_MULTILINE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>763262aaa38330aa748670bbd8e2f632</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MIN2</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>5bef305e85be632c82b08bf537980405</anchor>
      <arglist>(A, B)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAX2</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>f2f40e921a3a9b6c1561e279ec1efcf2</anchor>
      <arglist>(A, B)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MIN3</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>f13418fa243155f44069849d605cfc94</anchor>
      <arglist>(A, B, C)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAX3</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>1cdcf98f8cce7ad1fdec6e7b637de6f6</anchor>
      <arglist>(A, B, C)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>XSTR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>f91a7287bf813eb35774869d1eb7165e</anchor>
      <arglist>(s)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>STR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>39ad1cc4f304d3f70fb2ffa037998942</anchor>
      <arglist>(s)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>FILENAME_MAX_LENGTH</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>142d29359a127a2528b46afe773d16cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>FILENAME_ID_LENGTH</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>0025e374ed97a3c2e823a68579517f15</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>space</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>7261b88efae2d25e6121db2f28e86227</anchor>
      <arglist>(unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>xrealloc</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>d73da84e51499eea58b5d8b99996567d</anchor>
      <arglist>(void *p, unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>nrerror</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>0b4244d16c7814c329c79f78695829f8</anchor>
      <arglist>(const char message[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>warn_user</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a2b2c3e83d76b53455b7b3701b9f2fe9</anchor>
      <arglist>(const char message[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_rand</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>8dfeabde9f7b0a332554090d46073987</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>urn</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>80338abe87350131e92d705cada17655</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>int_urn</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>b17b6f164204698a8d126445e0a45c5d</anchor>
      <arglist>(int from, int to)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>time_stamp</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>0a3e429ef54025ba1bfd246e97eeae25</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>random_string</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>95cfb010006f3ca6b6bcb4dfe8b50059</anchor>
      <arglist>(int l, const char symbols[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>hamming</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>fcf5e47ee5c79af67af8992b7a251223</anchor>
      <arglist>(const char *s1, const char *s2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>hamming_bound</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>d16f32cd9d4fff5484b617542798c130</anchor>
      <arglist>(const char *s1, const char *s2, int n)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>get_line</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>534a8f548b75e7aa7189ae5c50f45ef4</anchor>
      <arglist>(FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_input_line</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>6070ebdaf319eda73d9ff617594436d3</anchor>
      <arglist>(char **string, unsigned int options)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>read_record</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>636f2b61332f86184f72acedf276d1c1</anchor>
      <arglist>(char **header, char **sequence, char ***rest, unsigned int options)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>pack_structure</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>512a80f073ad7dc3da7a640102496579</anchor>
      <arglist>(const char *struc)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>unpack_structure</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>066fddaf0f6cb0e0654aa7b4465e40d2</anchor>
      <arglist>(const char *packed)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>make_pair_table</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>f405ebc979e35f06a2834333c5d2532c</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>copy_pair_table</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>7864eec7dcf8be8dd0baf738c98ccb30</anchor>
      <arglist>(const short *pt)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>alimake_pair_table</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>5f89ae4d382b4f598db3c13367576b36</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>make_pair_table_snoop</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>d8dbc0de4bb04f0bfbc5213856e833f3</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>make_loop_index_pt</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a891517b86c9dae0aa73ee09db6bfee1</anchor>
      <arglist>(short *pt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_input_seq</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>acf2afdab7485c35a3295de0d8c0dc62</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_input_seq_str</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>3e798798d0f3431466ef0ef48029130d</anchor>
      <arglist>(const char *s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_constraint_full</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>163cab70a4be347d76df004e63c3e3a9</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_constraint</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>d4379c5248591e903c2d0061b610e4bc</anchor>
      <arglist>(unsigned int option)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>str_DNA2RNA</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>fbe81f36a52574b41919655f071aebde</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>str_uppercase</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>dabc7ac69c1ccbf96ee859b78c80592f</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>get_iindx</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>361bacacd3717eb1b96c0a8dbae2184d</anchor>
      <arglist>(unsigned int length)</arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>get_indx</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>89ddf44c570bb2b0893aa76525275fdd</anchor>
      <arglist>(unsigned int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>constrain_ptypes</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>9a4d418db794bd36f4f1681f0ae2b3b0</anchor>
      <arglist>(const char *constraint, unsigned int length, char *ptype, int *BP, int min_loop_size, unsigned int idx_type)</arglist>
    </member>
    <member kind="variable">
      <type>unsigned short</type>
      <name>xsubi</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>2b948e2a79c281706e88b57cca6e37f6</anchor>
      <arglist>[3]</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>1.8.4_epars.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/lib/</path>
    <filename>1_88_84__epars_8h</filename>
  </compound>
  <compound kind="file">
    <name>1.8.4_intloops.h</name>
    <path>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/lib/</path>
    <filename>1_88_84__intloops_8h</filename>
  </compound>
  <compound kind="page">
    <name>mp_parse</name>
    <title>Parsing and Comparing - Functions to Manipulate Structures</title>
    <filename>mp_parse</filename>
  </compound>
  <compound kind="page">
    <name>mp_example</name>
    <title>Example - A Small Example Program</title>
    <filename>mp_example</filename>
    <docanchor file="mp_example">utils_struc</docanchor>
    <docanchor file="mp_example">utils_dot</docanchor>
    <docanchor file="mp_example">toc</docanchor>
    <docanchor file="mp_example">utils_misc</docanchor>
    <docanchor file="mp_example">utils_ss</docanchor>
    <docanchor file="mp_example">utils_seq</docanchor>
    <docanchor file="mp_example">utils_aln</docanchor>
  </compound>
  <compound kind="group">
    <name>folding_routines</name>
    <title>RNA Secondary Structure Folding</title>
    <filename>group__folding__routines.html</filename>
    <subgroup>mfe_fold</subgroup>
    <subgroup>pf_fold</subgroup>
    <subgroup>subopt_fold</subgroup>
    <subgroup>cofold</subgroup>
    <subgroup>consensus_fold</subgroup>
    <subgroup>local_fold</subgroup>
    <subgroup>energy_parameters</subgroup>
    <subgroup>eval</subgroup>
    <subgroup>inverse_fold</subgroup>
    <subgroup>class_fold</subgroup>
    <subgroup>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/fold.h</subgroup>
    <subgroup>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/subopt.h</subgroup>
  </compound>
  <compound kind="group">
    <name>mfe_fold</name>
    <title>Calculating Minimum Free Energy (MFE) Structures</title>
    <filename>group__mfe__fold.html</filename>
    <subgroup>mfe_cofold</subgroup>
    <subgroup>consensus_mfe_fold</subgroup>
    <subgroup>local_mfe_fold</subgroup>
    <subgroup>kl_neighborhood_mfe</subgroup>
    <member kind="function">
      <type>float</type>
      <name>fold_par</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gf3d4bf3ba3c37140d643203b67112495</anchor>
      <arglist>(const char *sequence, char *structure, paramT *parameters, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>fold</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gd11bb36c69e7e75f0bbae81a6642a015</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>circfold</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>g85be769330de2b504a1c130fe55af166</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_arrays</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gb823f35d01c03e9c10823c0fd002c8f6</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_fold_params</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>g44d079d97f158f294c61acf250f22582</anchor>
      <arglist>(void)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>pf_fold</name>
    <title>Calculating Partition Functions and Pair Probabilities</title>
    <filename>group__pf__fold.html</filename>
    <subgroup>mea_fold</subgroup>
    <subgroup>centroid_fold</subgroup>
    <subgroup>pf_cofold</subgroup>
    <subgroup>up_cofold</subgroup>
    <subgroup>consensus_pf_fold</subgroup>
    <subgroup>local_pf_fold</subgroup>
    <subgroup>kl_neighborhood_pf</subgroup>
    <member kind="function">
      <type>float</type>
      <name>pf_fold_par</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gf5fd8cc57e35ddd713122c886f215e15</anchor>
      <arglist>(const char *sequence, char *structure, pf_paramT *parameters, int calculate_bppm, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>pf_fold</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g25c82afd90b779becea09444dbe39bcf</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>pf_circ_fold</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gc9027d62aebf1b8fa7c767fb221992be</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_pf_arrays</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gba5cac9edbc48e9d0e8c8fe83f4d0781</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_pf_params</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g72007173e035f9605597fd0b0f101984</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_pf_params_par</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g60d7c405d606c8b2756e0891b5b19d6e</anchor>
      <arglist>(int length, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_bppm</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g3a9f8c18a6043cfc8bbc30a8a37e3969</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>assign_plist_from_pr</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g4099fa99067719874f214e3d329fbeea</anchor>
      <arglist>(plist **pl, double *probs, int length, double cutoff)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_pf_arrays</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>g537854c3a777887c94bb0bed2baa26b0</anchor>
      <arglist>(short **S_p, short **S1_p, char **ptype_p, double **qb_p, double **qm_p, double **q1k_p, double **qln_p)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_distance</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ge6a7ca930ed3c671b4dbdd1bcd9717b0</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_distance_pr</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gb98e05fef617c0333942a089e93be018</anchor>
      <arglist>(int length, double *pr)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>mea_fold</name>
    <title>Compute the structure with maximum expected accuracy (MEA)</title>
    <filename>group__mea__fold.html</filename>
  </compound>
  <compound kind="group">
    <name>centroid_fold</name>
    <title>Compute the centroid structure</title>
    <filename>group__centroid__fold.html</filename>
    <member kind="function">
      <type>char *</type>
      <name>get_centroid_struct_pl</name>
      <anchorfile>group__centroid__fold.html</anchorfile>
      <anchor>g1a3ec2649a3c531d453181685f2936ac</anchor>
      <arglist>(int length, double *dist, plist *pl)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>get_centroid_struct_pr</name>
      <anchorfile>group__centroid__fold.html</anchorfile>
      <anchor>g19f026a960a09d083e0fe0445e4672bf</anchor>
      <arglist>(int length, double *dist, double *pr)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>subopt_fold</name>
    <title>Enumerating Suboptimal Structures</title>
    <filename>group__subopt__fold.html</filename>
    <subgroup>subopt_zuker</subgroup>
    <subgroup>subopt_wuchty</subgroup>
    <subgroup>subopt_stochbt</subgroup>
  </compound>
  <compound kind="group">
    <name>subopt_zuker</name>
    <title>Suboptimal structures according to Zuker et al. 1989</title>
    <filename>group__subopt__zuker.html</filename>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>zukersubopt</name>
      <anchorfile>group__subopt__zuker.html</anchorfile>
      <anchor>gec07305955df39b1d579611c0039ddb9</anchor>
      <arglist>(const char *string)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>zukersubopt_par</name>
      <anchorfile>group__subopt__zuker.html</anchorfile>
      <anchor>g38aa4496326af32ea9ec137a9bf71ffe</anchor>
      <arglist>(const char *string, paramT *parameters)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>subopt_wuchty</name>
    <title>Suboptimal structures within an energy band arround the MFE</title>
    <filename>group__subopt__wuchty.html</filename>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g31909e554f5083c4f59cf2b771fed238</anchor>
      <arglist>(char *seq, char *structure, int delta, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt_par</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g1f43921d8b66091434e458b8c556ac7f</anchor>
      <arglist>(char *seq, char *structure, paramT *parameters, int delta, int is_constrained, int is_circular, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt_circ</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g8851d3250b1abc95b8734795e459f5c5</anchor>
      <arglist>(char *seq, char *sequence, int delta, FILE *fp)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>subopt_sorted</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>gc6a9413cc2b66a5ed9f0a20756772ec1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>print_energy</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>g4dd8be3fb1adaa95baf6a7449aec1f68</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>subopt_stochbt</name>
    <title>Stochastic backtracking in the Ensemble</title>
    <filename>group__subopt__stochbt.html</filename>
    <subgroup>consensus_stochbt</subgroup>
    <subgroup>kl_neighborhood_stochbt</subgroup>
    <member kind="function">
      <type>char *</type>
      <name>pbacktrack</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>g2d831ab73b5b166053ae661bf2b0eeaa</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>pbacktrack_circ</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>g037f22eb5ea8de98ca02510bdaa9c893</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>st_back</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>g9e37de73e77ae3cd5542d188602e8b9c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>cofold</name>
    <title>Calculate Secondary Structures of two RNAs upon Dimerization</title>
    <filename>group__cofold.html</filename>
    <subgroup>mfe_cofold</subgroup>
    <subgroup>pf_cofold</subgroup>
    <subgroup>up_cofold</subgroup>
  </compound>
  <compound kind="group">
    <name>mfe_cofold</name>
    <title>MFE Structures of two hybridized Sequences</title>
    <filename>group__mfe__cofold.html</filename>
    <file>cofold.h</file>
    <member kind="function">
      <type>float</type>
      <name>cofold</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>g6de084332decd07c52ede3d20af76d95</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>cofold_par</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>g83d04c7fb08ac05039575183545e31df</anchor>
      <arglist>(const char *string, char *structure, paramT *parameters, int is_constrained)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_co_arrays</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>gf31159249b8aaadfb95aa3997d422541</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_cofold_params</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>g95d130f8c6f091db0d90816857001e21</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>export_cofold_arrays_gq</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>ga3206b80db2f82bbff8296f0673af991</anchor>
      <arglist>(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **fc_p, int **ggg_p, int **indx_p, char **ptype_p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>export_cofold_arrays</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>gcdfda1fefb12716b1bfb2a7cacd0e5dc</anchor>
      <arglist>(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **fc_p, int **indx_p, char **ptype_p)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>pf_cofold</name>
    <title>Partition Function for two hybridized Sequences</title>
    <filename>group__pf__cofold.html</filename>
    <file>part_func_co.h</file>
    <member kind="function">
      <type>cofoldF</type>
      <name>co_pf_fold</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gea46dc8f7374b7b93d51cbae3e162691</anchor>
      <arglist>(char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>cofoldF</type>
      <name>co_pf_fold_par</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gd0d429a4951be49455dcb0d20dd58643</anchor>
      <arglist>(char *sequence, char *structure, pf_paramT *parameters, int calculate_bppm, int is_constrained)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_co_bppm</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g98fd6158a5edc4911278b50e0ef94cb1</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_co_pf_arrays</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g91678df4db58bb0fd63ab9b98990394d</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_co_pf_params</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>geba781005bdd5793d554cc9ca9c56ce6</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_co_pf_params_par</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g2571a9e830f3f37d8027396c6efa807c</anchor>
      <arglist>(int length, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_probabilities</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gfefe156be07bf6f984dc1e165b77d3ce</anchor>
      <arglist>(double FAB, double FEA, double FEB, struct plist *prAB, struct plist *prA, struct plist *prB, int Alength)</arglist>
    </member>
    <member kind="function">
      <type>ConcEnt *</type>
      <name>get_concentrations</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gd1ca7df9c139484587aefaee6d1a3ce8</anchor>
      <arglist>(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>mirnatog</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g879522a85196a68fc0296c3cc9b10bcc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>F_monomer</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>g6229b16b4482944bebea9414a53e1d0c</anchor>
      <arglist>[2]</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>up_cofold</name>
    <title>Partition Function for two hybridized Sequences as a stepwise Process</title>
    <filename>group__up__cofold.html</filename>
    <file>part_func_up.h</file>
    <member kind="function">
      <type>pu_contrib *</type>
      <name>pf_unstru</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>g89dd9589f7d1f6582d0f9c428b9166ab</anchor>
      <arglist>(char *sequence, int max_w)</arglist>
    </member>
    <member kind="function">
      <type>interact *</type>
      <name>pf_interact</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>g3592ae8158ff790597f99a174df14408</anchor>
      <arglist>(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int max_w, char *cstruc, int incr3, int incr5)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_interact</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>gb273eebc8db290c0634fd769ade1d965</anchor>
      <arglist>(interact *pin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_pu_contrib_struct</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>gdd95a8c66d7a48219bef7ecc031770cf</anchor>
      <arglist>(pu_contrib *pu)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>consensus_fold</name>
    <title>Predicting Consensus Structures from Alignment(s)</title>
    <filename>group__consensus__fold.html</filename>
    <subgroup>consensus_mfe_fold</subgroup>
    <subgroup>consensus_pf_fold</subgroup>
    <subgroup>consensus_stochbt</subgroup>
    <subgroup>local_consensus_fold</subgroup>
    <member kind="function">
      <type>int</type>
      <name>get_mpi</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gcffae3b2d3a3f38b830e89488f2242c0</anchor>
      <arglist>(char *Alseq[], int n_seq, int length, int *mini)</arglist>
    </member>
    <member kind="function">
      <type>float **</type>
      <name>readribosum</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g50ef3f5b786c98b20169c7bb6dded8d3</anchor>
      <arglist>(char *name)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_alistruct</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g13133e191923ba173dc40c9358ed861f</anchor>
      <arglist>(const char **sequences, const char *structure, int n_seq, float *energy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>encode_ali_sequence</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g12ba60c4e7278aa0cead337071aebc49</anchor>
      <arglist>(const char *sequence, short *S, short *s5, short *s3, char *ss, unsigned short *as, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>alloc_sequence_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g393d56ec257e2ec4bcb26b27dd9df7c1</anchor>
      <arglist>(const char **sequences, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_sequence_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g160c03a96eb91c6dd35aabed49f215e8</anchor>
      <arglist>(unsigned int n_seq, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_alipf_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gd6f7cd96aeaa838e32ada56b350a22f3</anchor>
      <arglist>(short ***S_p, short ***S5_p, short ***S3_p, unsigned short ***a2s_p, char ***Ss_p, double **qb_p, double **qm_p, double **q1k_p, double **qln_p, short **pscore)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>cv_fact</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gd6fac190cfd3af07934dfbba7dfdf636</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>nc_fact</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>g74ff4bb41b72ed967db7dce00fa66edc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>consensus_mfe_fold</name>
    <title>MFE Consensus Structures for Sequence Alignment(s)</title>
    <filename>group__consensus__mfe__fold.html</filename>
    <member kind="function">
      <type>float</type>
      <name>alifold</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>g86c43a79871b43bbee379b7b8f39d564</anchor>
      <arglist>(const char **strings, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>circalifold</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>gc245f44507e3324a7cf04879a5a7f8c3</anchor>
      <arglist>(const char **strings, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_alifold_arrays</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>g2f8550298ee654e25e44c30bb242bbd1</anchor>
      <arglist>(void)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>consensus_pf_fold</name>
    <title>Partition Function and Base Pair Probabilities for Sequence Alignment(s)</title>
    <filename>group__consensus__pf__fold.html</filename>
    <member kind="function">
      <type>float</type>
      <name>alipf_fold_par</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>g4159d3d8bb0fbb5b0be41fd54ad69030</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl, pf_paramT *parameters, int calculate_bppm, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_fold</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gca4fc5a8337db6695f0fb9fd4f5f0530</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_circ_fold</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gd6ee46b5bd800138407994fd1efc2197</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_ali_bppm</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gc59b5cbb1810e9331bc91edfd388365a</anchor>
      <arglist>(void)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>consensus_stochbt</name>
    <title>Stochastic Backtracking of Consensus Structures from Sequence Alignment(s)</title>
    <filename>group__consensus__stochbt.html</filename>
    <member kind="function">
      <type>char *</type>
      <name>alipbacktrack</name>
      <anchorfile>group__consensus__stochbt.html</anchorfile>
      <anchor>g59ce07041a74911ca72ea5efa0bb3478</anchor>
      <arglist>(double *prob)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>local_fold</name>
    <title>Predicting Locally stable structures of large sequences</title>
    <filename>group__local__fold.html</filename>
    <subgroup>local_mfe_fold</subgroup>
    <subgroup>local_pf_fold</subgroup>
    <subgroup>local_consensus_fold</subgroup>
  </compound>
  <compound kind="group">
    <name>local_mfe_fold</name>
    <title>Local MFE structure Prediction and Z-scores</title>
    <filename>group__local__mfe__fold.html</filename>
    <member kind="function">
      <type>float</type>
      <name>Lfold</name>
      <anchorfile>group__local__mfe__fold.html</anchorfile>
      <anchor>g0e2866cb1e62e7dd232401f30150e5f2</anchor>
      <arglist>(const char *string, char *structure, int maxdist)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>Lfoldz</name>
      <anchorfile>group__local__mfe__fold.html</anchorfile>
      <anchor>g533fba4cfbbff061762d922fcffa55b4</anchor>
      <arglist>(const char *string, char *structure, int maxdist, int zsc, double min_z)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>local_pf_fold</name>
    <title>Partition functions for locally stable secondary structures</title>
    <filename>group__local__pf__fold.html</filename>
    <member kind="function">
      <type>void</type>
      <name>update_pf_paramsLP</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g7f109477f6756d03dbc780102b4079b2</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>pfl_fold</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>gb89faedaac434b9d540881f9d83c6bf9</anchor>
      <arglist>(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>pfl_fold_par</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g3c654120b40e7e5024521009680b34e8</anchor>
      <arglist>(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putoutpU_prob</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g6c82279376586dfe303f89bf473b2cb5</anchor>
      <arglist>(double **pU, int length, int ulength, FILE *fp, int energies)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putoutpU_prob_bin</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>g70499cf6be97458642c4d3941a37df6e</anchor>
      <arglist>(double **pU, int length, int ulength, FILE *fp, int energies)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>local_consensus_fold</name>
    <title>Local MFE consensus structures for Sequence Alignments</title>
    <filename>group__local__consensus__fold.html</filename>
    <member kind="function">
      <type>float</type>
      <name>aliLfold</name>
      <anchorfile>group__local__consensus__fold.html</anchorfile>
      <anchor>ged2aeb7e512657929eb777e791a735f7</anchor>
      <arglist>(const char **strings, char *structure, int maxdist)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>energy_parameters</name>
    <title>Change and Precalculate Energy Parameter Sets and Boltzmann Factors</title>
    <filename>group__energy__parameters.html</filename>
    <file>params.h</file>
    <subgroup>energy_parameters_rw</subgroup>
    <member kind="function">
      <type>paramT *</type>
      <name>scale_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g7046ab25440b97f1bcf686368c74942e</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>paramT *</type>
      <name>get_scaled_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gb423e01f73b60aaf952a79fcf5200796</anchor>
      <arglist>(double temperature, model_detailsT md)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_scaled_pf_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g99791e2d86bdbe49978022054b44263d</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_boltzmann_factors</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gceda7aaccb2d4719f4821fe16b93c8d3</anchor>
      <arglist>(double temperature, double betaScale, model_detailsT md, double pf_scale)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_boltzmann_factor_copy</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g8977a1422372354759395a0eb603b8ea</anchor>
      <arglist>(pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_scaled_alipf_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g6f060963d3357405f724e10e447cb1e0</anchor>
      <arglist>(unsigned int n_seq)</arglist>
    </member>
    <member kind="function">
      <type>PUBLIC pf_paramT *</type>
      <name>get_boltzmann_factors_ali</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>g3d89dad88ca4e64c7b137a6abc79526a</anchor>
      <arglist>(unsigned int n_seq, double temperature, double betaScale, model_detailsT md, double pf_scale)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>energy_parameters_rw</name>
    <title>Reading/Writing energy parameter sets from/to File</title>
    <filename>group__energy__parameters__rw.html</filename>
    <file>read_epars.h</file>
    <subgroup>energy_parameters_convert</subgroup>
    <member kind="function">
      <type>void</type>
      <name>read_parameter_file</name>
      <anchorfile>group__energy__parameters__rw.html</anchorfile>
      <anchor>gc450f1b828dd7af183d3c134f25b124d</anchor>
      <arglist>(const char fname[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_parameter_file</name>
      <anchorfile>group__energy__parameters__rw.html</anchorfile>
      <anchor>g0312c701ca5326156ab6ec28c8fe8364</anchor>
      <arglist>(const char fname[])</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>energy_parameters_convert</name>
    <title>Converting energy parameter files</title>
    <filename>group__energy__parameters__convert.html</filename>
    <file>convert_epars.h</file>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_ALL</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g80ca4b03c4e33f8224bd74cb65ea8f21</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g48f634f8a30d7a76e9f671216a8e0a41</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_STACK</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gcc7cee41b96c907b1661acce4a1bb20c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g758f179bbca896409b22df663fa72ae7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga7bc7f81c5188c4e9d7d72fefec8417a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT_1N</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gb3defd4125dac42a477b74ef967b416d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT_23</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g4812fbe4be82f8d3a34b08d5a809d868</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_MULTI</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g8c9ffd43ac2cf9dee4b5923c2367d5c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_EXT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g2bc3edd21d014b0891a567eaa6aff758</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DANGLE5</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gd3988d62d6ebe9e5d220aea78512c3e8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DANGLE3</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g6c3eab5261e2cd3fe6a30e51406bf9d2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_11</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g54502f72c1a16d08acd233d06242cc0f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_21</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga0b6db305805e42ca08488f8436e9536</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_22</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gfe4f9deeb34ccc2400203aed4f494685</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_BULGE</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g22004edd37d0efff575f28d0f79f363d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g8415ba389c7c6d0fe9621c395cfb9fe0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_ML</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g5356c521b724fad6b5272d22cb14adb9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MISC</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g08304cb51f65c09dae2fc4eeb5e6277e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_SPECIAL_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ge42ee1ec21bd62f632bb554faf336667</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_VANILLA</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g6a704a05e53ee38f37564ed36122134e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_NINIO</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga299c4d0a5ea47e7301a304d2d6aa9ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DUMP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g1aabc3c11923e5253d6f001cf19a8f95</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convert_parameter_file</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>g3855fb2de2943aed520512eab2f1c7c3</anchor>
      <arglist>(const char *iname, const char *oname, unsigned int options)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>eval</name>
    <title>Energy evaluation</title>
    <filename>group__eval.html</filename>
    <member kind="function">
      <type>float</type>
      <name>energy_of_structure</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g3532b3f697f33a33cf02c774b6f49c15</anchor>
      <arglist>(const char *string, const char *structure, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_struct_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>gfe8e6f53f3641396e718bb45269e09e6</anchor>
      <arglist>(const char *string, const char *structure, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_structure</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g164a4b025da054f7f4445eba3652fd4a</anchor>
      <arglist>(const char *string, const char *structure, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_struct_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g46602e8377fad2621b873654bfd976f1</anchor>
      <arglist>(const char *string, const char *structure, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_structure_pt</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>ge478bf8b55469a82b580ceb9d198849d</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_struct_pt_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g2d5a0f1e2797a05d09d1867a8c61554a</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>eos_debug</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>g05f1af959998612960ddcf1436fa1837</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>inverse_fold</name>
    <title>Searching Sequences for Predefined Structures</title>
    <filename>group__inverse__fold.html</filename>
    <file>inverse.h</file>
    <member kind="function">
      <type>float</type>
      <name>inverse_fold</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>geabaaebf882055ba2be4ee79bb9f7b68</anchor>
      <arglist>(char *start, const char *target)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>inverse_pf_fold</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>g5a6b5e7e690d82eb9f76fcdca4c9ba35</anchor>
      <arglist>(char *start, const char *target)</arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>symbolset</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>gac1b233b632dc5493d7db6faaabacea2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>final_cost</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>g324d08a838c433d88c61274296f2e430</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>give_up</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>g42a09fe340bc850eceecb1d250b82c4c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>inv_verbose</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>ge080aad5ddc389b259dc8ee473bc961e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>class_fold</name>
    <title>Classified Dynamic Programming</title>
    <filename>group__class__fold.html</filename>
    <subgroup>kl_neighborhood</subgroup>
    <subgroup>dos</subgroup>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood</name>
    <title>Distance based partitioning of the Secondary Structure Space</title>
    <filename>group__kl__neighborhood.html</filename>
    <subgroup>kl_neighborhood_mfe</subgroup>
    <subgroup>kl_neighborhood_pf</subgroup>
    <subgroup>kl_neighborhood_stochbt</subgroup>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood_mfe</name>
    <title>Calculating MFE representatives of a Distance Based Partitioning</title>
    <filename>group__kl__neighborhood__mfe.html</filename>
    <file>2Dfold.h</file>
    <member kind="function">
      <type>TwoDfold_vars *</type>
      <name>get_TwoDfold_variables</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>g490d7a00dc1a6188bb401438fd975dbe</anchor>
      <arglist>(const char *seq, const char *structure1, const char *structure2, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroy_TwoDfold_variables</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>g9837b781df6e7274d848c137fb438bb5</anchor>
      <arglist>(TwoDfold_vars *our_variables)</arglist>
    </member>
    <member kind="function">
      <type>TwoDfold_solution *</type>
      <name>TwoDfoldList</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>g131c9c01e5203a1d0f8ab14ab9eb8315</anchor>
      <arglist>(TwoDfold_vars *vars, int distance1, int distance2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDfold_backtrack_f5</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>gf453ac72055554d5057ddbb166a67999</anchor>
      <arglist>(unsigned int j, int k, int l, TwoDfold_vars *vars)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood_pf</name>
    <title>Calculate Partition Functions of a Distance Based Partitioning</title>
    <filename>group__kl__neighborhood__pf.html</filename>
    <file>2Dpfold.h</file>
    <member kind="function">
      <type>TwoDpfold_vars *</type>
      <name>get_TwoDpfold_variables</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>g674dad118a59dc002516fe02b6613912</anchor>
      <arglist>(const char *seq, const char *structure1, char *structure2, int circ)</arglist>
    </member>
    <member kind="function">
      <type>TwoDpfold_vars *</type>
      <name>get_TwoDpfold_variables_from_MFE</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>g1a04119ebae7e52b9f6e05be864cdb91</anchor>
      <arglist>(TwoDfold_vars *mfe_vars)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroy_TwoDpfold_variables</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>ge6a9242efcf43798acad404334bccad4</anchor>
      <arglist>(TwoDpfold_vars *vars)</arglist>
    </member>
    <member kind="function">
      <type>TwoDpfold_solution *</type>
      <name>TwoDpfoldList</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>g4b9ecb276d9ae3c1da9e5ebc04c61be0</anchor>
      <arglist>(TwoDpfold_vars *vars, int maxDistance1, int maxDistance2)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood_stochbt</name>
    <title>Stochastic Backtracking of Structures from Distance Based Partitioning</title>
    <filename>group__kl__neighborhood__stochbt.html</filename>
    <member kind="function">
      <type>char *</type>
      <name>TwoDpfold_pbacktrack</name>
      <anchorfile>group__kl__neighborhood__stochbt.html</anchorfile>
      <anchor>g0b1553be92cca053c3bb5aa1c2b2fd23</anchor>
      <arglist>(TwoDpfold_vars *vars, int d1, int d2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDpfold_pbacktrack5</name>
      <anchorfile>group__kl__neighborhood__stochbt.html</anchorfile>
      <anchor>g319e7822ddbe0a5a3c9a6568c2e2a33a</anchor>
      <arglist>(TwoDpfold_vars *vars, int d1, int d2, unsigned int length)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>dos</name>
    <title>Compute the Density of States</title>
    <filename>group__dos.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>density_of_states</name>
      <anchorfile>group__dos.html</anchorfile>
      <anchor>g18bf4f164623970156f3ac9f2986f437</anchor>
      <arglist>[MAXDOS+1]</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>parse</name>
    <title>Parsing and Comparing - Functions to Manipulate Structures</title>
    <filename>group__parse.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/2Dfold.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/2Dfold.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_22Dfold_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/2Dpfold.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/2Dpfold.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_22Dpfold_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/alifold.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/alifold.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2alifold_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/cofold.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/cofold.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2cofold_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/convert_epars.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/convert_epars.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2convert__epars_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/fold.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/fold.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2fold_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/inverse.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/inverse.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2inverse_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/Lfold.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/Lfold.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2Lfold_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/LPfold.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/LPfold.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2LPfold_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/params.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/params.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2params_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/part_func.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/part_func.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2part__func_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/part_func_co.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/part_func_co.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2part__func__co_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/part_func_up.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/part_func_up.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2part__func__up_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/read_epars.h</name>
    <title>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/read_epars.h</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2read__epars_8h.html</filename>
  </compound>
  <compound kind="group">
    <name>/scratch/2/miladim/sparse/installRep/ViennaRNA-2.1.2/H/subopt.h</name>
    <title>Enumerating Suboptimal Structures</title>
    <filename>group___2scratch_22_2miladim_2sparse_2installRep_2ViennaRNA-2_81_82_2H_2subopt_8h.html</filename>
  </compound>
  <compound kind="struct">
    <name>bondT</name>
    <filename>structbondT.html</filename>
  </compound>
  <compound kind="struct">
    <name>bondTEn</name>
    <filename>structbondTEn.html</filename>
  </compound>
  <compound kind="struct">
    <name>constrain</name>
    <filename>structconstrain.html</filename>
  </compound>
  <compound kind="struct">
    <name>COORDINATE</name>
    <filename>structCOORDINATE.html</filename>
  </compound>
  <compound kind="struct">
    <name>cpair</name>
    <filename>structcpair.html</filename>
  </compound>
  <compound kind="struct">
    <name>INTERVAL</name>
    <filename>structINTERVAL.html</filename>
  </compound>
  <compound kind="struct">
    <name>model_detailsT</name>
    <filename>structmodel__detailsT.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>dangles</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>1c7589cc408e93ac58a0cb690a1f4865</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>special_hp</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>e4d4db67a0145d70b75751ffc6d2c806</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noLP</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>ee27d2beefc1f2ae49c7ff671bb67962</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noGU</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>6ba8f2c32b78a2aa9fcce7d4e9a3df1c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noGUclosure</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>9957b67f4b639529b13e83d2e1880fca</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>logML</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>710fadca0325236747c2a9380b9f89d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>circ</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>b9d1c648f7a7ea1f5ce32e291dfba7b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>gquad</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>13e66d2b473781dbce41c52f7ffbce37</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>PAIR</name>
    <filename>structPAIR.html</filename>
  </compound>
  <compound kind="struct">
    <name>pair_info</name>
    <filename>structpair__info.html</filename>
    <member kind="variable">
      <type>unsigned</type>
      <name>i</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>9f2ef040f76a4834a5a494ed824d0c65</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned</type>
      <name>j</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>8eb41bf173e0e266c7e4da1f945c0498</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>p</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>af3c48c1b04d6c74f6a8b56129d53349</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>ent</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>fa2ce4e64d7255685b2ac3138b0340eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>short</type>
      <name>bp</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>d74e2926dcca5d4da9f0bc76587cfc2e</anchor>
      <arglist>[8]</arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>comp</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>238a46e209b2c05239d3eefc5ab5297b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>paramT</name>
    <filename>structparamT.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>temperature</name>
      <anchorfile>structparamT.html</anchorfile>
      <anchor>8c6f38b84e570417b9ab51d4e39e0d4e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>model_detailsT</type>
      <name>model_details</name>
      <anchorfile>structparamT.html</anchorfile>
      <anchor>7eb7d612ab288df899292b444be3c2bb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>pf_paramT</name>
    <filename>structpf__paramT.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>pf_scale</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>a382d5f2a3c8d2649756e6cf66929e49</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>temperature</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>ec0fbed8a1f03142a0822d06bf3da57d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>alpha</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>4fd69aa76ae72a02096f79adc7ff0527</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>model_detailsT</type>
      <name>model_details</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>216f1fa1d6cf5443c8845b02178df6c1</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>plist</name>
    <filename>structplist.html</filename>
  </compound>
  <compound kind="struct">
    <name>pu_contrib</name>
    <filename>structpu__contrib.html</filename>
    <member kind="variable">
      <type>double **</type>
      <name>H</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>ca0249c706a29de30691a88c3f323454</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>I</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>cd066700ec801ac42b8824fd602a5275</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>M</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>e29c878f854eaa644c4d283f6cc968f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>E</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>28c2f5659011c0029fecae31926bb6cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>length</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>bcff7c5902f42ea9448fd162c52c61ba</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>w</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>8f57fce7699e8422de1578eb0007f96e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>pu_out</name>
    <filename>structpu__out.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>len</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>599e5777f5c7d939b297c0e5b43cf417</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>u_vals</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>be2472ab482f7faf7fc31a1c6b4a9c03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>contribs</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>2530400532d13cd91b0549e35256091e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char **</type>
      <name>header</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>8dbeecdca124eeffdd2c56b904e9553d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>u_values</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>12b7cdfd95917ecfb85ccb7961c7da21</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sect</name>
    <filename>structsect.html</filename>
  </compound>
  <compound kind="struct">
    <name>SOLUTION</name>
    <filename>structSOLUTION.html</filename>
    <member kind="variable">
      <type>float</type>
      <name>energy</name>
      <anchorfile>structSOLUTION.html</anchorfile>
      <anchor>1e87633bd01f6e5b4142780b89d08f93</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>structure</name>
      <anchorfile>structSOLUTION.html</anchorfile>
      <anchor>2278baaeabb53027932b255452a9a3d7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TwoDfold_solution</name>
    <filename>structTwoDfold__solution.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>k</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>db18ac6f40243340d6cd2dee8349e777</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>l</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>87dd0e66be30b90aa84e88f4bd0bfc6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>en</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>8f2aab366c43cdbfbccf5829675bf0b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>s</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>195de64d19728613bf523b59031eb305</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TwoDfold_vars</name>
    <filename>structTwoDfold__vars.html</filename>
    <member kind="variable">
      <type>paramT *</type>
      <name>P</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>c25264d131f7ba1d8357c7dcaf813f75</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>do_backtrack</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>e5688e2289a5b2fc1c6cbca76b5649c0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>ptype</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>8682278e96b287af9796805cd82a1224</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>sequence</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>5dd94a828187c9132ae156697fee3323</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>short *</type>
      <name>S1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>c536790ed88070718c01fc9c69bde336</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>eea4f877be2b182cdc0b4dd79eb1c44e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD2</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>dbd36e27559b368fa92e5e180a04dcbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>cce1f441851d2dac059d585b8f11d2ea</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm2</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>42508480267f05f4d7cd011ed860f9a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>my_iindx</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>7f70ecfad22c555e686a201c0a40a0a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>c166ac3176c3f06633626614014f19da</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs2</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>e97cb71a77c4326e9c9d34268d8534be</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>bpdist</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>12c8334afc10d5318aaa17ef449297f2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TwoDpfold_solution</name>
    <filename>structTwoDpfold__solution.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>k</name>
      <anchorfile>structTwoDpfold__solution.html</anchorfile>
      <anchor>0015ca965033c2e93db2b29b423b6e2b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>l</name>
      <anchorfile>structTwoDpfold__solution.html</anchorfile>
      <anchor>e09ab415277eb4779b98b9208e5431ba</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>q</name>
      <anchorfile>structTwoDpfold__solution.html</anchorfile>
      <anchor>3bd11c5fc3f89c82d76dc02f8bea89c7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TwoDpfold_vars</name>
    <filename>structTwoDpfold__vars.html</filename>
    <member kind="variable">
      <type>char *</type>
      <name>ptype</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>5527d0790ccdc65c403703460f58bbbe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>sequence</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>ca21919d855ec911d40639cc9a02bd43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>short *</type>
      <name>S1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a69ca8bdfb0b6f72ed2574258acbf15c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>d232b1593acd705d1819912782c3099f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD2</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>4d9ee078e765936ff5d42b211ab3dcd8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>my_iindx</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>b4ce04b6fdfeb521baad78cd75b9f144</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>jindx</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>ac9e9f1db9f93cc84ccf45fb5a4c00c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>f29058d02e0ebca3a6d80752a02c00c0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs2</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>195de975b797aa53a24eb60fb04a2e0a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>bpdist</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>373d72528c408eeedc92bbf9805e3bfb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>65165b33c03582387472866beff21e34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm2</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>ca12b3381c0686b0b1c6780192432079</anchor>
      <arglist></arglist>
    </member>
  </compound>
</tagfile>
