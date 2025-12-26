#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

doc: |
    This workflow will perform preprocessing steps on VCFs for the OxoG/Variantbam/Annotation workflow.

$namespaces:
 s: https://schema.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:author:
    s:name: "Solomon Shorser"
    s:email: "solomon.shorser@oicr.on.ca"

requirements:
    - class: SchemaDefRequirement
      types:
          - $import: PreprocessedFilesType.yaml
    - class: ScatterFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: MultipleInputFeatureRequirement
    - class: InlineJavascriptRequirement
      expressionLib:
        - { $include: ./preprocess_util.js }
        # Shouldn't have to *explicitly* include vcf_merge_util.js but there's
        # probably a bug somewhere that makes it necessary
        - { $include: ./vcf_merge_util.js }
    - class: SubworkflowFeatureRequirement

inputs:
    - id: vcfdir
      type: Directory
      doc: "The directory where the files are"
    - id: filesToPreprocess
      type: string[]
      doc: "The files to process"
    - id: ref
      type: File
      doc: "Reference file, used for normalized INDELs"
    - id: out_dir
      type: string
      doc: "The name of the output directory"

# There are three output sets:
# - The merged VCFs.
# - The VCFs that are cleaned and normalized.
# - The SNVs that were extracted from INDELs (if there were any - usually there are none).
outputs:
    preprocessedFiles:
        type: "PreprocessedFilesType.yaml#PreprocessedFileset"
        outputSource: populate_output_record/output_record

steps:
    # TODO: Exclude MUSE files from PASS-filtering. MUSE files still need to be cleaned, but should
    # not be PASS-filtered.
    pass_filter:
      doc: "Filter out non-PASS lines from the VCFs in filesToProcess."
      in:
        vcfdir: vcfdir
        filesToFilter:
            source: [ filesToPreprocess ]
            valueFrom: |
                ${
                    var VCFs = []
                    for (var i in self)
                    {
                        if (self[i].toLowerCase().indexOf("muse") == -1)
                        {
                            VCFs.push(self[i]);
                        }
                    }
                    return VCFs;
                }
      run: pass-filter.cwl
      out: [output]

    clean:
      doc: "Clean the VCFs."
      run: clean_vcf.cwl
      scatter: [vcf]
      in:
        vcf: pass_filter/output
      out: [clean_vcf]

    gather_muse_snvs_for_cleaning:
      in:
        vcfdir:
            source: vcfdir
        vcfs:
            source: filesToPreprocess
      out: [snvs_for_cleaning]
      run:
        class: ExpressionTool
        inputs:
          vcfs: string[]
          vcfdir: Directory
        outputs:
          snvs_for_cleaning: File[]
        expression: |
            $({
                snvs_for_cleaning: (  filterFor("MUSE","snv_mnv", inputs.vcfs ).map(function(e) {
                    e = { "class":"File", "location":inputs.vcfdir.location+"/"+e }
                    return e;
                } )  )
            })

    clean_muse:
      doc: "Clean the MUSE VCFs."
      run: clean_vcf.cwl
      scatter: [vcf]
      in:
        vcf: gather_muse_snvs_for_cleaning/snvs_for_cleaning
      out: [clean_vcf]


    filter_for_indel:
      doc: "Filters the input list and selects the INDEL VCFs."
      in:
        in_vcf: clean/clean_vcf
      out: [out_vcf]
      run:
        class: ExpressionTool
        inputs:
          in_vcf: File[]
        outputs:
          out_vcf: File[]
        expression: |
            $({ out_vcf: filterForIndels(inputs.in_vcf) })

    normalize:
      doc: "Normalize the INDEL VCFs."
      run: normalize.cwl
      scatter: normalize/vcf
      in:
        vcf:
          source: filter_for_indel/out_vcf
        ref: ref
      out: [normalized-vcf]

    extract_snv:
      doc: "Extract SNVs from normalized INDELs"
      run: extract_snv.cwl
      scatter: extract_snv/vcf
      in:
          vcf: normalize/normalized-vcf
      out: [ extracted_snvs ]

    # Remove "null" elements from the array.
    null_filter_extracted_snvs:
      in:
        extracted_snvs:
          source: extract_snv/extracted_snvs
        muse_files:
          source: clean_muse/clean_vcf
      run:
          class: ExpressionTool
          inputs:
              extracted_snvs:
                  type:
                      type: array
                      items: [ File, "null" ]
          outputs:
              cleaned_extracted_snvs: File[]
          expression: |
            $(
                { cleaned_extracted_snvs: inputs.extracted_snvs.filter(function(n){return n != null}) }
            )
      out: [cleaned_extracted_snvs]

    #############################################
    # Gather SNVs on a per-workflow basis
    #############################################

    gather_sanger_snvs:
      in:
        clean_vcfs:
            source: clean/clean_vcf
        extracted_snvs:
            source: null_filter_extracted_snvs/cleaned_extracted_snvs
      out: [snvs_for_merge]
      run:
        class: ExpressionTool
        inputs:
          clean_vcfs: File[]
          extracted_snvs: File[]
        outputs:
          snvs_for_merge: File[]
        expression: |
            $({
                snvs_for_merge: ( (filterFor("svcp","snv_mnv",inputs.clean_vcfs)).concat(filterFor("svcp","snv_mnv",inputs.extracted_snvs)) )
            })

    gather_dkfz_embl_snvs:
      in:
        clean_vcfs:
            source: clean/clean_vcf
        extracted_snvs:
            source: null_filter_extracted_snvs/cleaned_extracted_snvs
      out: [snvs_for_merge]
      run:
        class: ExpressionTool
        inputs:
          clean_vcfs: File[]
          extracted_snvs: File[]
        outputs:
          snvs_for_merge: File[]
        expression: |
            $({
                snvs_for_merge: ( (filterFor("dkfz-snvCalling","snv_mnv",inputs.clean_vcfs)).concat(filterFor("dkfz-snvCalling","snv_mnv",inputs.extracted_snvs)) )
            })

    gather_broad_snvs:
      in:
        clean_vcfs:
            source: clean/clean_vcf
        extracted_snvs:
            source: null_filter_extracted_snvs/cleaned_extracted_snvs
      out: [snvs_for_merge]
      run:
        class: ExpressionTool
        inputs:
          clean_vcfs: File[]
          extracted_snvs: File[]
        outputs:
          snvs_for_merge: File[]
        expression: |
            $({
                snvs_for_merge: ( (filterFor("broad-mutect","snv_mnv",inputs.clean_vcfs)).concat(filterFor("broad-mutect","snv_mnv",inputs.extracted_snvs)) )
            })

    gather_muse_snvs:
      in:
        clean_vcfs:
            source: clean_muse/clean_vcf
        extracted_snvs:
            # MUSE INDELs will not be provided so there will *never* be any extracted SNVs for MUSE.
            source: null_filter_extracted_snvs/cleaned_extracted_snvs
      out: [snvs_for_merge]
      run:
        class: ExpressionTool
        inputs:
          clean_vcfs: File[]
          extracted_snvs: File[]
        outputs:
          snvs_for_merge: File[]
        expression: |
            $({
                snvs_for_merge: ( filterFor("MUSE","snv_mnv",inputs.clean_vcfs)).concat(filterFor("MUSE","snv_mnv",inputs.extracted_snvs))
            })

    #############################################
    # Gather INDELs on a per-workflow basis
    #############################################

    gather_sanger_indels:
      in:
        normalized_vcfs: normalize/normalized-vcf
      out: [indels_for_merge]
      run:
        class: ExpressionTool
        inputs:
          normalized_vcfs: File[]
        outputs:
          indels_for_merge: File[]
        expression: |
            $({ indels_for_merge: filterFor("svcp","indel",inputs.normalized_vcfs) })

    gather_dkfz_embl_indels:
      in:
        normalized_vcfs: normalize/normalized-vcf
      out: [indels_for_merge]
      run:
        class: ExpressionTool
        inputs:
          normalized_vcfs: File[]
        outputs:
          indels_for_merge: File[]
        expression: |
            $({ indels_for_merge: filterFor("dkfz-indelCalling","indel",inputs.normalized_vcfs) })

    gather_broad_indels:
      in:
        normalized_vcfs: normalize/normalized-vcf
      out: [indels_for_merge]
      run:
        class: ExpressionTool
        inputs:
          normalized_vcfs: File[]
        outputs:
          indels_for_merge: File[]
        expression: |
            $({ indels_for_merge: filterFor("broad-snowman","indel",inputs.normalized_vcfs) })

    gather_smufin_indels:
      in:
        normalized_vcfs: normalize/normalized-vcf
      out: [indels_for_merge]
      run:
        class: ExpressionTool
        inputs:
          normalized_vcfs: File[]
        outputs:
          indels_for_merge: File[]
        expression: |
            $({ indels_for_merge: filterFor("smufin","indel",inputs.normalized_vcfs) })


    #############################################
    # Gather SVs on a per-workflow basis
    #############################################

    gather_sanger_svs:
      in:
        in_vcf: pass_filter/output
      out: [svs_for_merge]
      run:
        class: ExpressionTool
        inputs:
          in_vcf: File[]
        outputs:
          svs_for_merge: File[]
        expression: |
            $({ svs_for_merge: filterFor("svfix",".sv.",inputs.in_vcf) })

    gather_broad_svs:
      in:
        in_vcf: pass_filter/output
      out: [svs_for_merge]
      run:
        class: ExpressionTool
        inputs:
          in_vcf: File[]
        outputs:
          svs_for_merge: File[]
        expression: |
            $({ svs_for_merge: filterFor("broad-dRanger_snowman",".sv.",inputs.in_vcf) })

    gather_dkfz_embl_svs:
      in:
        in_vcf: pass_filter/output
      out: [svs_for_merge]
      run:
        class: ExpressionTool
        inputs:
          in_vcf: File[]
        outputs:
          svs_for_merge: File[]
        expression: |
            $({ svs_for_merge: filterFor("embl-delly",".sv.",inputs.in_vcf) })


    ###############################################
    # Do the VCF Merge
    ###############################################

    merge_vcfs:
      doc: "Perform VCF merge by type (SNV, SV, INDEL) across workflows."
      run:
        vcf_merge.cwl
      in:
        sanger_snv:
            source: gather_sanger_snvs/snvs_for_merge
        de_snv:
            source: gather_dkfz_embl_snvs/snvs_for_merge
        broad_snv:
            source: gather_broad_snvs/snvs_for_merge
        muse_snv:
            source: gather_muse_snvs/snvs_for_merge
        sanger_indel:
            source: gather_sanger_indels/indels_for_merge
        de_indel:
            source: gather_dkfz_embl_indels/indels_for_merge
        broad_indel:
            source: gather_broad_indels/indels_for_merge
        smufin_indel:
            source: gather_smufin_indels/indels_for_merge
        sanger_sv:
            source: gather_sanger_svs/svs_for_merge
        de_sv:
            source: gather_dkfz_embl_svs/svs_for_merge
        broad_sv:
            source: gather_broad_svs/svs_for_merge
        out_dir: out_dir
      out:
          [output]


    populate_output_record:
        in:
            mergedVcfs : merge_vcfs/output
            extractedSnvs : null_filter_extracted_snvs/cleaned_extracted_snvs
            normalizedVcfs: normalize/normalized-vcf
            cleanedVcfs: clean/clean_vcf
        out:
            [output_record]
        run:
            class: ExpressionTool
            inputs:
                mergedVcfs: File[]
            outputs:
              output_record: "PreprocessedFilesType.yaml#PreprocessedFileset"
            expression: |
                    $(
                        {output_record: {
                            "mergedVcfs": inputs.mergedVcfs,
                        }}
                    )
