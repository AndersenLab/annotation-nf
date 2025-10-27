process BEDTOOLS_HDR_ANNOTATIONS {
    tag "$strain_name"

    input:
    tuple val(meta), path(strain_matrix)
    tuple val(meta1), path(no_mt_vcf_bed)
    path hdr
    val strain_name

    output:
    path "${strain_name}.geno", emit: geno
    path "versions.yml",        emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Subset of master HDR BED file to only contain HDRs for the strain that is being iterated on
    if [[ "${hdr}" == *.gz ]]; then
        zcat ${hdr} | grep -w ${strain_name} > temp_${strain_name}.bed || touch temp_${strain_name}.bed
    else
        grep -w ${strain_name} ${hdr} > temp_${strain_name}.bed || touch temp_${strain_name}.bed
    fi

    # Check if temp bed file has content
    if [[ -s temp_${strain_name}.bed ]]; then
        # Identify variants within HDRs for the specific strain being iterated on
        bedtools intersect -a ${no_mt_vcf_bed} -b temp_${strain_name}.bed -wa | cut -f1,2 > ${strain_name}_intersect.bed
        
        # Generate genotype data for each variant position
        awk -v intersect_bed="${strain_name}_intersect.bed" -v strain="${strain_name}" '
            BEGIN {
                # Read intersect bed file into array
                while ((getline line < intersect_bed) > 0) {
                    split(line, fields, "\\t")
                    arr[fields[1] "\\t" fields[2]] = 1
                }
                close(intersect_bed)
            }
            {
                if (strain == "N2") {
                    print 0
                } else if ((\$1 "\\t" \$2) in arr) {
                    print 1
                } else {
                    print 0
                }
            }' ${strain_matrix} > ${strain_name}.geno
        
        rm ${strain_name}_intersect.bed
    else
        echo "No HDR calls found for strain ${strain_name} - generating all zeros"
        
        # Count lines in strain matrix and generate zeros
        N=\$(wc -l < ${strain_matrix})
        for i in \$(seq 1 \${N}); do
            echo "0"
        done > ${strain_name}.geno
    fi

    # Clean up temporary file
    rm -f temp_${strain_name}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | cut -f2 -d" ")
    END_VERSIONS
    """

    stub:
    """
    touch ${strain_name}.geno

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | cut -f2 -d" ")
    END_VERSIONS
    """
}