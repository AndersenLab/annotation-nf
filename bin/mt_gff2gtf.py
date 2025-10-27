#!/usr/bin/env python

import sys


def main():
    gff_fname, gtf_fname = sys.argv[1:3]
    exons = {}
    cds = {}
    with open(gtf_fname, 'w') as output:
        for line in open(gff_fname):
            chrom, lab, biotype, start, stop, _, strand, frame, attributes = line.rstrip().split("\t")
            if chrom != "MtDNA":
                continue
            attributes = {field.split("=")[0]: field.split("=")[1] for field in attributes.split(';') if field != "NA"}
            gtf_info = {
                "gene_biotype": attributes["biotype"],
                "gene_name": attributes.get('sequence_name', ""),
            }
            if biotype == "gene":
                gtf_info["gene_id"] = attributes['ID'].split(":")[-1]
            elif biotype == "mRNA":
                gtf_info["gene_id"] = attributes['Parent'].split(":")[-1]
                biotype = "transcript"
                gtf_info['transcript_id'] = attributes['ID'].split(':')[-1]
                exons[gtf_info['transcript_id']] = [gtf_info['gene_id'], 0]
            elif biotype == "exon":
                gtf_info['transcript_id'] = attributes['Parent'].split(':')[-1]
                exons[gtf_info['transcript_id']][1] += 1
                gene_id, exon_num = exons[gtf_info['transcript_id']]
                gtf_info['gene_id'] = gene_id
                gtf_info['exon_number'] = str(exon_num)
            elif biotype in ["CDS", "three_prime_UTR", "five_prime_UTR", "stop_codon", "start_codon"]:
                gtf_info['transcript_id'] = attributes['Parent'].split(':')[-1]
                gene_id, exon_num = exons[gtf_info['transcript_id']]
                gtf_info['gene_id'] = gene_id
                gtf_info['exon_number'] = str(exon_num)
                if biotype == "CDS":
                    gtf_info['protein_id'] = attributes['Parent'].split(':')[-1]
            elif biotype == "tRNA":
                gtf_info["gene_id"] = attributes['Parent'].split(":")[-1]
                gtf_info['transcript_id'] = attributes['ID'].split(':')[-1]
                del gtf_info['gene_name']
            if not biotype in ["CDS", "stop_codon", "start_codon", "exon"]:
                frame = "."
            gtf_info = "; ".join([f'{k} "{v}"' for k, v in gtf_info.items()])
            output.write(f"{chrom}\t{lab}\t{biotype}\t{int(float(start))}\t{int(float(stop))}\t.\t{strand}\t{frame}\t{gtf_info};\n")

main()

            