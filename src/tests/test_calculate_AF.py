from src.calculate_AF import calculate_af

def test():
    ref_fa = "/data2/mito_lineage/BWA-Primers-MT/MT_genome/MT.fasta"
    maxBP=16571
    concat_wt_dir="/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/scPileup_concat_200"
    coverage_wt_dir="/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/CD34_mt_PolydT_scPileup_200"
    calculate_af(coverage_wt_dir, concat_wt_dir, ref_fasta=ref_fa,
                 AF_F=None, maxBP=maxBP, topN=500, min_cells=100,
                 min_reads=10)

    return
