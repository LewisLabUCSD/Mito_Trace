from src.calculate_AF_by_cell import calculate_af

def test_calculate_AF_by_cell():
    ref_fa = "/data2/mito_lineage/BWA-Primers-MT/MT_genome/MT.fasta"
    maxBP=16571
    concat_dir="/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/scPileup_concat_200"
    coverage_dir="/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/CD34_mt_PolydT_scPileup_200"

    AF_F =  "{results}/{sample}/mapq_{mapq}/{sample}_scPileup_{num_read}/test_af.csv"
    calculate_af(coverage_dir, concat_dir, maxBP=maxBP, af_f=None, ref_fasta=ref_fa,
                 topn=500, min_cells=10,
                 min_reads=10, coverage_thresh=-1,
                 het_thresh=0.01, min_het_cells=10)
    return

test_calculate_AF_by_cell()