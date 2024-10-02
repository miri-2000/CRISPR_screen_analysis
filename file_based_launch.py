from start_program import CRISPRScreenAnalysis


class Args:
    input_file = (r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Datasets to be analysed"
                  r"\dataset\input_data\7234_all_Brunello_library_target_genes-req.txt")
    essential_genes = (r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Program data"
                       r"\list_essentials_roderick.csv")
    non_essential_genes = (r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Program data"
                           r"\non_essentials_roderick.csv")
    library_file = (r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Program data"
                    r"\broadgpp-brunello-library-contents.txt")
    target_samples = "t2_2d_tr,t2_3d_tr,t1_3d,t1_3d2d,t2_3d2d_ut,t2_3d2d_tr,t2_2d_tr,t2_2d_ut,t2_3d_tr,t2_3d_ut"
    reference_samples = "t2_2d_ut,t2_3d_ut,t1_2d,t1_3d,t2_3d_ut,t2_3d_tr,t1_2d,t1_2d,t1_3d,t1_3d"
    unwanted_columns = "guide_mm1_mismatch1,mismatch1_,nohit_cols,guide_mm1_nohit"
    unwanted_rows = ""
    unwanted_row_substrings = ":mismatch"
    threshold_reads = 0
    x_axis = "normZ"
    threshold_fdr = 0.25
    top = 15
    distribution_condition1 = "t1_2d"
    distribution_condition2 = "t0"
    replicate_type = "biological"
    working_dir = r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test"


if __name__ == "__main__":
    CRISPRScreenAnalysis(Args)
