This folder stores deprecated, intermediate, or reference scripts/files.

wrong_order_fastp_before_split/
    Scripts generated before confirming that snTAPS should follow split -> fastp.
    These outputs should not be used for formal analysis.

failed_test_scripts/
    Early shell-based read structure checking scripts. Replaced by Python version.

reference_scATAC_scripts/
    Original scATAC reference scripts. Kept only for provenance.

index_intermediate/
    Intermediate I5 files used during barcode dictionary normalization.
    Formal analysis uses 03_scripts/split/index/I5, I7, T5, T7.
