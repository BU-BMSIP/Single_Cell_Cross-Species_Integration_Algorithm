python3 ../SATURN/train-saturn.py --in_data=data/zebrafish_frog_run.csv \
                              --in_label_col=cell_type --ref_label_col=cell_type \
                              --num_macrogenes=2000     --hv_genes=8000          \
                              --score_adata --ct_map_path=data/frog_zebrafish_cell_type_map.csv \
                              --work_dir=. \
                              --device_num=0
