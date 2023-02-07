import os

for i in range(500):
    # os.system('python /Users/xwyan/Dropbox/LAB/github/2021-10_ecDNA/ecDNA_napari-env/analysis/15_20230119_model_clustering/05_sample_generation_arg.py -c 0 -cr 5 -r %s' % i)
    # os.system('python /Users/xwyan/Dropbox/LAB/github/2021-10_ecDNA/ecDNA_napari-env/analysis/15_20230119_model_clustering/13_sample_generation_different-r-nuclear.py -c 5000 -cr 5 -r %s' % i)
    os.system(
        'python /Users/xwyan/Dropbox/LAB/github/2021-10_ecDNA/ecDNA_napari-env/analysis/15_20230119_model_clustering/20_sample_generation_different-cp_fixed-cen_r_by-ac.py -c 5000 -cr 5 -cen_r 10 -r %s' % i)