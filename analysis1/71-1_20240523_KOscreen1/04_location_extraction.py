import pandas as pd
from aspose.cells import Workbook

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B3'
total_fov = 16
workbook = Workbook("%s/xml/%s.xml" % (data_dir, sample))
workbook.save("%s/%s/04_%s_xml.txt" % (output_dir, sample, sample))
location = pd.read_csv("%s/%s/04_%s_xml.txt" % (output_dir, sample, sample), na_values=['.'], sep='\t')
pd_loc = pd.DataFrame(columns=['fov', 'x', 'y'])
for i in range(total_fov):
    x = location['Unnamed: %s' % (15 * i + 12)][3]
    y = location['Unnamed: %s' % (15 * i + 14)][3]
    pd_loc.loc[len(pd_loc.index)] = [i, x, y]
pd_loc.to_csv('%s/%s/04_%s_location.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")





