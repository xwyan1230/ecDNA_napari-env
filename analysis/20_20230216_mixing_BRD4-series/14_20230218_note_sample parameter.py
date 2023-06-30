# for confoal DNAFISH

sample1 = 'DM-Ctrl_mix_mCh-BRD4'
figure_name = sample1
pos_threshold = 21000
neg_threshold = 16000
sample1_pos = 'DM H2B-mCherry BRD4ko'
sample1_neg = 'DM'
new_seg = 11000
bg_neg = 4945.8102602298695
bg_pos = 5961.066231416957

sample1 = 'DM-Ctrl_mix_mCh-Ctrl'
figure_name = sample1
pos_threshold = 15000
neg_threshold = 12000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM'
new_seg = 11000
bg_neg = 4538.9979142022785
bg_pos = 4604.3040743596775

sample1 = 'DM-Ctrl_mix_mCh-BRD2'
figure_name = sample1
pos_threshold = 20000
neg_threshold = 14500
sample1_pos = 'DM H2B-mCherry BRD2ko'
sample1_neg = 'DM'
new_seg = 12000

sample1 = 'DM-Ctrl_mix_mCh-AHCTF1'
figure_name = sample1
pos_threshold = 12000
neg_threshold = 9000
sample1_pos = 'DM H2B-mCherry AHCTF1ko'
sample1_neg = 'DM'
new_seg = 12000

sample1 = 'DM-Ctrl_mix_mCh-BRD1'
figure_name = sample1
pos_threshold = 12500
neg_threshold = 8500
sample1_pos = 'DM H2B-mCherry BRD1ko'
sample1_neg = 'DM'
new_seg = 10000

sample1 = 'DM-Ctrl_mix_mCh-POLR3D'
figure_name = sample1
pos_threshold = 20000
neg_threshold = 11000
sample1_pos = 'DM H2B-mCherry POLR3Dko'
sample1_neg = 'DM'
new_seg = 10000


sample1 = 'DM-Ctrl_mix_mCh-BRD3'
figure_name = sample1
pos_threshold = 20000
neg_threshold = 11000
sample1_pos = 'DM H2B-mCherry BRD3ko'
sample1_neg = 'DM'
new_seg = 11000

sample1 = 'DM-Ctrl_mix_mCh-Ctrl_forBRD3'
figure_name = sample1
pos_threshold = 20000
neg_threshold = 11000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM'
new_seg = 12000


# for confocal IF
# DM-Ctrl_newMYC
img_MYC = np.concatenate([np.zeros(shape=[10, 3144]), img_MYC], axis=0)[:3144, :3144]
