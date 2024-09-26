import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240325_create_SDC_position/"
interval = 350
z = 4050

"""file1.write('<?xml version="1.0" encoding="UTF-16"?>')
file1.write('<variant version="1.0">')
file1.write('<no_name runtype="CLxListVariant">')
file1.write('<bIncludeZ runtype="bool" value="true"/>')
file1.write('<bPFSEnabled runtype="bool" value="false"/>')"""

x_lst = np.arange(40500, -45000, -9000)
y_lst = np.arange(-22900, 25000, 9060)
xy_lst = []
for j in range(len(y_lst)):
    for i in range(len(x_lst)):
        if j % 2 == 0:
            xy_lst.append([x_lst[i], y_lst[j]])
        else:
            xy_lst.append([x_lst[len(x_lst)-i-1], y_lst[j]])

x_25pos_lst = [1.5, 0.5, -0.5, -1.5]
y_25pos_lst = [-1.5, -0.5, 0.5, 1.5]
xy_25pos_lst = []
position_lst = []
for k in range(len(xy_lst)):
    for j in range(len(y_25pos_lst)):
        for i in range(len(x_25pos_lst)):
            if j % 2 == 0:
                xy_25pos_lst.append([int(xy_lst[k][0]+x_25pos_lst[i]*interval), int(xy_lst[k][1]+y_25pos_lst[j]*interval)])
            else:
                xy_25pos_lst.append([int(xy_lst[k][0] + x_25pos_lst[len(x_25pos_lst)-i-1] * interval), int(xy_lst[k][1] + y_25pos_lst[j] * interval)])

for k in range(len(xy_25pos_lst)):
    if k < 10:
        position_lst.append('0000%s' % k)
    elif k < 100:
        position_lst.append('000%s' % k)
    elif k < 1000:
        position_lst.append('00%s' % k)
    elif k < 10000:
        position_lst.append('0%s' % k)

file1 = open("%s/test_960pos.txt" % master_folder, 'w')
for k in range(len(xy_25pos_lst)):
    print(k)
    file1.write('<Point%s runtype="NDSetupMultipointListItem">' % position_lst[k])
    file1.write('<bChecked runtype="bool" value="true"/>')
    file1.write('<strName runtype="CLxStringW" value=""/>')
    file1.write('<dXPosition runtype="double" value="%s"/>' % xy_25pos_lst[k][0])
    file1.write('<dYPosition runtype="double" value="%s"/>' % xy_25pos_lst[k][1])
    file1.write('<dZPosition runtype="double" value="%s"/>' % z)
    file1.write('<dPFSOffset runtype="double" value="-1.000000000000000"/>')
    file1.write('<baUserData runtype="CLxByteArray" value=""/>')
    file1.write('</Point%s>' % position_lst[k])

file1.close()
