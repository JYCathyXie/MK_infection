import matplotlib.pyplot as plt
import numpy as np
from turtle import *
from PIL import Image
import cv2

# Read the defined map of MKs and boundaries
mk = cv2.imread('PATH\\mk1250x1250.png') # 1250 pix = 400 um
mk_gray = cv2.cvtColor(mk,cv2.COLOR_BGR2GRAY)
cv2.imshow('mk_gray',mk_gray)

# Define the index of background and CXCR4 low and high MKs
index_hi = np.argwhere(mk_gray == 0) # CXCR4 high MKs
index_lo = np.argwhere(mk_gray == 128) # CXCR4 low MKs
index_bg = np.argwhere(mk_gray == 255) # background - white
len(index_hi)
len(index_lo)
len(index_bg)

# Caculate the distance of a point and the closet target cell
def dist_min(point,target):
  if ((point == target).all(1).any()):
    dist_m = 0
    #print("YES",point)
  else:
    dist = np.sqrt(np.sum(np.asarray(target - point) ** 2, axis=1)) # distance of a point to targets
    dist_m = min(dist) 
  return dist_m

# Perform 500 times simulation
hi_perc_arr = []
lo_perc_arr = []
hi_mean_arr = []
lo_mean_arr = []
dist_data_arr = ()
for j in range(500):
  # Caculate min dist of points to CXCR4 low and high MKs
  dist_to_low = []
  dist_to_hig = []
  r = 2.5 # the radius of positioned myeloid cells
  hi_cnt = 0
  hi_arr = []
  lo_arr = []
  # Place randomly positioned myeloid cells on a 1250 pix x 1250 pix map
  n = 200 # 200 positioned cells
  def rand_data():
    return np.random.randint(low = 0, high = 1250, size = (n,))
  x1,y1 = [rand_data() for i in range(2)]
  # Caculate percentage of myeloid cells closer to CXCR4 high or CXCR4 low MKs, 
  # and mean dist of randomly placed cells to CXCR4 high or low MKs
  for i in range(n):
    pt = [x1[i],y1[i]]
    lo = dist_min(point = pt,target = index_lo) - r
    hi = dist_min(point = pt,target = index_hi) - r
    if (lo < 0):
      lo = 0
    if (hi < 0):
      hi = 0
    if (lo >= hi):
      hi_cnt += 1
      hi_arr.append(hi)
    else:
      lo_arr.append(lo)
    dist_to_low.append(lo)
    dist_to_hig.append(hi)
  print(j)
  hi_perc = hi_cnt/n*100 # percentage of myeloid cells closer to CXCR4 high
  lo_perc = 100 - hi_perc # percentage of myelodi cells closer to CXCR4 low
  hi_mean = np.mean(hi_arr) # mean dist to CXCR4 high
  lo_mean = np.mean(lo_arr) # mean dist to CXCR4 low
  hi_perc_arr.append(hi_perc)
  lo_perc_arr.append(lo_perc)
  hi_mean_arr.append(hi_mean)
  lo_mean_arr.append(lo_mean)
  dist_data = np.array([dist_to_low,dist_to_hig]) 
  dist_data_a = dist_data.transpose()
  dist_data_arr = dist_data_arr+(dist_data_a,)

np.mean(lo_perc_arr) # 63.27 % randomly placed cells are closer to CXCR4 low MKs
np.mean(hi_perc_arr) # 36.73 % randomly placed cells are closer to CXCR4 high MKs
np.mean(lo_mean_arr)/1250*400 # 1250 pix = 400 um
np.mean(hi_mean_arr)/1250*400 # 1250 pix = 400 um

# Output
np.savetxt('E:\\Docum\\data\\MK_immune\\Cd11b_MK_distribution\\random_distribution\\closer_high_percentage.txt',hi_perc_arr)
np.savetxt('E:\\Docum\\data\\MK_immune\\Cd11b_MK_distribution\\random_distribution\\closer_low_percentage.txt',lo_perc_arr)
np.savetxt('E:\\Docum\\data\\MK_immune\\Cd11b_MK_distribution\\random_distribution\\closer_high_mean_dist.txt',hi_mean_arr)
np.savetxt('E:\\Docum\\data\\MK_immune\\Cd11b_MK_distribution\\random_distribution\\closer_low_mean_dist.txt',lo_mean_arr)

exmp = dist_data_arr[4] # a presentative example of random distribution (63.3% cells closer to CXCR4 low MKs and 36.7% closer to CXCR4 high MKs) 
np.savetxt('PATH\\example_low_high_dist.txt',exmp)
