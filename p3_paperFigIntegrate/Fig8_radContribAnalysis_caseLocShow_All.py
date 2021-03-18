'''
Author       : LYC
Date         : 2021-01-24 23:17:56
LastEditTime : 2021-02-06 11:11:20
LastEditors  : LYC
Description  : 
FilePath     : /code/p3_paperFigIntegrate/Fig8_radContribAnalysis_caseLocShow_All.py
 
'''
from __future__ import print_function
import os
import numpy as np
from numpy import nan
import Ngl,Nio
import scipy.io as sio





# ---Start the graphics section 创建画板
wks_type = "eps"
figPath="/home/liuyc/Research/P02.Ts_change_research/figure/proj3_PaperFig/v0.0/Fig8_caseShow"
wks = Ngl.open_wks(wks_type, figPath)

#
#  Create the plot.
#

#
#  Map resources.
#
res = Ngl.Resources()
res.nglFrame            = False   # Don't draw plot or advance frame
res.nglDraw             = False   # until we add shapefile outlines.

res.mpDataSetName         = "Earth..4"   # Change database
res.mpDataBaseVersion     = "MediumRes"  # Medium resolution database


res.mpFillOn                   = True
res.mpFillBoundarySets         = "NoBoundaries"
res.mpFillAreaSpecifiers       = ["land","water"]
res.mpSpecifiedFillColors      = ["gray","SlateGray2"]

# res.mpAreaMaskingOn            = True
# res.mpMaskAreaSpecifiers       = ["China"]
res.mpOutlineBoundarySets      = "National"
# res.mpOutlineSpecifiers        = ["China"]

# grid draw set
res.mpLimitMode           = "LatLon"        # Limit map via lat/lon
res.mpGridMaskMode        = "MaskFillArea"  # Draw no grid.
# res.mpOceanFillColor      = "SlateGray2"#-1
# res.mpGridLineDashPattern = 7

res.mpProjection          = "CylindricalEquidistant"
res.mpMinLatF             = -75
res.mpMaxLatF             =  75
res.mpMinLonF             =-150
res.mpMaxLonF             = 150

res.mpGeophysicalLineThicknessF = 2
res.mpNationalLineThicknessF = 2
res.mpCountyLineThicknessF = 3

#
# TickMark resources.
#
res.tmBorderThicknessF    = 5

# res.tmXBLabelFontAspectF = 2
res.tmXBLabelFontHeightF = 0.015
res.tmXBLabelFontThicknessF =100
res.tmXBMajorThicknessF = 5
res.tmXBMajorLengthF = 0.01

res.tmYLLabelFontThicknessF = 100
res.tmYLMajorThicknessF = 5

plot = Ngl.map(wks,res)                  # Create the map plot, but don't advance frame.

# 
# add specific out line
# 
# attribute of polyline
lnres = Ngl.Resources()
lnres.gsLineColor = "red"
lnres.gsLineThicknessF = 6.0            # 3x as thick
lnres.gsLineDashPattern = 0

# china east
lonCN=[112, 122];
latCN=[22, 37];
cnEast_bndry_lon=np.ravel([lonCN[0], lonCN[0], lonCN[1], lonCN[1], lonCN[0]]);
cnEast_bndry_lat=np.ravel([latCN[0], latCN[1], latCN[1], latCN[0], latCN[0]]);
lnid = Ngl.add_polyline(wks, plot, cnEast_bndry_lon, cnEast_bndry_lat, lnres)

# USA east 
lonUSA=[-90, -80];
latUSA=[30, 45];
USAEast_bndry_lon=np.ravel([lonUSA[0], lonUSA[0], lonUSA[1], lonUSA[1], lonUSA[0]]);
USAEast_bndry_lat=np.ravel([latUSA[0], latUSA[1], latUSA[1], latUSA[0], latUSA[0]]);
lnid = Ngl.add_polyline(wks, plot, USAEast_bndry_lon, USAEast_bndry_lat, lnres)

# EUR west
lonEUR=[-1, 14];
latEUR=[44, 54];
EURwest_bndry_lon=np.ravel([lonEUR[0], lonEUR[0], lonEUR[1], lonEUR[1], lonEUR[0]]);
EURwest_bndry_lat=np.ravel([latEUR[0], latEUR[1], latEUR[1], latEUR[0], latEUR[0]]);
lnid = Ngl.add_polyline(wks, plot, EURwest_bndry_lon, EURwest_bndry_lat, lnres)


Ngl.draw(plot)
Ngl.frame(wks) # 翻页操作
Ngl.destroy(wks)


Ngl.end()


'''
  res.mpDataSetName              = "./database/Earth..4"
  res.mpDataBaseVersion          = "MediumRes" # or "Ncarg4_1"
  
  res.mpFillOn                   = True
  res.mpFillBoundarySets         = "NoBoundaries"
  res.mpFillAreaSpecifiers       = ["land","water"]
  res.mpSpecifiedFillColors      = ["white","white"]
  
  res.mpAreaMaskingOn            = True
  res.mpMaskAreaSpecifiers       = ["China"]
  
  res.mpOutlineBoundarySets      = "NoBoundaries"
  res.mpOutlineSpecifiers        = ["China","China:Provinces"]
'''

