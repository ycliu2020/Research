from __future__ import print_function
import os
import numpy as np
from numpy import nan
import Ngl
import Nio
import scipy.io as sio
import matlab
import matlab.engine
#  matlab engine start!
eng = matlab.engine.start_matlab()

# load global parm
for exmNum in [1,2]:
    # [readme, Experiment, level, tLin, mPlev, vars]
    glbParm = eng.cmipParameters(exmNum, nargout=6)
    mdlName = glbParm[2]['model2']
    for mdlNum in range(1,len(mdlName)+1): # len(mdlName)+1
        # mdlNum = 1
        
        # function [lon_f,lat_f,cc_global,headLineTxt,figPath,colorLab] = Pyn_figure1_2(exmNum, mdlNum)
        figData = eng.Pyn_figure4(exmNum, mdlNum, nargout=6)
        figPath = np.array(figData[4])
        figPath = str(figPath)
        
        if figPath=='/no':
            continue
        else: 
            # ---Generate some dummy lat/lon data
            lat = np.squeeze(np.array(figData[1], 'f'))
            lon = np.squeeze(np.array(figData[0], 'f'))

            nlat = lat.size
            nlon = lon.size

            # ---Start the graphics section 创建画板
            wks_type = "eps"
            wks = Ngl.open_wks(wks_type, figPath)


            # # One color map per row of plots
            colormap = np.array(figData[5])

            # ---Values to use for contour labelbar

            varMin = -0.8     # Data mins
            varMax = 0.8     # Data maxs
            varInt = 0.2     # Data spacing
            levels = np.arange(varMin, varMax+0.01, varInt)
            # list(np.arange(varMin,varMax,varInt))
            nlevs = len(levels)  # -- number of levels
            # -- convert list of floats to list of strings
            labels = ['{:.2f}'.format(x) for x in levels]


            # Create resource list for customizing contour over maps
            res = Ngl.Resources()

            res.nglMaximize = False
            res.nglFrame = False
            res.nglDraw = False

            res.sfXArray = lon
            res.sfYArray = lat
            res.sfMissingValueV = -99999


            # Title
            res.tiMainString = figData[3]
            res.tiMainFont = "Helvetica-Bold"
            res.tiMainFontHeightF = 0.015
            res.tiMainOffsetYF = 0.02

            # Map
            # res                             =  Ngl.Resources()
            res.nglDraw = False  # -- turn off plot draw and frame advance. We will
            res.nglFrame = False  # -- do it later after adding subtitles.

            res.mpProjection = "Mollweide"
            res.mpOutlineOn = True

            res.mpPerimOn = False
            # res.mpGridAndLimbDrawOrder="PreDraw"
            res.mpGridAndLimbOn = True
            res.mpLimbLineThicknessF = 2.
            res.mpGridLatSpacingF = 360               # spacing for lat lines
            res.mpGridLonSpacingF = 360               # spacing for lon lines
            res.mpGridLineColor = -1
            # res.mpOutlineOn                 = True
            # res.mpOutlineBoundarySets       = "National"

            # res.mpGridLineThicknessF = 5.
            res.mpGeophysicalLineThicknessF = 2.
            res.pmTickMarkDisplayMode = "Never"
            res.mpCenterLonF = 180
            # res.mpOceanFillColor = "white"
            # #-- create only a map
            # map = Ngl.map(wks,res)
            # Ngl.draw(map)

            # 等值线图的相关设置
            res.cnFillOn = True
            res.cnFillMode = 'CellFill'
            # res.cnCellFillMissingValEdgeColor = "gray50"  #-- missing value edges color
            # res.cnMissingValFillColor = "gray50"          #-- missing value fill color

            res.cnFillPalette = colormap     # Set the color map to use.

            res.cnLinesOn = False
            res.cnLineLabelsOn = False
            res.cnLevelSelectionMode = "ManualLevels"
            res.cnMinLevelValF = varMin
            res.cnMaxLevelValF = varMax
            res.cnLevelSpacingF = varInt

            # color bar set
            # res.lbLabelBarOn           = False
            res.lbOrientation = "Horizontal"   # Draw labelbar horizontally.
            res.lbLabelStrings = labels


            # find
            plotData = figData[2]
            plotData = np.transpose(np.squeeze(np.array(plotData)[:, :, 1]))

            # -- define _FillValue and missing_value if not existing
            missing = nan
            plotDataMask = np.ma.array(plotData, mask=np.equal(
                plotData, missing))  # -- mask array with missing values
            plotData[np.isnan(plotData)] = -99999

            plot = Ngl.contour_map(wks, plotData, res)

            # # Change some map resources and apply to existing contour/map plot.
            # res2 = Ngl.Resources()
            # res2.mpLandFillColor             = "transparent"  # Make sure land doesn't get filled again
            # res2.mpOceanFillColor            = "white"        # Fill water areas in white.
            # res2.mpInlandWaterFillColor      = "white"
            # Ngl.set_values(plot,res2)

            # # Change one resource to the contours in the existing plot, making them transparent.
            # res3 = Ngl.Resources()
            # res3.cnFillOpacityF = 0.0
            # Ngl.set_values(plot.contour,res3)

            Ngl.draw(plot)

            Ngl.frame(wks)
            Ngl.destroy(wks)
            # Ngl.end()
