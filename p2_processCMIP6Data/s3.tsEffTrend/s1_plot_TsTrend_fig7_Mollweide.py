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
for exmNum in [1,2,3,4]:  # [1,2,3,4]:
    # [readme, Experiment, level, tLin, mPlev, vars]
    glbParm = eng.cmipParameters(exmNum, nargout=6)
    mdlName = glbParm[2]['model2']
    for mdlNum in range(1,len(mdlName)+1):  # range(1,len(mdlName)+1): # len(mdlName)+1
        # mdlNum = 1

        # function [lon_f, lat_f, trendyr,figTitle, figPath, colorLab] = Pyn_figure3(exmNum, mdlNum)
        figData = eng.Pyn_figure7(exmNum, mdlNum, nargout=6)
        figPath = np.array(figData[4])
        figPath = str(figPath)

        if figPath == '/no':
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
            res.tiMainFont = "Helvetica-Bold"
            res.tiMainFontHeightF = 0.015
            # res.tiMainOffsetYF = 0.02

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
            res.lbLabelBarOn           = False
            res.lbOrientation = "Horizontal"   # Draw labelbar horizontally.
            res.lbLabelStrings = labels

            plotData = figData[2]  # 144 72 6
            # plotData = np.transpose(np.squeeze(np.array(plotData)[:, :, 1]))


            #
            # Loop 2 times and create 2 dummy plots;
            #
            nplots = 6
            plots = []
            for n in range(nplots):
                print("plot #{}".format(n))
                res.tiMainString = figData[3]['subTitle'][n] # 子图标题

                subplotData=np.transpose(np.squeeze(np.array(plotData)[:, :, n]))
                # -- define _FillValue and missing_value if not existing
                subplotData[np.isnan(subplotData)] = -99999

                plots.append(Ngl.contour_map(wks,subplotData,res))
  
            # Resources for panelling
            pres                  = Ngl.Resources() 
            pres.nglFrame         = False
            pres.nglPanelLabelBar = True

            # Calculate start Y position for first row of plots
            height = 0.25            # we know this will be height of small plots
            extra  = 1.0-(3*height)
            top    = 1.0-(extra/2.)  #子图第一行所在的高度

            # Draw a title before we draw plots
            title               = figData[3]['headLineTxt']
            txres               = Ngl.Resources()
            txres.txJust        = "BottomCenter"
            txres.txFontHeightF = 0.02
            Ngl.text_ndc(wks,title,0.5,top+0.01,txres)

            # Loop across plots and panel them on one page
            for n in range(0,1):
                # Define location in a unit square for each set of plots.
                pres.nglPanelTop    = top-(0*height)
                pres.nglPanelBottom = top-((2+1)*height)

                Ngl.panel(wks,plots[0:7],[3,2],pres)


            # plot = Ngl.contour_map(wks, plotData, res)

            # Ngl.draw(plot)

            Ngl.frame(wks)
            Ngl.destroy(wks)
            # Ngl.end()
