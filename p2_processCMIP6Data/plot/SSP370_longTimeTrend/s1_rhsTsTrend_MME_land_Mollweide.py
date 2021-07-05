'''
Author       : LYC
Date         : 2021-05-11 14:31:40
LastEditors  : Please set LastEditors
Description  :
FilePath     : /code/p2_processCMIP6Data/plot/SSP370_longTimeTrend/s1_rhsTsTrend_MME_land_Mollweide.py
symbol_custom_string_obkoro1:
'''
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

MMEType = "MME1"
plotType=4
for exmNum in [4]:  # [1,2,3,4]:

    # glbParm = eng.cmipParameters(exmNum, nargout=6)

    figData = eng.ssp370_rhsTsTrend_MME(MMEType, exmNum, plotType, nargout=6)
    figPath = np.array(figData[4])
    figPath = str(figPath)
    
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
    if plotType==6:
        varMin = -0.5     # Data mins
    elif plotType==4:        
        varMin = -1.5     # Data mins

    varMax = 4     # Data maxs
    varInt = 0.5     # Data spacing
    levels = np.arange(varMin, varMax+0.0001, varInt)
    # list(np.arange(varMin,varMax,varInt))
    nlevs = len(levels)  # -- number of levels
    # -- convert list of floats to list of strings
    labels = ['{:.1f}'.format(x) for x in levels]
    labels = str(labels)
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
    res.lbLabelBarOn = False
    res.lbOrientation = "Horizontal"   # Draw labelbar horizontally.
    res.lbLabelFont = "complex_roman"
    res.lbLabelStrings = labels

    plotData = figData[2]  # 144 72 6
     # plotData = np.transpose(np.squeeze(np.array(plotData)[:, :, 1]))

     #
     # Loop 2 times and create 2 dummy plots;
     #
    nplots = plotType
    if plotType==4:
            f_row = 2    
    elif plotType==6:        
            f_row = 3

    f_col = 2 
    plots = []
    for n in range(nplots):
        print("plot #{}".format(n))
        # 子图标题
        res.tiMainString = figData[3]['subTitle'][n]  
        res.tiMainJust = "TopCenter"
        res.tiMainSide = "Top"
        res.tiMainFontHeightF = 0.0225
        res.tiMainFont = "Helvetica"
        
        # data
        subplotData = np.transpose(
            np.squeeze(np.array(plotData)[:, :, n]))
        # -- define _FillValue and missing_value if not existing
        subplotData[np.isnan(subplotData)] = -99999
        plots.append(Ngl.contour_map(wks, subplotData, res))

    # Resources for panelling
    pres = Ngl.Resources()
    pres.nglFrame = False
    pres.nglPanelLabelBar = True

    # Calculate start Y position for first row of plots
    height = 0.275            # we know this will be height of small plots
    extra = 1.0-(f_row*height)
    top = 1.0-(extra/f_row)  # 子图第一行所在的高度

    # Draw a title before we draw plots
    title = figData[3]['headLineTxt']
    txres = Ngl.Resources()
    txres.txJust = "BottomCenter"
    txres.txFontHeightF = 0.02
    Ngl.text_ndc(wks, title, 0.5, top+0.01, txres)

    # # Loop across plots and panel them on one page
    # for n in range(0,1):
    #     # Define location in a unit square for each set of plots.
    # add 子图序号
    pres.tiMainFont = "Helvetica-Bold"
    pres.nglPanelFigureStringsPerimOn     = False
    pres.nglPanelFigureStringsJust = "TopLeft"
    # pres.nglPanelFigureStringsFont = "Helvetica-Bold"
    pres.nglPanelFigureStringsOrthogonalPosF = -0.66
    pres.nglPanelFigureStringsParallelPosF= -0.525
    pres.nglPanelFigureStringsFontHeightF = 0.017
    pres.nglPanelFigureStrings = figData[3]['seqNum']
    
    pres.nglPanelTop = top-(0*height) # leave room for title
    pres.nglPanelBottom = top-(f_row*height)
    pres.nglPanelYWhiteSpacePercent = 3

    Ngl.panel(wks, plots[0:nplots+1], [f_row, f_col], pres)
    # plot = Ngl.contour_map(wks, plotData, res)

    # Ngl.draw(plot)

    Ngl.frame(wks)
    Ngl.destroy(wks)
    # Ngl.end()
