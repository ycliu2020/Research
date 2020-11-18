#
#  File:
#    newcolor4.py
#
#  Synopsis:
#    Illustrates new color capabilities in PyNGL 1.5.0.
#
#  Categories:
#    contour plots
#
#  Author:
#    Mary Haley, based on NCL example
#  
#  Date of initial publication:
#    November 2012
#
#  Description:
#    This example shows how to draw nine plots with three different colormaps
#
#  Effects illustrated:
#    o  Using the new "cnFillPalette" resource.
#    o  Panelling plots.
#    o  Drawing an Aitoff map
#    o  Generating dummy data
# 
#  Output:
#     A single visualization with nine plots
#
from __future__ import print_function
import numpy,os
import Ngl,Nio

#---Start the graphics section
wks_type = "png"
wks = Ngl.open_wks (wks_type,"newcolor4")

#---Generate some dummy lat/lon data
nlat      =  64
nlon      = 128
lat       = Ngl.fspan(-90,90,nlat)
lon       = Ngl.fspan(-178.5,178.5,nlon)

#---Values to use for contour labelbar
dmins = [-20.,-10.,-21.]     # Data mins
dmaxs = [ 16., 10., 17.]     # Data maxs
dspas = [  4.,  1.,  3.]     # Data spacing

# One color map per row of plots
colormaps = ["wgne15","StepSeq25","BlueDarkRed18"]

# Create resource list for customizing contour over maps
res                        = Ngl.Resources()

res.nglMaximize            = False
res.nglFrame               = False
res.nglDraw                = False

res.mpProjection           = "Aitoff"
res.mpOutlineOn            = True

res.mpPerimOn              = False
res.mpGridAndLimbOn        = False
res.pmTickMarkDisplayMode  = "Never"

res.cnFillOn               = True
res.cnLinesOn              = False
res.cnLineLabelsOn         = False
res.cnLevelSelectionMode   = "ManualLevels"
res.lbLabelBarOn           = False

res.sfXArray               = lon
res.sfYArray               = lat

#
# Loop 9 times and create 9 dummy plots; each group
# of 3 plots has the same color map.
#
nplots = 9
plots  = []
for n in range(nplots):
  print("plot #{}".format(n))
  dmin = dmins[n//3]
  dmax = dmaxs[n//3]
  dspa = dspas[n//3]

  mstart = numpy.random.uniform(10,25,1).astype(int)
  mend   = numpy.random.uniform(10,25,1).astype(int)
  xstart = numpy.random.uniform(dmin,dmin+2,1)
  xend   = numpy.random.uniform(dmax-2,dmax,1)

#---This is a new resource added in PyNGL 1.5.0
  res.cnFillPalette          = colormaps[n//3]

  res.cnMinLevelValF         = dmin
  res.cnMaxLevelValF         = dmax
  res.cnLevelSpacingF        = dspa

  data = Ngl.generate_2d_array([nlat,nlon],mstart[0],mend[0],xstart[0],xend[0])

  plots.append(Ngl.contour_map(wks,data,res))

# Resources for panelling
pres                  = Ngl.Resources() 
pres.nglFrame         = False
pres.nglPanelLabelBar = True

# Calculate start Y position for first row of plots
height = 0.15            # we know this will be height of small plots
extra  = 1.0-(3*height)
top    = 1.0-(extra/2.)

# Draw a title before we draw plots
title               = "Multiple panels on one page, 3 different colormaps"
txres               = Ngl.Resources()
txres.txJust        = "BottomCenter"
txres.txFontHeightF = 0.02
Ngl.text_ndc(wks,title,0.5,top+0.01,txres)

# Loop across plots and panel them on one page
for n in range(0,3):
# Define location in a unit square for each set of plots.
  pres.nglPanelTop    = top-(n*height)
  pres.nglPanelBottom = top-((n+1)*height)

  Ngl.panel(wks,plots[n*3:n*3+3],[1,3],pres)
  
Ngl.frame(wks)
Ngl.end()


#
#  File:
#    ngl09p.py
#
#  Synopsis:
#    Illustrates animation of contours over a map and using masked
#    arrays.
#
#  Category:
#    Contours over maps
#
#  Author:
#    Fred Clare (based on a code of Mary Haley)
#  
#  Date of initial publication:
#    September, 2005
#
#  Description:
#    This example utilizes masked arrays handle missing data values
#
#  Effects illustrated:
#    o  masked arrays to handle missing values.
#    o  Adding a cyclic point.
#    o  Looping to produce an animation.
#    o  Using a stereographic map projection.
#    o  Many map resource settings via the resource file ngl09p.res.
#    o  Label bar resource settings via the resource file ngl09p.res.
# 
#  Output:
#    The example produces twelve contour plots over a map using
#    a stereographic map projection.
#
#  Notes:
#    This example requires the resource file ngl09p.res.
#     

from __future__ import print_function
import numpy
#
#  Import the masked array module.
#
from numpy import ma as MA
import os

#
#  Import Ngl support functions.
#
import Ngl

#
#  Import Nio for a NetCDF reader.
#
import Nio

#
#  To use the ScientificPython module to read in the netCDF file,
#  comment out the above "import" command, and uncomment the 
#  import line below.
#
# from Scientific.IO.NetCDF import NetCDFFile

#
#  Open the netCDF files, get variables.
#
data_dir = Ngl.pynglpath("data")
ice1     = Nio.open_file(os.path.join(data_dir,"cdf","fice.nc"),"r")

#
#  This is the ScientificPython method for opening a netCDF file.
#
# ice1     = NetCDFFile(data_dir + "/cdf/fice.nc","r")

#
#  Create a masked array to accommodate missing values in the fice variable.
#
fice = ice1.variables["fice"]  # fice[120,49,100]
ficea = fice[:,:,:]
fill_value = None
if (hasattr(fice,"missing_value")):
  fill_value = fice.missing_value
elif (hasattr(fice,"_FillValue")):
  fill_value = fice._FillValue
fice_masked = MA.transpose(MA.masked_values(ficea,fill_value),(1,2,0))

hlat = ice1.variables["hlat"]  # hlat[49]
hlon = ice1.variables["hlon"]  # hlon[100]


dimf     = fice.shape  # Define an array to hold long-term monthly means.
ntime    = fice.shape[0]
nhlat    = fice.shape[1]
nhlon    = fice.shape[2]

nmo    = 0
month  = nmo+1

icemon = MA.zeros((nhlat,nhlon),dtype=float)
for i in range(fice_masked.shape[0]):
  for j in range(fice_masked.shape[1]):
    icemon[i,j] = MA.average(fice_masked[i,j,0:ntime:12])

#
#  Fill the places where icemon is zero with the fill value.
#
icemon = MA.masked_values(icemon,0.,rtol=0.,atol=1.e-15)
icemon = MA.filled(icemon,fill_value)

                       # Calculate the January (nmo=0) average.


nsub = 16 # Subscript location of northernmost hlat to be plotted.

cmap = numpy.array([                                           \
         [1.00,1.00,0.50], [0.00,0.00,0.50], [0.50,1.00,1.00], \
         [0.50,0.00,0.00], [1.00,0.00,1.00], [0.00,1.00,1.00], \
         [1.00,1.00,0.00], [0.00,0.00,1.00], [0.00,1.00,0.00], \
         [1.00,0.00,0.00] ],dtype=float)

wks_type = "png"
wks = Ngl.open_wks(wks_type,"ngl09p")

resources = Ngl.Resources()

# Add a longitude cyclic point
icemonnew,hlonnew = Ngl.add_cyclic(icemon[0:nsub+1,:],hlon[:])

resources.sfMissingValueV = fill_value
resources.cnFillPalette   = cmap
resources.sfXArray        = hlonnew   # Necessary for overlay on a map.
resources.sfYArray        = hlat[0:nsub+1]
resources.tiMainString    = "CSM Y00-99 Mean Ice Fraction Month ={}".format(month)

resources.pmTickMarkDisplayMode = "Never"

map = Ngl.contour_map(wks,icemonnew,resources) # Draw a contour
                                               # over a map.

nmos = 12    # Specify the number of months in the loop (max 120).
for nmo in range(1,nmos): 
  month  = nmo+1
  for i in range(fice_masked.shape[0]):
    for j in range(fice_masked.shape[1]):
      icemon[i,j] = MA.average(fice_masked[i,j,nmo:ntime:12])
  icemon = MA.masked_values(icemon,0.,rtol=0.,atol=1.e-15)
  icemon = MA.filled(icemon,fill_value)

  resources.tiMainString = "CSM Y00-99 Mean Ice Fraction Month ={}".format(month)
  map = \
    Ngl.contour_map(wks,Ngl.add_cyclic(icemon[0:nsub+1,:]),resources)

del icemon       # Clean up.
del icemonnew 
del map 

Ngl.end()