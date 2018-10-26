import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.lines import Line2D
import scipy.ndimage.filters as filters
import scipy.interpolate as interpolate

class plotCLASS: 

  def __init__(self):
    self.NANVAL = 1.234567e-10
    self.colors = ['k', 'b', 'r', 'g', 'c', 'y', 'm', 'orangered', 'navy']

  def getShape(self, fileName):
    shape = []
    ind1 = fileName.find("[")
    ind2 = fileName.find(",", ind1)
    while ind2 is not -1:
      print("num", fileName[ind1+1:ind2])
      shape.append(int(fileName[ind1+1:ind2]))
      ind1 = ind2
      ind2 = fileName.find(",", ind1+1)
    ind2 = fileName.find("]", ind1)
    print("inds", ind1, ind2)
    shape.append(int(fileName[ind1+1: ind2]))

    return shape




  def importImage(self, fileName, getShape=True):
    print("file", fileName)
    image = np.fromfile(fileName, dtype = np.double)
    if getShape:
      shape = self.getShape(fileName)
    else:
      shape = None
    if (shape is not None) and (len(shape) > 1):
      image = np.reshape(image, shape)

    return image, shape


  def printHist(self, fileName, Nbins, outputName, 
      options=None, returnPlt=False):

    fig, ax = plt.subplots()
    handles = []

    vals,_ = self.importImage(fileName, False)
    print("size", vals.shape)
    if vals.shape[0] is 0:
      return
    ax.hist(vals, Nbins, color="k")

    if options is not None:
      fig, ax = self.beautify(fig, ax, options, handles)
    if returnPlt:
      return fig,ax
    else:
      fig.savefig(outputName + ".png")
      plt.close()



  def print1d(self, 
      inpImages, outputName, 
      xRange=None, normalize=None, 
      scale=None, isFile=True,
      options=None):

    Nimages = 0
    if isFile:
      Nimages = len(inpImages)
    else:
      if len(inpImages.shape) is 1:
        Nimages = 1
        inpImages = np.reshape(inpImages, (1,-1))
      else:
        Nimages = inpImages.shape[0]


    handles = []
    fig, ax = plt.subplots()
    for i in range(Nimages):
      if isFile:
        image,_ = self.importImage(inpImages[i], False)
      else:
        image = inpImages[i,:]
      X = np.arange(image.shape[0])
      if xRange is not None:
        X = xRange[0] + X*(xRange[1] - xRange[0])/float(image.shape[0])
      if normalize is not None:
        if normalize is "max":
          if "normalize" in options:
            image /= np.amax(image[options["normalize"][0]:options["normalize"][1]])
          else:
            print("in max", np.amax(image))
            image /= np.amax(image)
        elif normalize is "abs":
          image = np.abs(image)
        elif normalize is "0min":
          image -= np.amin(image)
        else:
          image -= np.mean(image[7:25])
          image /= max(image[7:20])
      if scale is not None:
        image *= scale[i]

      h, = ax.plot(X, image, color=self.colors[i], linestyle='-')
      handles.append(h)

    ax.set_xlim(X[0], X[-1])

    if options is not None:
      fig, ax = self.beautify(fig, ax, options, handles)
    fig.savefig(outputName + ".png")
    plt.close()




  def printLineOut(self, fileName, axis, inds, outputName, 
      X=None, xRange=None, samePlot=True, options=None):

    image,shape = self.importImage(fileName)
    print(image.shape)

    if X is None:
      if axis is 0:
        xLen = shape[1]
      else:
        xLen = shape[0]

      X = np.arange(xLen)
      if xRange is not None:
        X = xRange[0] + X*(xRange[1] - xRange[0])/float(xLen)
 
    handles = []
    fig, ax = plt.subplots()
    for i,ind in enumerate(inds):

      if not samePlot:
        plt.close()
        fig, ax = plt.subplots()

      print("ind",ind)
      if axis is 0:
        inp = np.reshape(image[ind,:], (-1))
      elif axis is 1:
        inp = np.reshape(image[:,ind], (-1))
      else:
        print("ERROR: Does not support axis = {}".format(axis))
        sys.exit(0)

      if options is not None:
        if "Qscale" in options:
          inp *= options["Qscale"]
        if "smooth" in options:
          inp = filters.gaussian_filter1d(inp, options["smooth"])

      if samePlot:
        print("color", self.colors[i])
        print("shapes", X.shape, inp.shape)
        h, = ax.plot(X, inp, color=self.colors[i], linestyle='-')
        handles.append(h)
      else:
        h, = ax.plot(X, inp, linestyle='-')
        ax.set_xlim(X[0], X[-1])
        if options is not None:
          fig, ax = self.beautify(fig, ax, options, handles)
        fig.savefig(outputName + "_" + str(ind) + ".png")

    if samePlot:
      ax.set_xlim(X[0], X[-1])

      if options is not None:
        fig, ax = self.beautify(fig, ax, options, handles)
      fig.savefig(outputName + ".png")
      plt.close()




  #subplot(aspect='equal')
  def print2d(self, 
      inpImage, outputName, 
      X=None, xRange=None, 
      Y=None, yRange=None, 
      xRebin=None, yRebin=None, 
      scale=1, isFile=True,
      options=None):

    if isFile:
      imageO,shape = self.importImage(inpImage)
      if "Diff" in inpImage:
        image = np.zeros((imageO.shape[0], 555/5), dtype=float)
        for i in range(555/5):
          image[:,i] = np.mean(imageO[:,i*5:(i+1)*5], axis=1)
      else :
        image = imageO
    else:
      image = inpImage

    shape = np.array(image.shape)
    image *= scale

    
    print("shape",shape)

    if X is None:
      X = np.arange(image.shape[0])
      if xRange is not None:
        X = xRange[0] + X*(xRange[1] - xRange[0])/float(shape[0])
    if Y is None:
      Y = np.arange(image.shape[1])
      if yRange is not None:
        Y = yRange[0] + Y*(yRange[1] - yRange[0])/float(shape[1])

    X,Y = np.meshgrid(X,Y)
    print("meshgrid")

    fig, ax = plt.subplots()

    cNorm = None
    if options is not None:
      if "colorSTDrange" in options:
        mean = np.mean(image[:,int(0.2*shape[1]):int(0.7*shape[1])])
        std = np.std(image[:,int(0.2*shape[1]):int(0.7*shape[1])])
        if mean > 0:
          vRange = mean + options["colorSTDrange"]*std
        else:
          vRange = options["colorSTDrange"]*std - mean
        vMin = -1*vRange
        vMax = vRange
      elif "colorRange" in options:
        vMin = options["colorRange"][0]
        vMax = options["colorRange"][1]
      else:
        vMin = np.amin(image)
        vMax = np.amax(image)

      if "colorMap" in options:
        cMap = options["colorMap"]
      else:
        cMap = 'jet'

      if "colorNorm" in options:
        if options["colorNorm"] is "log":
          print("VMM ",vMin, vMax)
          cNorm = colors.LogNorm(vMin, vMax)

      if "interpolate" in options:
        print("shapes", X.shape, Y.shape, image.shape)
        print(X[0,0], X[0,-1])
        tck = interpolate.bisplrep(X[:,:-1], Y[:,:-1], image.transpose(), s=0.01)
        X = np.linspace(X[0,0], X[0,-1], options["interpolate"][0])
        Y = np.linspace(Y[0,0], Y[-1,0], options["interpolate"][1])
        X,Y = np.meshgrid(X,Y)
        print("shapes2 ",X.shape, Y.shape)
        print(X[0,:], Y[:,-1])
        image = interpolate.bisplev(X[0,:], Y[:,0], tck)

    else:
      vMin = np.amin(image)
      vMax = np.amax(image)
      cMap = 'jet'

    print("plotting", X.shape, Y.shape, image.shape)
    plot = ax.pcolor(X, Y, image.transpose(), 
              norm=cNorm,
              cmap=cMap, 
              vmin=vMin, vmax=vMax)

    print("plotted")
    ax.set_xlim(X[0,0], X[0,-1])
    ax.set_ylim(Y[0,0], Y[-1,0])
    cbar = fig.colorbar(plot)

    if options is not None:
      fig, ax = self.beautify(fig, ax, options)

    plt.tight_layout()
    fig.savefig(outputName + ".png")
    plt.close()




  def beautify(self, fig, ax, options, handles=None):
    if "yLim" in options:
      ax.set_ylim(options["yLim"])
    if "xLim" in options:
      ax.set_xlim(options["xLim"])
    if "yTitle" in options:
      ax.set_ylabel(options["yTitle"])
    if "xTitle" in options:
      ax.set_xlabel(options["xTitle"])
    if "xSlice" in options:
      ax.set_xlim(options["xSlice"])
    if "ySlice" in options:
      ax.set_ylim(options["ySlice"])
    if "xLog" in options:
      if options["xLog"]:
        ax.set_xscale("log", nonposx='clip')
    if "yLog" in options:
      if options["yLog"]:
        ax.set_yscale("log", nonposy='clip')
    if "labels" in options:
      if handles is None:
        print("ERROR: plotting the legend requirese handles!!!")
        sys.exit(0)
      if "legOpts" in options:
        ax.legend(tuple(handles), tuple(options["labels"]), **options["legOpts"])
      else:
        ax.legend(tuple(handles), tuple(options["labels"]))
    if "line" in options:
      l = Line2D(options["line"][0], options["line"][1], 
          color=options["line"][2], linewidth=options["line"][3])
      ax.add_line(l)

    return fig, ax


