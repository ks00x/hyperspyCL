import h5py as h5
import hyperspy.api as hs
import numpy as np

# switch of the too many warnings...
import warnings
warnings.filterwarnings('ignore')

def odemis_to_hyperspy(filename='sampledata/cltest.h5',specbin=1) :
    """
    odemis_to_hyperspy(filename='cltest.h5',specbin=1)
    open h5 files from Delmic Odemis and convert them to a hyperspy object
    reads in the spectral data (only) from hyperspectral maps and panchromatic
    CL maps using the PM
    1.arg: filename
    2.arg: binning for hyperspectral maps

    """

    f=h5.File(filename,'r')
    shome = 'Acquisition2//ImageData/'
    x = f[shome + 'Image']
    cdesc =f['Acquisition2/PhysicalData/ChannelDescription'].value[0].decode('utf-8')
    #print(cdesc)

    cltype = None
    if 'Spectrum' in cdesc :
        cltype = 'spectrum'
    elif 'CL intensity' in cdesc:
        cltype = 'panchrom'

    print('<' + filename + '> original shape :' ,x.shape, cltype)

    # strip unused dimensions and transpose/ reverse index order
    if cltype == 'panchrom' :
        xx=x[0,0,0,:,:].transpose((1,0))
        # just an image..
    else :
        xx=x[:,0,0,:,:].transpose((2,1,0))

    if cltype == 'spectrum' :
        #interpolate data to linearize the wavelength scale
        w  = f[shome + 'DimensionScaleC'].value *1e9
        wx = np.linspace(w.min(),w.max(),w.size)
        for i in np.arange(xx.shape[0]) :
            for k in np.arange(xx.shape[1]) :
                xx[i,k,:] = np.interp(wx,w,xx[i,k,:])

        wslope = wx[1]-wx[0]
        woffset = wx.min()
        #wx = np.arange(w.size)
        #wslope,woffset=np.polyfit(wx,w,1)
        s = hs.signals.Signal1D(xx)

    elif cltype == 'panchrom' :
        s = hs.signals.Signal2D(xx)
    else :
        print('unknown type')

    print('hyperspy shape :' ,s.data.shape)


    s.metadata.General.title = 'Odemis: ' + cdesc
    s.metadata.General.original_filename = filename
    s.metadata.General.notes = cltype
    s.axes_manager[0].name = 'pos x'
    s.axes_manager[0].scale = f[shome + 'DimensionScaleX'].value * 1e6
    s.axes_manager[0].offset = f[shome + 'XOffset'].value * 1e6
    s.axes_manager[0].units = 'um'


    s.axes_manager[1].name = 'pos y'
    s.axes_manager[1].scale = f[shome + 'DimensionScaleX'].value * 1e6
    s.axes_manager[1].offset = f[shome + 'YOffset'].value * 1e6
    s.axes_manager[1].units = 'um'

    if cltype == 'spectrum' :
        s.axes_manager[2].name = 'wavelength'
        s.axes_manager[2].units = 'nm'
        s.axes_manager[2].offset = woffset
        s.axes_manager[2].scale = wslope
        s.metadata.signal_type = 'CL'

    f.close()
    if (specbin > 1) and (cltype == 'spectrum'):
        return( s.rebin(scale=[1,1,specbin]) )
    else :
        return( s )
    #end odemis_to_hyperspy
    #######################



#odemis_to_hyperspy()


def odemisSEM_to_hyperspy(filename='sampledata/cltest.h5') :
    """
    odemis_to_hyperspy(filename='cltest.h5',specbin=1)
    open h5 files from Delmic Odemis and convert them to a hyperspy object
    reads in the spectral data (only) from hyperspectral maps and panchromatic
    CL maps using the PM
    1.arg: filename
    2.arg: binning for hyperspectral maps

    """

    f=h5.File(filename,'r')
    shome = 'Acquisition1//ImageData/'
    x = f[shome + 'Image']
    cdesc =f['Acquisition1/PhysicalData/ChannelDescription'].value[0].decode('utf-8')
    #print(cdesc)


    print('<' + filename + '> original shape :' ,x.shape)
    # strip unused dimensions and transpose/ reverse index order
    xx=x[0,0,0,:,:].transpose((1,0))

    s = hs.signals.Signal2D(xx)
    print('hyperspy shape :' ,s.data.shape)


    s.metadata.General.title = 'Odemis: ' + cdesc
    s.metadata.General.original_filename = filename
    #s.metadata.General.notes = cltype
    s.axes_manager[0].name = 'pos x'
    s.axes_manager[0].scale = f[shome + 'DimensionScaleX'].value * 1e6
    s.axes_manager[0].offset = f[shome + 'XOffset'].value * 1e6
    s.axes_manager[0].units = 'um'


    s.axes_manager[1].name = 'pos y'
    s.axes_manager[1].scale = f[shome + 'DimensionScaleX'].value * 1e6
    s.axes_manager[1].offset = f[shome + 'YOffset'].value * 1e6
    s.axes_manager[1].units = 'um'


    f.close()
    return( s )
    #end odemisSEM_to_hyperspy
    #######################
