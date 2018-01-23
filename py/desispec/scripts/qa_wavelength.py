# Script for checking wavelength solutions
from __future__ import absolute_import, division

from desiutil.log import get_logger
import argparse

from desispec.qa import __offline_qa_version__

def parse(options=None):
    parser = argparse.ArgumentParser(description="Generate Wavelength QA [v{:s}]".format(__offline_qa_version__))
    parser.add_argument('--expid', type = int, required=True, help='Exposure ID')
    #parser.add_argument('--channels', type=str, help="List of channels to include. Default = b,r,z]")
    parser.add_argument('--reduxdir', type = str, default = None, metavar = 'PATH',
                        help = 'Override default path ($DESI_SPECTRO_REDUX/$SPECPROD) to processed data.')
    parser.add_argument('--sky', default = False, action='store_true',
                        help = 'Calculate shift in wavelengths from archived sky?')


    args = None
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args) :

    import numpy as np
    from astropy.table import Table

    from desispec.io import meta
    from desispec.flexure import flex_shift
    #from desispec.qa.qa_plots import exposure_sky_wave
    from desispec.io import get_files, find_exposure_night, read_frame
    from desispec.spectrum import Spectrum
    log=get_logger()

    log.info("starting")
    if args.reduxdir is None:
        specprod_dir = meta.specprod_root()
    else:
        specprod_dir = args.reduxdir
    #if args.channels is None:
    #    channels = ['b','r','z']
    #else:
    #    channels = [iarg for iarg in args.channels.split(',')]

    # Fiber QA
    if args.sky:
        # Loop on frames to generate the QA
        sky_QA = {}
        night = find_exposure_night(args.expid)
        frames_dict = get_files(filetype = str('frame'), night=night,
                                expid=args.expid)#, specprod_dir = self.specprod_dir)
        for camera,frame_fil in frames_dict.items():
            print(camera)
            sky_QA[camera] = {}
            # Read
            frame = read_frame(frame_fil)
            fibermap = frame.fibermap
            # Sky fibers
            skyfibers = np.where(frame.fibermap['OBJTYPE'] == 'SKY')[0]
            # Save X,Y
            sky_QA[camera]['x'] = fibermap['X_TARGET'][skyfibers]
            sky_QA[camera]['y'] = fibermap['Y_TARGET'][skyfibers]
            sky_QA[camera]['shifts'] = []
            # Loop to get offsets
            for skyfiber in skyfibers:
                spec = Spectrum(frame.wave, frame.flux[skyfiber,:], frame.ivar[skyfiber,:])
                flex_dict = flex_shift(camera[0], spec)#, debug=True)
                sky_QA[camera]['shifts'].append(flex_dict['shift'])
        # Generate a Table and print to screen
        tbl = Table()
        cameras = []
        shifts = []
        for camera in sky_QA.keys():
            cameras += [camera]*len(sky_QA[camera]['shifts'])
            shifts += sky_QA[camera]['shifts']
        tbl['camera'] = cameras
        tbl['shifts'] = shifts
        tbl.more()
        #
        import pdb; pdb.set_trace()



