import numpy as np
import json
import yaml
from desispec.io import findfile
import os,sys
from desispec.quicklook import qlexceptions,qllogger

class Config(object):
    """ 
    A class to generate Quicklook configurations for a given desi exposure. 
    expand_config will expand out to full format as needed by quicklook.setup
    """

    def __init__(self, configfile, night, camera, expid, singqa, amps=True,rawdata_dir=None,specprod_dir=None, outdir=None,qlf=False):
        """
        configfile: a configuration file for QL eg: desispec/data/quicklook/qlconfig_dark.yaml
        night: night for the data to process, eg.'20191015'
        camera: which camera to process eg 'r0'
        expid: exposure id for the image to be processed 
        amps: for outputing amps level QA
        Note:
        rawdata_dir and specprod_dir: if not None, overrides the standard DESI convention       
        """
  
        #- load the config file and extract
        self.conf = yaml.load(open(configfile,"r"))
        self.night = night
        self.expid = expid
        self.camera = camera
        self.singqa = singqa
        self.amps = amps
        self.rawdata_dir = rawdata_dir 
        self.specprod_dir = specprod_dir
        self.outdir = outdir
        self.dumpintermediates = self.conf["WriteIntermediatefiles"]
        self.writepreprocfile = self.conf["WritePreprocfile"]
        self.writeskymodelfile = self.conf["WriteSkyModelfile"]
        self.writestaticplots = self.conf["WriteStaticPlots"]
        self.usesigma = self.conf["UseResolution"]
        try:
            self.flexure = self.conf["Flexure"]    
        except:
            self.flexure = False
        self.pipeline = self.conf["Pipeline"]
        if not self.flexure and "Flexure" in self.pipeline:
            self.pipeline.remove("Flexure")
        self.algorithms = self.conf["Algorithms"]
        self._palist = Palist(self.pipeline,self.algorithms)
        self.pamodule = self._palist.pamodule
        self.qamodule = self._palist.qamodule
        if "BoxcarExtract" in self.algorithms.keys():
            if "wavelength" in self.algorithms["BoxcarExtract"].keys():
                self.wavelength = self.algorithms["BoxcarExtract"]["wavelength"][self.camera[0]]
        else: self.wavelength = None
        self._qlf=qlf
        qlog=qllogger.QLLogger(name="QLConfig")
        self.log=qlog.getlog()
        self._qaRefKeys={"Check_HDUs":"HDUs_OK","Trace_Shifts":"TRACE_REF","Bias_From_Overscan":"BIAS_AMP", "Get_RMS":"NOISE_AMP", "Count_Pixels":"LITFRAC_AMP", "Calc_XWSigma":"XWSIGMA", "CountSpectralBins":"NGOODFIB", "Sky_Peaks":"PEAKCOUNT", "Sky_Continuum":"SKYCONT", "Integrate_Spec":"DELTAMAG_TGT", "Sky_Residual":"MED_RESID", "Calculate_SNR":"FIDSNR_TGT"}

    @property
    def mode(self):
        """ what mode of QL, online? offline?
        """
        return self._palist.mode

    @property
    def qlf(self):
        return self._qlf

    @property
    def palist(self): 
        """ palist for this config
            see :class: `Palist` for details.
        """
        return self._palist.palist

    @property
    def qalist(self):
        """ qalist for the given palist
        """
        return self._palist.qalist

    @property
    def paargs(self,psfspfile=None):
        """
        Many arguments for the PAs are taken default. Some of these may need to be variable
        psfspfile is for offline extraction case
        """
        wavelength=self.wavelength
        if self.wavelength is None:
            #- setting default wavelength for extraction for different cam
            if self.camera[0] == 'b':
                self.wavelength='3570,5730,0.8'
            elif self.camera[0] == 'r':
                self.wavelength='5630,7740,0.8'
            elif self.camera[0] == 'z':
                self.wavelength='7420,9830,0.8'

        #- Make kwargs less verbose using '%%' marker for global variables. Pipeline will map them back
        paopt_initialize={'camera': self.camera}

        if self.writepreprocfile:
            preprocfile=self.dump_pa("Preproc")
        else: 
            preprocfile = None
        paopt_preproc={'camera': self.camera,'dumpfile': preprocfile}

        if self.dumpintermediates:
            framefile=self.dump_pa("BoxcarExtract")
            fframefile=self.dump_pa("ApplyFiberFlat_QL")
            sframefile=self.dump_pa("SkySub_QL")
        else:
            framefile=None
            fframefile=None
            sframefile=None

        if self.conf["Flavor"] == 'arcs':
            arcimg=findfile('preproc',night=self.night,expid=self.expid,camera=self.camera,specprod_dir=self.specprod_dir)
            flatimg=findfile('preproc',night=self.night,expid=self.conf["FiberflatExpid"],camera=self.camera,specprod_dir=self.specprod_dir)
            bootfile=findfile('psfboot',expid=self.expid,night=self.night,camera=self.camera,specprod_dir=self.specprod_dir)
            psffile=findfile('psf',expid=self.expid,night=self.night,camera=self.camera,specprod_dir=self.specprod_dir)
        else:
            arcimg=None
            flatimg=None
            bootfile=None
            psffile=None

        if self.flexure:
            preproc_file=findfile('preproc',self.night,self.expid,self.camera,specprod_dir=self.specprod_dir)
            inputpsf=findfile(self.conf["PSFType"],self.night,self.psfexpid,self.camera,specprod_dir=self.specprod_dir)
            outputpsf=findfile('psf',self.night,self.expid,self.camera,specprod_dir=self.specprod_dir)
        else:
            preproc_file=None
            inputpsf=None
            outputpsf=None

        paopt_flexure={'preprocFile':preproc_file, 'inputPSFFile': inputpsf, 'outputPSFFile': outputpsf}

        paopt_bootcalib={'ArcLampImage':arcimg, 'FlatImage':flatimg, 'outputFile':bootfile}

        paopt_extract={'BoxWidth': 2.5, 'FiberMap': self.fibermap, 'Wavelength': self.wavelength, 'Nspec': 500, 'PSFFile': self.psf,'usesigma': self.usesigma, 'dumpfile': framefile}

        paopt_resfit={'PSFbootfile':bootfile, 'PSFoutfile': psffile, 'usesigma': self.usesigma}

        if self.conf["Flavor"] == 'flat':
            fiberflatfile=findfile('fiberflat',night=self.night,expid=self.conf["FiberflatExpid"],camera=self.camera,specprod_dir=self.specprod_dir)
        else:
            fiberflatfile=None

        paopt_comflat={'outputFile': fiberflatfile}

        paopt_apfflat={'FiberFlatFile': self.fiberflat, 'dumpfile': fframefile}

        if self.writeskymodelfile:
            outskyfile = findfile('sky',night=self.night,expid=self.expid, camera=self.camera, rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir,outdir=self.outdir)
        else:
            outskyfile=None       
        paopt_skysub={'Outskyfile': outskyfile, 'dumpfile': sframefile, 'Apply_resolution': self.usesigma}

        paopts={}
        defList={
            'Initialize':paopt_initialize,
            'Preproc':paopt_preproc,
            'Flexure':paopt_flexure,
            'BootCalibration':paopt_bootcalib,
            'BoxcarExtract':paopt_extract,
            'ResolutionFit':paopt_resfit,
            'ComputeFiberflat_QL':paopt_comflat,
            'ApplyFiberFlat_QL':paopt_apfflat,
            'SkySub_QL':paopt_skysub
        }

        def getPAConfigFromFile(PA,algs):
            def mergeDicts(source,dest):
                for k in source:
                    if k not in dest:
                        dest[k]=source[k]
            userconfig={}
            if PA in algs:
                fc=algs[PA]
                for k in fc: #do a deep copy leave QA config out
                    if k != "QA":
                        userconfig[k]=fc[k]
            defconfig={}
            if PA in defList:
                defconfig=defList[PA]
            mergeDicts(defconfig,userconfig)
            return userconfig

        for PA in self.palist:
            paopts[PA]=getPAConfigFromFile(PA,self.algorithms)
        #- Ignore intermediate dumping and write explicitly the outputfile for 
        self.outputfile=self.dump_pa(self.palist[-1]) 

        return paopts 
        
    def dump_pa(self,paname):
        """
        dump the PA outputs to respective files. This has to be updated for fframe and sframe files as QL anticipates for dumpintermediate case.
        """
        pafilemap={'Preproc': 'preproc', 'Flexure': None, 'BootCalibration': 'psfboot', 'BoxcarExtract': 'frame', 'ResolutionFit': None, 'ComputeFiberflat_QL': 'fiberflat', 'ApplyFiberFlat_QL': 'fframe', 'SkySub_QL': 'sframe'}
        
        if paname in pafilemap:
            filetype=pafilemap[paname]
        else:
            raise IOError("PA name does not match any file type. Check PA name in config") 

        pafile=None
        if filetype is not None:
            pafile=findfile(filetype,night=self.night,expid=self.expid,camera=self.camera,rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir,outdir=self.outdir)

        return pafile

    def dump_qa(self): 
        """ 
        yaml outputfile for the set of qas for a given pa
        Name and default locations of files are handled by desispec.io.meta.findfile
        """

        #- both PA level and QA level outputs
        qa_pa_outfile = {}
        qa_pa_outfig = {}
        qa_outfile = {}
        qa_outfig = {}
        for PA in self.palist:
            #- pa level outputs
            if self.dumpintermediates:
                qa_pa_outfile[PA] = self.io_qa_pa(PA)[0]
                qa_pa_outfig[PA] = self.io_qa_pa(PA)[1]
            else:
                qa_pa_outfile[PA] = None
                qa_pa_outfig[PA] = None
            #- qa_level output
            for QA in self.qalist[PA]:
                qa_outfile[QA] = self.io_qa(QA)[0]
                qa_outfig[QA] = self.io_qa(QA)[1]
                
                #- make path if needed
                path = os.path.normpath(os.path.dirname(qa_outfile[QA]))
                if not os.path.exists(path):
                    os.makedirs(path)

        return ((qa_outfile,qa_outfig),(qa_pa_outfile,qa_pa_outfig))

    @property
    def qaargs(self):

        qaopts = {}
        referencemetrics=[]        

        for PA in self.palist:
            for qa in self.qalist[PA]: #- individual QA for that PA
                if self.writestaticplots:
                    qaplot = self.dump_qa()[0][1][qa]
                else:
                    qaplot = None

                pa_yaml = PA.upper()
                params=self._qaparams(qa)
                qaopts[qa]={'night' : self.night, 'expid' : self.expid,
                            'camera': self.camera, 'paname': PA, 'PSFFile': self.psf,
                            'amps': self.amps, 'qafile': self.dump_qa()[0][0][qa],
                            'qafig': qaplot, 'FiberMap': self.fibermap,
                            'param': params, 'qlf': self.qlf, 'refKey':self._qaRefKeys[qa],
                            'singleqa' : self.singqa}
                if self.singqa is not None:
                    qaopts[qa]['rawdir']=self.rawdata_dir
                    qaopts[qa]['specdir']=self.specprod_dir
                    if qa == 'Sky_Residual':
                        skyfile = findfile('sky',night=self.night,expid=self.expid, camera=self.camera, rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir,outdir=self.outdir)
                        qaopts[qa]['SkyFile']=skyfile

                if self.reference != None:
                    refkey=qaopts[qa]['refKey']
                    for padict in range(len(self.reference)):
                        pa_metrics=self.reference[padict].keys()
                        if refkey in pa_metrics:
                            qaopts[qa]['ReferenceMetrics']={'{}'.format(refkey): self.reference[padict][refkey]}
        return qaopts

    def _qaparams(self,qa):
            
        params={}
        if self.algorithms is not None:
            for PA in self.palist:
                if qa in self.qalist[PA]:
                    params[qa]=self.algorithms[PA]['QA'][qa]['PARAMS']

        else:
            if qa == 'Count_Pixels':
                params[qa]= dict(
                                CUTLO = 100,
                                 CUTHI = 500
                                )
            elif qa == 'CountSpectralBins':
                params[qa]= dict(
                                 CUTLO = 100,   # low threshold for number of counts
                                 CUTMED = 250,
                                 CUTHI = 500
                                )
            elif qa == 'Sky_Residual':
                params[qa]= dict(
                                 PCHI_RESID=0.05, # P(Chi^2) limit for bad skyfiber model residuals
                                 PER_RESID=95.,   # Percentile for residual distribution
                                 BIN_SZ=0.1,) # Bin size for residual/sigma histogram
            else:
                params[qa]= dict()
        
        return params[qa]

    def io_qa_pa(self,paname):
        """
        Specify the filenames: json and png of the pa level qa files"
        """
        filemap={'Initialize': 'initial',
                 'Preproc': 'preproc',
                 'Flexure': 'flexure',
                 'BootCalibration': 'bootcalib',
                 'BoxcarExtract': 'boxextract',
                 'ResolutionFit': 'resfit',
                 'ComputeFiberflat_QL': 'computeflat',
                 'ApplyFiberFlat_QL': 'fiberflat',
                 'SkySub_QL': 'skysub'
                 }

        if paname in filemap:
            outfile=findfile('ql_file',night=self.night,expid=self.expid, camera=self.camera, rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir,outdir=self.outdir)
            outfile=outfile.replace('qlfile',filemap[paname])
            outfig=findfile('ql_fig',night=self.night,expid=self.expid, camera=self.camera, rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir,outdir=self.outdir)
            outfig=outfig.replace('qlfig',filemap[paname])
        else:
            raise IOError("PA name does not match any file type. Check PA name in config for {}".format(paname))

        return (outfile,outfig)


    def io_qa(self,qaname):
        """
        Specify the filenames: json and png for the given qa output
        """
        filemap={'Check_HDUs':'checkHDUs',
                 'Trace_Shifts':'trace',
                 'Bias_From_Overscan': 'getbias',
                 'Get_RMS' : 'getrms',
                 'Count_Pixels': 'countpix',
                 'Calc_XWSigma': 'xwsigma',
                 'CountSpectralBins': 'countbins',
                 'Sky_Continuum': 'skycont',
                 'Sky_Peaks': 'skypeak',
                 'Sky_Residual': 'skyresid',
                 'Integrate_Spec': 'integ',
                 'Calculate_SNR': 'snr'
                 }

        if qaname in filemap:
            outfile=findfile('ql_file',night=self.night,expid=self.expid, camera=self.camera, rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir,outdir=self.outdir)
            outfile=outfile.replace('qlfile',filemap[qaname])
            outfig=findfile('ql_fig',night=self.night,expid=self.expid, camera=self.camera, rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir,outdir=self.outdir)
            outfig=outfig.replace('qlfig',filemap[qaname])
        else:
            raise IOError("QA name does not match any file type. Check QA name in config for {}".format(qaname))

        return (outfile,outfig)

    def expand_config(self):
        """
        config: desispec.quicklook.qlconfig.Config object
        """

        self.log.debug("Building Full Configuration")

        self.program = self.conf["Program"]
        self.flavor = self.conf["Flavor"]
        self.debuglevel = self.conf["Debuglevel"]
        self.period = self.conf["Period"]
        self.timeout = self.conf["Timeout"]
        self.fiberflatexpid = self.conf["FiberflatExpid"]
        self.psfexpid = self.conf["PSFExpid"]
        self.psftype = self.conf["PSFType"]
        self.templateexpid = self.conf["TemplateExpid"]
        self.templatenight = self.conf["TemplateNight"]

        #- some global variables:
        self.rawfile=findfile("raw",night=self.night,expid=self.expid,camera=self.camera,rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir)

        self.fibermap=findfile("fibermap", night=self.night,expid=self.expid,camera=self.camera,rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir)
 
        self.fiberflat=findfile("fiberflat",night=self.night,expid=self.fiberflatexpid,camera=self.camera,rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir) #- TODO: Assuming same night for calibration files (here and psf)
        
        self.psf=findfile(self.psftype,night=self.night,expid=self.psfexpid,camera=self.camera,rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir)

        #- Get reference metrics from template json file
        template=findfile('ql_mergedQA_file',night=self.templatenight,expid=self.templateexpid,camera=self.camera,rawdata_dir=self.rawdata_dir,specprod_dir=self.specprod_dir)
        self.reference=None
        if os.path.isfile(template):
            try:
                with open(template) as reference:
                    refdict=json.load(reference)
                    tasks=refdict["TASKS"]
                    tasklist=[]
                    for task in tasks.keys():
                        tasklist.append(task)
                    ref_met=[]
                    for ttask in tasklist:
                        ref_met.append(tasks[ttask]['METRICS'])
                    self.reference=ref_met
            except:
                self.log.warning("WARNING, template file is malformed %s"%template)
        else:
            self.log.warning("WARNING, can't open template file %s"%template)

        outconfig={}

        outconfig['Night'] = self.night
        outconfig['Program'] = self.program
        outconfig['Flavor'] = self.flavor
        outconfig['Camera'] = self.camera
        outconfig['Expid'] = self.expid
        outconfig['DumpIntermediates'] = self.dumpintermediates
        outconfig['FiberMap']=self.fibermap
        outconfig['FiberFlatFile'] = self.fiberflat
        outconfig['PSFFile'] = self.psf
        outconfig['Period'] = self.period

        pipeline = []
        for ii,PA in enumerate(self.palist):
            pipe={'OutputFile': self.dump_qa()[1][0][PA]} #- integrated QAs for that PA. 
            pipe['PA'] = {'ClassName': PA, 'ModuleName': self.pamodule, 'kwargs': self.paargs[PA]}
            pipe['QAs']=[]
            for jj, QA in enumerate(self.qalist[PA]):
                pipe_qa={'ClassName': QA, 'ModuleName': self.qamodule, 'kwargs': self.qaargs[QA]}
                pipe['QAs'].append(pipe_qa)
            pipe['StepName']=PA
            pipeline.append(pipe)

        outconfig['PipeLine'] = pipeline
        outconfig['RawImage'] = self.rawfile
        outconfig['OutputFile'] = self.outputfile
        outconfig['singleqa'] = self.singqa
        outconfig['Timeout'] = self.timeout

        #- Check if all the files exist for this QL configuraion
        check_config(outconfig,self.singqa)

        return outconfig

def check_config(outconfig,singqa):
    """
    Given the expanded config, check for all possible file existence etc....
    """
    if singqa is None:
        qlog=qllogger.QLLogger(name="QLConfig")
        log=qlog.getlog()
        log.info("Checking if all the necessary files exist.")

        calib_flavors=['arcs','dark','bias']
        if outconfig["Flavor"]=='science':
            files = [outconfig["RawImage"], outconfig["FiberMap"], outconfig["FiberFlatFile"], outconfig["PSFFile"]]
            for thisfile in files:
                if not os.path.exists(thisfile):
                    sys.exit("File does not exist: {}".format(thisfile))
                else:
                    log.info("File check: Okay: {}".format(thisfile))
        elif outconfig["Flavor"] in calib_flavors:
            files = [outconfig["RawImage"], outconfig["FiberMap"]]
            for thisfile in files:
                if not os.path.exists(thisfile):
                    sys.exit("File does not exist: {}".format(thisfile))
                else:
                    log.info("File check: Okay: {}".format(thisfile))
        elif outconfig["Flavor"]=="flat":
            files = [outconfig["RawImage"], outconfig["FiberMap"], outconfig["PSFFile"]]
            for thisfile in files:
                if not os.path.exists(thisfile):
                    sys.exit("File does not exist: {}".format(thisfile))
                else:
                    log.info("File check: Okay: {}".format(thisfile))
        log.info("All necessary files exist for {} configuration.".format(outconfig["Flavor"]))

    return 


class Palist(object):
    
    """
    Generate PA list and QA list for the Quicklook Pipeline for the given exposure
    """
    def __init__(self,thislist=None,algorithms=None,flavor=None,mode=None):
        """
        thislist: given list of PAs
        algorithms: Algorithm list coming from config file: e.g desispec/data/quicklook/qlconfig_dark.yaml
        flavor: only needed if new list is to be built.
        mode: online offline?
        """
        self.flavor=flavor
        self.mode=mode
        self.thislist=thislist
        self.algorithms=algorithms
        self.palist=self._palist()
        self.qalist=self._qalist()
        qlog=qllogger.QLLogger(name="QLConfig")
        self.log=qlog.getlog()
        
    def _palist(self):
        
        if self.thislist is not None:
            pa_list=self.thislist
        else: #- construct palist
            if self.flavor == 'arcs':
                pa_list=['Initialize','Preproc','BootCalibration','BoxcarExtract','ResolutionFit'] #- class names for respective PAs (see desispec.quicklook.procalgs)
            elif self.flavor == "flat":
                pa_list=['Initialize','Preproc','BoxcarExtract','ComputeFiberflat_QL']
            elif self.flavor == 'bias' or self.flavor == 'dark':
                pa_list=['Initialize','Preproc']
            elif self.flavor == 'science':
                pa_list=['Initialize','Preproc','BoxcarExtract', 'ApplyFiberFlat_QL','SkySub_QL']
            else:
                self.log.warning("Not a valid flavor. Use a valid flavor type to build a palist. Exiting.")
                sys.exit(0)
        self.pamodule='desispec.quicklook.procalgs'
        return pa_list       
    

    def _qalist(self):

        if self.thislist is not None:
            qalist={}
            for PA in self.thislist:
                qalist[PA]=self.algorithms[PA]['QA'].keys()
        else:
            if self.flavor == 'arcs':
                QAs_initial=['Bias_From_Overscan']
                QAs_preproc=['Get_RMS','Count_Pixels']
                QAs_bootcalib=['Calc_XWSigma']
                QAs_extract=['CountSpectralBins']
                QAs_resfit=[]
            elif self.flavor =="flat":
                QAs_initial=['Bias_From_Overscan']
                QAs_preproc=['Get_RMS','Count_Pixels']
                QAs_extract=['CountSpectralBins']
                QAs_computeflat=[]
            elif self.flavor == 'bias' or self.flavor == 'dark':
                QAs_initial=['Bias_From_Overscan']
                QAs_preproc=['Get_RMS','Count_Pixels']
            elif self.flavor =="science":
                QAs_initial=['Bias_From_Overscan']
                QAs_preproc=['Get_RMS','Count_Pixels','Calc_XWSigma']
                QAs_extract=['CountSpectralBins']
                QAs_apfiberflat=['Sky_Continuum','Sky_Peaks']
                QAs_SkySub=['Sky_Residual','Integrate_Spec','Calculate_SNR']
        
            qalist={}
            for PA in self.palist:
                if PA == 'Initialize':
                    qalist[PA] = QAs_initial
                elif PA == 'Preproc':
                    qalist[PA] = QAs_preproc
                elif PA == 'BootCalibration':
                    qalist[PA] = QAs_bootcalib
                elif PA == 'BoxcarExtract':
                    qalist[PA] = QAs_extract
                elif PA == 'ResolutionFit':
                    qalist[PA] = QAs_resfit
                elif PA == 'ComputeFiberflat_QL':
                    qalist[PA] = QAs_computeflat
                elif PA == 'ApplyFiberFlat_QL':
                    qalist[PA] = QAs_apfiberflat
                elif PA == 'SkySub_QL':
                    qalist[PA] = QAs_SkySub
                else:
                    qalist[PA] = None #- No QA for this PA
        self.qamodule='desispec.qa.qa_quicklook'
        return qalist


