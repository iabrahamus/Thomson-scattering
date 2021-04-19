#!/usr/bin/python

from Pegasus.DAX3 import *
from optparse import OptionParser
import sys
import textwrap
import os
import math
import time
import stat


class ThomsonScattering:
    """Contains methods to generate a Pegasus workflow for high-energy scattering of photons from non-linear Thomson scattering"""

    def __init__(self):
        
        self.input_filenames = {"ephase" : "ephase_space.dat",
                                "params" : "parameters.dat"
                                     }
        
        self.output_file_basenames = {"powspec_inc" : "powspec_inc.dat",
					                   "e_final_traj" : "e_final_tra.dat",
					                   "powspec_coh" : "powspec_coh.dat"}
        
        
        self.output_file_format = {"powspec_inc" : {'rows': 1000, 'columns' : 2},
                                   "powspec_coh" : {'rows': 1000, 'columns' : 7},
                                   "efield" : {'rows' : 30000, 'columns' : 5},
                                   "traj" : {'rows' : 30000, 'columns' : 5}
                                   }
        
        
        self.default_options = {"electrons_per_job" : 25,
                                "total_electrons" : 100000, 
                                "main_exe" : "traject.x",
                                "combine_exe" : "combine.x",
                                "output_dir" : "OUTPUT",
                                "temp_dir" : "tmp"
                                }

        self.param_names = {"numf" : "Number_of_frequencies_VAR_NUMF",
                       "numtheta" : "Number_of_polar_angles_VAR_NTHETA",
                       "numphi" : "Number_azimuthal_angles_VAR_NPHI"
                       }
        
        self.param_values = {"numf" : "",
                             "numtheta" : "",
                             "numphi" : ""
        		    }

        self.site_tmpdirs = {  "local" : {"shared-scratch" : "work"},
				"local-pbs" : {"shared-scratch" : "shared-scratch"}
			    }


	self.dax_name = "thomson_scattering.dax"
        self.sites_xml_name = "sites.xml"
        self.pegasusrc_name = "pegasusrc"
        self.submitscript_name = "submit.sh"
	

    def _parseOptions(self):
        usage = "usage: %prog [options] working_directory"
        description = """\
                    Generates Pegasus workflow for high-energy scattering of photons from non-linear Thomson scattering.
                    The working directory must contain the compiled executables and input files.
                    Required input files are: """ + ', '.join([infile for infile in self.input_filenames.values()]) + '.'
        
        parser = OptionParser(usage=usage,description=textwrap.dedent(description))
        
        parser.add_option("-p","--perjob",dest="electrons_per_job",help="Number of electrons per job (defaults to %default).",default=self.default_options['electrons_per_job'],type="int")
        parser.add_option("-e","--totalelectrons",dest="total_electrons",help="Total number of electrons (defaults to %default).",default=self.default_options['total_electrons'],type="int")
        parser.add_option("-m","--mainexe",dest="main_exe",help="Name of the main executable (defaults to %default).",default=self.default_options['main_exe'])
        parser.add_option("-c","--combineexe",dest="combine_exe",help="Name of the executable to combine the results (defaults to %default).",default=self.default_options['combine_exe'])
        parser.add_option("-o","--output",help="Directory where final job output will be placed (defaults to '%default' in working directory).",dest="output_dir",type="string",default=self.default_options['output_dir'])
        parser.add_option("-t","--temp",dest="temp_dir",help="Directory used by Pegasus to store workflow-related temporary files (defaults to '%default' in working directory).",default=self.default_options['temp_dir'])
        
        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit(0)
            
        (options,args) = parser.parse_args()        
        
        if len(args) !=1:
            parser.print_help()
            sys.exit(1)
            
        working_dir = args[0]
        
        return (working_dir,options,args)
    
 
    
    def _verifyOptions(self,working_dir,options,args): 
        if os.path.isdir(os.path.normpath(working_dir)):
            self.run_home = os.path.abspath(working_dir)    
        else:
            print "Working directory '%s' does not exist" % working_dir
            sys.exit(1)

            
        if os.path.isfile(os.path.join(self.run_home,options.main_exe)):
            self.main_exe = os.path.join(self.run_home,options.main_exe)
        else:
            print "'%s' not found in directory '%s'." % (options.main_exe,self.run_home)
            sys.exit(1)
        
        if os.path.isfile(os.path.join(self.run_home,options.combine_exe)):
            self.combine_exe = os.path.join(self.run_home,options.combine_exe)
        else:
            print "'%s' not found in directory '%s'." % (options.combine_exe,self.run_home)
            sys.exit(1)
         
        for (key,filename) in self.input_filenames.items():
            if not os.path.isfile(os.path.join(self.run_home,filename)):
                print "Input file '%s' not found in '%s'." % (filename,self.run_home)
                sys.exit(1)
        
        if options.output_dir == self.default_options['output_dir']:
            self.output_dir = os.path.join(self.run_home,self.default_options['output_dir'])
        else:
            self.output_dir = os.path.abspath(options.output_dir)
        
        if not os.path.isdir(self.output_dir):
            try:
                os.makedirs(self.output_dir)
            except:
                print "Cannot make output directory '%s'." % self.output_dir 
        
        if options.temp_dir == self.default_options['temp_dir']:
            self.temp_dir = os.path.join(self.run_home,self.default_options['temp_dir'])
        else:
            self.temp_dir = os.path.abspath(options.temp_dir)

	self.work_dir = os.path.join(self.temp_dir,self.site_tmpdirs['local']['shared-scratch'])
 	self.sharedscratch_dir = os.path.join(self.temp_dir,self.site_tmpdirs['local-pbs']['shared-scratch'])
        
        if not os.path.isdir(self.temp_dir):
            try:
                os.makedirs(self.temp_dir)
            except:
                print "Cannot make temp directory '%s'." % self.temp_dir 

	if not os.path.isdir(self.work_dir):
	    try:
		os.makedirs(self.work_dir)
	    except:
		print "Cannot make work directory '%s'." % self.work_dir
	
	if not os.path.isdir(self.sharedscratch_dir):
	    try:
		os.makedirs(self.sharedscratch_dir)
	    except:
		print "Cannot make shared scratch directory '%s'." % self.sharedscratch_dir	

 
        self.electrons_per_job = options.electrons_per_job
        self.total_electrons = options.total_electrons
        self.num_jobs = int(math.ceil(float(self.total_electrons)/float(self.electrons_per_job)))
        self.job_id = "%s_%d" % (os.path.basename(self.run_home),int(time.time()))
    
    def _setupRun(self,dag):      
        self.combine_jobs = {}
        #=======================================================================
        # for key in self.output_file_basenames.keys():
        #    self.combine_jobs[key] = Job(name=self.exes['combine'],namespace=self.job_id)
        #    self.combine_jobs[key].addArguments(str(self.num_elements),self.files['listfiles'][key].name,self.output_file_basenames[key])
        #    self.combine_jobs[key].uses(self.files['listfiles'][key],link=Link.INPUT)
        #    self.combine_jobs[key].uses(self.output_file_basenames[key],link=Link.OUTPUT,transfer=True)
        #    dag.addJob(self.combine_jobs[key])            
        #=======================================================================
        
        self.combine_jobs['powspec_inc'] = Job(name=self.exes['combine'],namespace=self.job_id)
        self.combine_jobs['powspec_inc'].addArguments('incoherent',str(self.num_elements),self.files['listfiles']['powspec_inc'].name,self.output_file_basenames['powspec_inc'])
        self.combine_jobs['powspec_inc'].uses(self.files['listfiles']['powspec_inc'],link=Link.INPUT)
        self.combine_jobs['powspec_inc'].uses(self.output_file_basenames['powspec_inc'],link=Link.OUTPUT,transfer=True)
        dag.addJob(self.combine_jobs['powspec_inc'])
        
        self.combine_jobs['powspec_coh'] = Job(name=self.exes['combine'],namespace=self.job_id)
        self.combine_jobs['powspec_coh'].addArguments('coherent',str(self.num_elements),self.files['listfiles']['powspec_coh'].name,self.output_file_basenames['powspec_coh'])
        self.combine_jobs['powspec_coh'].uses(self.files['listfiles']['powspec_coh'],link=Link.INPUT)
        self.combine_jobs['powspec_coh'].uses(self.output_file_basenames['powspec_coh'],link=Link.OUTPUT,transfer=True)
        dag.addJob(self.combine_jobs['powspec_coh'])        
        
        self.combine_jobs['e_final_traj'] = Job(name=self.exes['combine'],namespace=self.job_id)
        self.combine_jobs['e_final_traj'].addArguments('trajectory',str(self.electrons_per_job),self.files['listfiles']['e_final_traj'].name,self.output_file_basenames['e_final_traj'])
        self.combine_jobs['e_final_traj'].uses(self.files['listfiles']['e_final_traj'],link=Link.INPUT)
        self.combine_jobs['e_final_traj'].uses(self.output_file_basenames['e_final_traj'],link=Link.OUTPUT,transfer=True)
        dag.addJob(self.combine_jobs['e_final_traj'])        
        
        
        self.main_jobs = {}
        for job_number in range(1,self.num_jobs):
            start = str(((job_number-1)*self.electrons_per_job)+1)
            stop = str(job_number*self.electrons_per_job)
            
            self.main_jobs[job_number]= Job(name=self.exes['main'],namespace=self.job_id)
            self.main_jobs[job_number].addArguments(start,stop)
            self.main_jobs[job_number].uses(self.files['ephase'],link=Link.INPUT)
            self.main_jobs[job_number].uses(self.files['params'],link=Link.INPUT)
            dag.addJob(self.main_jobs[job_number])
            for key in self.output_file_basenames.keys():
                self.main_jobs[job_number].uses(self.files[key][job_number],link=Link.OUTPUT,transfer=False)
                self.combine_jobs[key].uses(self.files[key][job_number],link=Link.INPUT)
                dag.depends(parent=self.main_jobs[job_number],child=self.combine_jobs[key])
            
            
        self.main_jobs[self.num_jobs] = Job(name=self.exes['main'],namespace=self.job_id)
        self.main_jobs[self.num_jobs].addArguments(str(int(stop)+1),str(self.total_electrons))
        self.main_jobs[self.num_jobs].uses(self.files['ephase'],link=Link.INPUT)
        self.main_jobs[self.num_jobs].uses(self.files['params'],link=Link.INPUT)
        dag.addJob(self.main_jobs[self.num_jobs])
        for key in self.output_file_basenames.keys():
                self.main_jobs[self.num_jobs].uses(self.files[key][self.num_jobs],link=Link.OUTPUT,transfer=False)
                self.combine_jobs[key].uses(self.files[key][self.num_jobs],link=Link.INPUT)
                dag.depends(parent=self.main_jobs[self.num_jobs],child=self.combine_jobs[key])
        
        
      
            
    def _addFiles(self,dag):
        self.files = {}
        for key in self.input_filenames.keys():
            self.files[key] = File(self.input_filenames[key])
            self.files[key].addPFN(PFN("file://"+os.path.join(self.run_home,self.input_filenames[key]),"local-pbs"))
            dag.addFile(self.files[key])
            
        
        self.files['listfiles'] = {}
        for key in self.output_file_basenames.keys():
            self.files[key] = {}
            self.files['listfiles'][key] = File("%s.list" % key)
            self.files['listfiles'][key].addPFN(PFN("file://"+os.path.join(self.run_home,'%s.list' % key),"local-pbs"))

        for job_number in range(1,self.num_jobs):
            start = str(((job_number-1)*self.electrons_per_job)+1)
            stop = str(job_number*self.electrons_per_job)
            for key in self.output_file_basenames.keys():
                self.files[key][job_number] = File("%s.%s.%s" % (self.output_file_basenames[key],start,stop))
        
        for key in self.output_file_basenames.keys():
            self.files[key][self.num_jobs] = File("%s.%s.%s" % (self.output_file_basenames[key],str(int(stop)+1),self.total_electrons))
            
        for key1 in self.output_file_basenames.keys():
            dag.addFile(self.files['listfiles'][key1])
            for key2 in self.files[key1].keys():
                dag.addFile(self.files[key1][key2])

    def _parseParams(self,run_home):
        for line in open(os.path.join(run_home,self.input_filenames['params'])):
            for key in self.param_names.keys():
                if line.split()[1] == self.param_names[key]:
                    self.param_values[key] = int(float(line.split()[0]))
                    
        self.num_elements = 1
        for key in self.param_values.keys():
            print "%s : %d" % (key,self.param_values[key])
            self.num_elements *= self.param_values[key]
                    
        print "total elements is %d" % (self.num_elements)
        
    def _createDag(self):
        thomson_wf = ADAG(self.job_id)
        return thomson_wf
        
        
    def _addExecutables(self,dag):
        self.exes = {}
        self.exes['main'] = Executable(name=os.path.basename(self.main_exe),namespace=self.job_id,os='linux',arch='x86_64',installed=False)
        self.exes['main'].addPFN(PFN("file://"+self.main_exe,site="local-pbs"))
        
        self.exes['combine'] = Executable(name=os.path.basename(self.combine_exe),namespace=self.job_id,os='linux',arch='x86_64',installed=False)
        self.exes['combine'].addPFN(PFN("file://"+self.combine_exe,site="local-pbs"))
    
        for key in self.exes.keys():
            dag.addExecutable(self.exes[key])

    def _writeListFiles(self):
        for name in self.files['listfiles'].keys():
            list_fh = open(os.path.join(self.run_home,self.files['listfiles'][name].name),'w') 
            for name2 in self.files[name].keys():
                list_fh.write('%s\n' % self.files[name][name2].name)
            list_fh.close()



    def _createSiteFile(self,run_home,work_dir,output_dir,sitesxml_name,sharedscratch_dir):
        sitesxmltxt = """\
        <?xml version="1.0" encoding="UTF-8"?>
        <sitecatalog xmlns="http://pegasus.isi.edu/schema/sitecatalog" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://pegasus.isi.edu/schema/sitecatalog http://pegasus.isi.edu/schema/sc-4.0.xsd" version="4.0">
        <site handle="local" arch="x86_64" os="LINUX">
                <directory type="shared-scratch" path="%s">
                        <file-server operation="all" url="file://%s" />
                </directory>
                <directory type="local-storage" path="%s">
                        <file-server operation="all" url="file://%s" />
                </directory>
        </site>

        <site  handle="condorpool" arch="x86_64" os="LINUX">
            <profile namespace="pegasus" key="style" >condor</profile>
            <profile namespace="condor" key="requirements" >(Target.Arch == "X86_64")</profile>
        </site>
	
        <site  handle="local-pbs" arch="x86_64" os="LINUX">
                <directory type="shared-scratch" path="%s">
                        <file-server operation="all" url="file://%s"/>
                </directory>

        <profile namespace="env" key="PEGASUS_HOME">/usr</profile>

        <profile namespace="pegasus" key="style" >glite</profile>
        <profile namespace="pegasus" key="change.dir">true</profile>

        <profile namespace="condor" key="grid_resource">pbs</profile>
        <profile namespace="condor" key="batch_queue">guest</profile>
        <profile namespace="globus" key="maxwalltime">240</profile>
        </site>
        </sitecatalog>""" \
        % (work_dir,work_dir,output_dir,output_dir,sharedscratch_dir,sharedscratch_dir)
        
        wrapper_fh = open(os.path.join(run_home,sitesxml_name),'w')
        wrapper_fh.write(textwrap.dedent(sitesxmltxt))
        wrapper_fh.close()

    def _createPegasusRc(self,run_home,rcfilename,sitesxml_name):
        pegasusrctxt = """\
        pegasus.catalog.site=XML
        pegasus.catalog.site.file=%s
        pegasus.dir.storage.deep=false
        pegasus.condor.arguments.quote  false
        pegasus.data.configuration=sharedfs """ % (sitesxml_name)
        
        wrapper_fh = open(os.path.join(run_home,rcfilename),'w')
        wrapper_fh.write(textwrap.dedent(pegasusrctxt))
        wrapper_fh.close()

    def _createSubmitScript(self,run_home,submit_filename,rcfilename,temp_dir,dax_name):
        submit_txt = """\
        #!/bin/bash
        pegasus-plan --conf %s --sites local-pbs --dir %s --output-site local --dax %s --randomdir --submit
        """ % (rcfilename,temp_dir,dax_name)
        
        wrapper_fh = open(os.path.join(run_home,submit_filename),'w')
        wrapper_fh.write(textwrap.dedent(submit_txt))
        wrapper_fh.close()
        os.chmod(os.path.join(run_home,submit_filename),0755)
        
        
        


    
if __name__ == '__main__':
    job = ThomsonScattering()
    (working_dir,options,args) = job._parseOptions()
    job._verifyOptions(working_dir, options, args)
    job._parseParams(job.run_home)
   
    print "My options are:"
    print "\tElectrons per job: %s" % job.electrons_per_job
    print "\tTotal electrons: %s" % job.total_electrons
    print "\tTotal number of jobs: %s" % job.num_jobs
    #print "\tMain exe: %s" % job.main_exe
    #print "\tCombine exe: %s" % job.combine_exe
    print "\tOutput dir: %s" % job.output_dir
    print "\tJob ID: %s" % job.job_id
    
    thomson_wf = job._createDag()
    job._addExecutables(thomson_wf)
    job._addFiles(thomson_wf)
    job._setupRun(thomson_wf)
    job._writeListFiles()
    job._createSiteFile(job.run_home, job.work_dir, job.output_dir, job.sites_xml_name, job.sharedscratch_dir)
    job._createPegasusRc(job.run_home,job.pegasusrc_name,job.sites_xml_name)
    job._createSubmitScript(job.run_home,job.submitscript_name,job.pegasusrc_name,job.work_dir,job.dax_name)
    thomson_wf.writeXMLFile(os.path.join(job.run_home,job.dax_name))
    
    print
    print "Job setup complete."
    print "Change directory to '%s' and run './%s'" % (job.run_home,job.submitscript_name)
    
    
