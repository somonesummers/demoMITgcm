# TODO: consolidate with helper function in run_config

import os
import MITgcmutils
from IPython.core.debugger import set_trace
#import helper_functions as hf
import re
import imp
import subprocess
import shutil
import numpy as np
import matplotlib.pylab as plt
from datetime import date
from decimal import Decimal
import copy
import pickle 

#imp.reload(hf)
debug_here = set_trace

secsInYr = 3600*24*365

home_dir = '/home/users/earlew/research/MITgcm_python/idealized_channel'

# get exp_dir
with open(os.path.join(home_dir, 'exp_dir_loc.txt'), 'r') as f:
    exp_dir = f.read() #<--update to exp_home_dir

#def update_data():

def backup_file(fpath, cleanup=False):
    
    
    fname = fpath.split('/')[-1]
    fname_dir = "/".join(fpath.split('/')[:-1])
    
    if len(fname.split('.'))==1:
        # for file names with no extension
        fname_root = fname
        ext = ''
    else:
        fname_root = "".join(fname.split('.')[:-1])
        ext = '.' + fname.split('.')[-1]
        #assert len(fname_root)<=1
    
     
    
    print("Backing up %s..." %fname)
    # remove all backup copies except the first one
    dir_list = os.listdir(fname_dir)
    if cleanup:
        #print("Removing up old backups. Original backup will be saved")
        old_back_ups = [f for f in dir_list if f.startswith(fname_root + '_old_v')]
        for fname in old_back_ups:
            if fname not in [fname_root + '_old_v001' + ext]:
                #debug_here()
                os.remove(os.path.join(fname_dir, fname)) 
    
    ver = 1
    
    new_fname = fname_root + '_old_v%03d'%ver + ext
    while new_fname in dir_list:
        ver+=1
        new_fname = fname_root + '_old_v%03d'%ver + ext
  
    #debug_here()
    new_fpath = os.path.join(fname_dir, new_fname)
    shutil.copyfile(fpath, new_fpath) 
    
    #print("file backed up at %s" %new_fpath)

def getLastEndTime(exp_name):
    
    params,_ = hf.get_exp_params_py(exp_name)
    datafile_path = os.path.join(params['input_path'], 'data')
    with open(datafile_path, 'r') as f:
        for line in f:
            if line.startswith(' endTime'):
                r = re.compile(r"\d+")
                LastEndTime = int(r.search(line).group())
                print("LastEndTime: %s secs or %s years" %(LastEndTime, LastEndTime/secsInYr))
                
                return LastEndTime 
                                  
        print("Error: endtime was not found!")
        
        
def getLastIter(exp_name, results_path=''):
    
    #params,_ = hf.get_exp_params_py(exp_name)
    
    if len(results_path)==0:
        results_path = os.path.join(exp_dir, exp_name, 'results')
    
    # get most recent pickup file
    results_files = os.listdir(results_path)
    r = re.compile(r"^(pickup)\.\d+\.data")
    pickup_files = sorted(list(filter(r.search, results_files)))
    
    # get latest iter number
    r2 = re.compile(r"\d+")
    iterN_str = r2.search(pickup_files[-1]).group()
    
    return iterN_str

def getSavedIters(exp_name, vname, results_path=''):
    
    #params,_ = hf.get_exp_params_py(exp_name)
    
    if len(results_path)==0:
        results_path = os.path.join(exp_dir, exp_name, 'results')

    results_files = os.listdir(results_path)
    vble_files = [f for f in results_files if f.startswith(vname) and f.endswith('.data')]
    vble_files.sort()
    time_stmps = np.array([os.path.getmtime(os.path.join(results_path, f)) for f in vble_files])
    time_stmps = time_stmps-time_stmps[0]

    # get iters 
    iters = [int(f.split('.')[1]) for f in vble_files]

    # get model time 
    iters = np.array([int(f.split('.')[1]) for f in vble_files])
    
    return iters

def getClosestIter(exp_name, time_yr, results_path=''):
    
    params,_ = hf.get_exp_params_py(exp_name)
    
    if len(results_path)==0:
        results_path = params['results_path']
    
    # get most recent pickup file
    results_files = os.listdir(results_path)
    
    # method 1 (the right way)
    r = re.compile(r"^(pickup)\.\d+\.data")
    pickup_files1 = sorted(list(filter(r.search, results_files)))
    
    # method 2 (TODO: remove)
    pickup_files2 = sorted([f for f in results_files if f.startswith('pickup.0') and f.endswith('.data')])
    
    
    assert len(pickup_files1) == len(pickup_files2), "lengths don't match. %s vs %s" %(len(pickup_files1), len(pickup_files2))
    
    # extract iters
    iters_strs = [f.split('.')[1] for f in pickup_files1]
    iters = np.array(iters_strs).astype(int)
    
#     print(len(iters_strs[0]))
#     print(iters)
    
    pickup_yrs =  iters*params['deltaT']/secsInYr
    
    idx = np.argmin(np.abs(pickup_yrs-time_yr))
    
    
    time_diff = np.abs(pickup_yrs[idx]-time_yr)
    assert time_diff <1, "No pickup files found within 1 year of year %s. Closest is year %s."%(time_yr, pickup_yrs[idx])
    
    print("Closest pickup_yr=%s, iter_str=%s" %(pickup_yrs[idx], iters_strs[idx]))

    return iters_strs[idx]


def getSavedYrs(exp_name, vname_root='THETA_inst', results_path=''):
    
    
    params,_ = hf.get_exp_params_py(exp_name)
    
    if len(results_path)==0:
        results_path = params['results_path']
    
    # get most recent pickup file
    results_files = os.listdir(results_path)
    
    r = re.compile(r"^(%s)\.\d+\.data"%vname_root)
    results_files = sorted(list(filter(r.search, results_files)))
    
    # extract iters
    iters_strs = [f.split('.')[1] for f in results_files]
    iters = np.array(iters_strs).astype(int)
    
    
    saved_yrs =  iters*params['deltaT']/secsInYr
    
    
    return saved_yrs

    

def createPerturbDir(branch_exp_name, new_exp_name='', branch_year=-1, overwrite=False):
    
    # create new_exp_name if necessary
    params,_ = hf.get_exp_params_py(branch_exp_name)
    exp_dir_root = "/".join(params['exp_path'].split('/')[:-1])
    exp_list = os.listdir(exp_dir_root)
    if len(new_exp_name)==0:
        ver = 1
        new_exp_name = branch_exp_name + '_branch%03d'%ver
        while new_exp_name in exp_list and overwrite==False:
            ver+=1
            new_exp_name = branch_exp_name + '_branch%03d'%ver
            debug_here()
        
    # create new exp directory directory
    new_exp_path = os.path.join(exp_dir, new_exp_name)
    
    # check if directory already 
    if os.path.exists(new_exp_path):
        if overwrite:
            response = input("WARNING %s already exists.\nContinue overwite [yes|no]"%new_exp_path)
            if response=='yes':
                
                # delete existing exp directory
                shutil.rmtree(new_exp_path)
                
                # create new directory
                os.makedirs(new_exp_path, exist_ok=True)
                os.makedirs(os.path.join(new_exp_path, 'results/'), exist_ok=True)
            else:
                print("Exiting...")
                return 
        else:
            print("Error: %s already exists" %new_exp_path)
            print("Supply a new experiment directory name (recommended) or set overwrite=True if that's is desired.")
            print("Exiting...")
            return
    else:
        os.makedirs(new_exp_path, exist_ok=True)
        os.makedirs(os.path.join(new_exp_path, 'results/'), exist_ok=True)
        
    # copy relevant files over to new experiment directory
    print("copying/linking files to new exp directory. Make take a few seconds...")  
    
    shutil.copytree(params['input_path'], os.path.join(new_exp_path, 'input/'))
    shutil.copytree(os.path.join(params['exp_path'], 'build/'),  os.path.join(new_exp_path, 'build/'))
    #os.symlink(os.path.join(params['exp_path'], 'code/'), os.path.join(new_exp_path, 'code/'))
    shutil.copytree(os.path.join(params['exp_path'], 'code/'),  os.path.join(new_exp_path, 'code/'))
    
    if branch_year==-1:
        iter_str = getLastIter(branch_exp_name)
    else:
        iter_str = getClosestIter(branch_exp_name, branch_year)
        
    
    # copy over script files
#     script_files = ['download_from_cluster.sh', 'update_input.sh', 'upload_to_cluster.sh']
#     for fname in script_files:
#         old_fpath = os.path.join(params['exp_path'], fname)
#         new_fpath = os.path.join(new_exp_path, fname)
#         shutil.copyfile(old_fpath, new_fpath)
    
    # copy over results data
    old_results_path = os.path.join(params['exp_path'], 'results')
    new_results_path = os.path.join(new_exp_path, 'results')
    #os.makedirs(new_results_path)
    results_files = ['run.sh', 'run_mitgcm', 
                     'DRF.data', 'hFacS.data', 'hFacW.data', 'hFacC.data',
                     'DRF.meta', 'hFacS.meta', 'hFacW.meta', 'hFacC.meta']
    for fname in results_files:
        old_fpath = os.path.join(old_results_path, fname)
        new_fpath = os.path.join(new_results_path, fname)
        shutil.copyfile(old_fpath, new_fpath) 

    # link datafiles associated with branch_year 
    old_results_files = os.listdir(old_results_path)
    files2link = [f for f in old_results_files if iter_str in f.lower()]
    for fname in files2link:
        old_fpath = os.path.join(old_results_path, fname)
        new_fpath = os.path.join(new_results_path, fname)
        os.symlink(old_fpath, new_fpath)
     
    # TODO: create text file that details where files came from
    source_fname = os.path.join(new_exp_path, 'source_info.txt')
    with open(source_fname, 'w') as f:
        time_stamp = date.today().strftime("%Y-%m-%d-%H:%M")
        info_str = "source_dir: %s\nBranch_iter: %s\nTime_copied: %s" %(params['exp_path'], iter_str, time_stamp)
        f.write(info_str)
        
    print("Done!")
    
    return new_exp_name
        
def extendRun(exp_name, newEndYr, ncpus, dump_config={}, new_diags={}, updateFiles=False, cleanup=False, 
              printNewRunScript=False, newRunName='', results_path='', overwrite_diag=False,
              wall_mins_per_sim_yr=60, update_params_file=True, max_walltime_hrs=48):
    
#     if restartModel==False:
#         print('------------------------------------------------------------------------------------')
#         print("Warning: This is a dry run. New params are copied to temporary files. Model will not restart.")
#         print("To submit run to cluster set restartModel=True")
#         print('------------------------------------------------------------------------------------')
    
    if updateFiles==False:
        print('------------------------------------------------------------------------------------')
        print("Warning: New params and data are copied to temporary files.")
        print("To replace param and data files. set updateFiles=True")
        print('------------------------------------------------------------------------------------')
    
#     if restartModel and updateFiles==False:
#         print("Warning: restartModel==True but updateFiles==False")
#         print("Setting updateFiles==True...")
#         updateFiles=True
    
    newEndTime = newEndYr*secsInYr 
    
    params,_ = hf.get_exp_params_py(exp_name)
    
    if len(newRunName)==0:
        newRunName = exp_name
    
    if len(results_path)==0:
        results_path = params['results_path']
        
    # backup old datafiles and save updated ones.
    backup_file(os.path.join(params['input_path'], 'data'), cleanup=cleanup)
    backup_file(os.path.join(params['input_path'], 'params.p'), cleanup=cleanup)
    backup_file(os.path.join(params['input_path'], 'data.diagnostics'), cleanup=cleanup)
    
    # get most recent pickup file
    results_files = os.listdir(results_path)
    r = re.compile(r"^(pickup)\.\d+\.data")
    pickup_files = sorted(list(filter(r.search, results_files)))
    
    # get latest iter number
    lastIter = getLastIter(exp_name) # this method actually checks the last pickup that's on disk
    newStartIter = int(lastIter)
    oldEndTime = int(newStartIter)*params['deltaT']
    
    #----------------update data file----------------------#
    # update iter0 in data file
    datafile_path = os.path.join(params['input_path'], 'data')
    datafile_path_tmp = os.path.join(params['input_path'], 'data_tmp')
    with open(datafile_path_tmp,'w') as new_data:
        with open(datafile_path, 'r') as old_data:
            for line in old_data:
                if line.startswith(' nIter0=')==False and line.startswith(' endTime=')==False:
                    new_data.write(line)
                else:
                    if line.startswith(' nIter0='):
                        new_data.write(' nIter0=%s,\n'%newStartIter)
                        #debug_here()
                    if line.startswith(' endTime='):
                        #debug_here()
                        r3 = re.compile(r"\d+")
                        #oldEndTime = int(r3.search(line).group())
                        assert newEndTime>=oldEndTime, "NewEndTime=%s must >= oldEndTime=%s"%(newEndTime, oldEndTime)
                        new_data.write(' endTime=%s,\n'%newEndTime)

    #----------------update data.diagnostics dump frequency----------------------#
    diagnostics_dict = parseDataDiagnosticFile(params)
    
    # add new diagnostics
    if len(new_diags)>0 and len(new_diags['fieldName'])>0: 
        jjj = 0
        for ii, fname in enumerate(new_diags['fieldName']):
            
            # check if fileName already exists
            if fname in diagnostics_dict:
                if overwrite_diag==False:
                    print("WARNING: %s is an existing diagnostic. Skipping." %fname)
                    if jjj==0:
                        print("To overwrite set overwrite_diag=True.")
                        jjj =+ 1
                    continue
                else:
                    pass
                    print("Under construction!")
                    # TODO: write code to remove existing file entries
                    debug_here()

            new_diag = {'fieldName': fname, 'fileName': fname, 'frequency': 
                        new_diags['frequency'][ii], 'timePhase': new_diags['timePhase'][ii]}
            diagnostics_dict = append_diag_dict(diagnostics_dict, new_diag)
    
    # update dump_freq (if necessary)
    if len(dump_config)==2:
        #dump_config['fnames'] = ['UVEL_inst', 'VVEL_inst'] 
        #dump_config['freq_d'] = 5
        dump_freq = dump_config['freq_d']*24*3600
        dump_freq_str = "%.8e" %Decimal(dump_freq)
        for fname in dump_config['fnames']:
            for key in diagnostics_dict[fname]:
                if key.startswith('frequency'):
                    freq_str = key
            #print(diagnostics_dict[fname])
            diagnostics_dict[fname][freq_str] = dump_freq_str
            #print(diagnostics_dict[fname])
        
    #debug_here()
    writeDataDiagnosticFile(params, diagnostics_dict)
    
    #----------------update params file----------------------#
    if update_params_file:
        with open(os.path.join(params['input_path'], 'params_tmp.p'), 'wb') as f_new:  
            with open(os.path.join(params['input_path'], 'params.p'), 'rb') as f_old:
                new_params = copy.deepcopy(params)
                new_params['mitgcm_params']['data'][2]['nIter0'] = newStartIter
                new_params['mitgcm_params']['data'][2]['endTime'] = newEndTime
                
                pickle.dump(new_params, f_new)
                
                # TODO: add lines to update dump frequency 
                
                #new_params[
#                 for line in old_params:
#                     if line.startswith('nIter0=')==False and line.startswith('endTime=')==False:
#                         new_params.write(line)
#                     else:
#                         if line.startswith('nIter0='):
#                             new_params.write('nIter0=0;\n') # fix to zero
#                         if line.startswith('endTime='):
#                             new_params.write('endTime=%s;\n'%newEndTime)
                        
    
    #----------------update sbatch file----------------------#
    # update wall clock time
    backup_file(os.path.join(results_path, 'run_mitgcm'), cleanup=cleanup)
    lastIter = getLastIter(exp_name)
    pickupYear = int(lastIter)*params['deltaT']/(3600*24*365)
    nyrs = newEndYr - pickupYear
    walltime_hrs = wall_mins_per_sim_yr*nyrs/60
    
    # update params file
    
    # open run config file
    with open(os.path.join(params['input_path'], 'params.p'), 'rb') as f:
        run_config = pickle.load(f)
        
    params01 = run_config['mitgcm_params']['data'][2]
    params01['nIter0'] = newStartIter
    params01['endTime'] = newEndTime
    
    params_fpath = os.path.join(run_config['input_dir'], 'params.p')

    #ru.backup_file(params_fpath, cleanup=True) # already backed up

    print("saving new params file...")
    with open(params_fpath, 'wb') as f:
        pickle.dump(run_config, f)
    
    

    # pad walltime by 5%
    walltime_hrs = np.ceil(walltime_hrs*1.05)
    print("Estimated runtime: %.1f hrs" %walltime_hrs)
    
    walltime_hrs_str = getWalltimeStr(walltime_hrs, max_walltime_hrs=max_walltime_hrs)

    paramNames = ['#SBATCH -t', 'cd /central', '#SBATCH -J']
    newLines = ['#SBATCH -t %s  # run time (hh:mm:ss)\n' %walltime_hrs_str, '#cd %s\n' %results_path, 
               '#SBATCH -J %s\n' %newRunName] 
 
    runscript_fpath = os.path.join(params['results_path'], 'run_mitgcm')
    #updateFileName = restartModel
    updateLineInText(runscript_fpath, paramNames, newLines, updateFileName=updateFiles, printNewFile=printNewRunScript)
    
    if updateFiles:
        print("updating data and param files!")
        print("Ready to submit...")
        # update data, data.diagnostics, and param files 
        os.rename(os.path.join(params['input_path'], 'data_tmp'), os.path.join(params['input_path'], 'data'))
        os.rename(os.path.join(params['input_path'], 'params_tmp.p'), os.path.join(params['input_path'], 'params.m'))
        if len(dump_config)==2:
            os.rename(os.path.join(params['input_path'], 'data.diagnostics.tmp'), 
                      os.path.join(params['input_path'], 'data.diagnostics'))
    
#     if restartModel:  
#         response = input("Update run files and submit job to cluster? [y]")
#         if response.lower() in ['yes', 'y']:
            
#             # restart model
#             pwd = os.getcwd()
#             os.chdir(results_path)
#             os.chmod('run.sh', 509) # ensure script is executable
#             rc = subprocess.call("./run.sh", shell=True)
#             os.chdir(pwd)
#             if rc!=0:
#                 print("Something bad happened. The job may not have been submitted")
#             else:
#                 print("Job successfully submitted!")
#         else:
#             print("Files were not updated. Run was not submitted cluster.")

def getWalltimeStr(walltime_hrs, min_walltime_hrs=1/60, max_walltime_hrs=48):
    
    if walltime_hrs<min_walltime_hrs:
        print("rounding walltime up to %s minutes" %(min_walltime_hrs*60))
        walltime_hrs = min_walltime_hrs # enforce minimum  
        
    if walltime_hrs>max_walltime_hrs:
        print("Requested walltime exceeds maximum. Capping walltime to %s hours" %max_walltime_hrs)
        walltime_hrs = max_walltime_hrs # enforce maximum 
        
    print('--------------------\n')
    # create walltime str (hh:mm:ss) 
    walltime_str = "%02d:%02d:%02d" %(walltime_hrs//1, (walltime_hrs%1)*60, (((walltime_hrs%1)*60)%1)*60)
    
    return walltime_str
            

def writeDataDiagnosticFile(params, diagnostics_dict):
    
    #backup_file(os.path.join(params['input_path'], 'data.diagnostic'), cleanup=cleanup)
    datafile_path_old = os.path.join(params['input_path'], 'data.diagnostics')
    datafile_path_new = os.path.join(params['input_path'], 'data.diagnostics.tmp')
    
    header=[]
    num_header_lines = 7
    
    footer = [' &\n', '\n', '# Statistics package choices\n', ' &DIAG_STATIS_PARMS\n', ' &\n', '\n']
    
    with open(datafile_path_old, 'r') as f_old:
        for ii, line in enumerate(f_old):
            if ii<num_header_lines: 
                header.append(line)
       
  
    with open(datafile_path_new, 'w') as f_new:
        # write header
        for line in header:
            f_new.write(line)
        
        # write updated diagnostics
        for fname in diagnostics_dict:
            for key in diagnostics_dict[fname]:
                if key.startswith('fields') or key.startswith('fileName'):
                    f_new.write(" %s='%s',\n" %(key, diagnostics_dict[fname][key]))
                else:
                    f_new.write(" %s=%s,\n" %(key, diagnostics_dict[fname][key]))
                    
        # write footer
        for line in footer:
            f_new.write(line)         


def parseDataDiagnosticFile(params):
    
    datafile_path = os.path.join(params['input_path'], 'data.diagnostics')
    
    diagnostic_dict = {}
    #diagnostic_dict['field_names'] = {}

    # get filenames
    with open(datafile_path, 'r') as f: 
        for line in f:
            #print(line.strip())
            line = line.strip()
            if line.strip().startswith("#") or len(line)==0 or '=' not in line:
                continue
            elif line.startswith('fileName'):
                field_name = line.split('=')[-1][1:-2]   
                diagnostic_dict[field_name] = {}
                #debug_here()
            else:
                pass
     
    # get filename data
    field_data = ['fields', 'fileName', 'frequency', 'timePhase']
    with open(datafile_path, 'r') as f: 
        xx = 0
        tmp_dict = {}
        for line in f:
            line = line.strip()
            
            if line.strip().startswith("#") or len(line)==0 or '=' not in line:
                continue
                
            if xx>=4:
                diagnostic_dict[fileName] = tmp_dict
                xx=0
                tmp_dict = {}
            
            for ff in field_data:
                if line.startswith(ff):
                    xx = xx + 1
                    ff2 = line.split('=')[0]
                    if ff in field_data[:2]:
                        tmp_dict[ff2] = line.split('=')[-1][1:-2]  
                        if ff=='fileName':
                            fileName=line.split('=')[-1][1:-2]  
                    else:
                        tmp_dict[ff2] = line.split('=')[-1][:-1] 
     
    # add last line
    diagnostic_dict[fileName] = tmp_dict
    
#     for key in diagnostic_dict:
#         print("%s: %s" %(key,diagnostic_dict[key]))

    return diagnostic_dict

def append_diag_dict(diagnostics_dict, new_diag):
    
    diag_num = len(diagnostics_dict)+1        
     
    diagnostics_dict[new_diag['fileName']] = {}
    diagnostics_dict[new_diag['fileName']]['fields(1, %s)'%diag_num] = new_diag['fieldName']
    diagnostics_dict[new_diag['fileName']]['fileName(%s)'%diag_num] = new_diag['fileName']
    diagnostics_dict[new_diag['fileName']]['frequency(%s)'%diag_num] = new_diag['frequency']
    diagnostics_dict[new_diag['fileName']]['timePhase(%s)'%diag_num] = new_diag['timePhase']
    
    return diagnostics_dict

def updateLineInText(fpath, paramNames, newLines, updateFileName=False, printNewFile=False):
    
    fname = fpath.split('/')[-1]
    fname_dir = "/".join(fpath.split('/')[:-1])
    line_updated = False
    tmp_fpath = os.path.join(fname_dir, 'tmp_file.txt')
    paramNames_tup = tuple(paramNames)
    
    with open(tmp_fpath, 'w') as new_file:  
        with open(fpath, 'r') as old_file:
            for line in old_file:
                if line.startswith(paramNames_tup)==False:
                    new_file.write(line)
                else:
                    # find the specific parameter 
                    xx = np.array([line.startswith(param) for param in paramNames])
                    idx = np.arange(len(paramNames))[xx][0]
                    #debug_here()
                    new_file.write(newLines[idx])
                    line_updated=True
                    
    if line_updated:
        print("Line updated!")
    else:
        print("WARNING: line was not found")
      
    
    if printNewFile:
        with open(tmp_fpath, 'r') as f:
            print('---------------------------------')
            print(f.read())
        
    if updateFileName:
        os.rename(tmp_fpath, fpath)

    return
        
def modify_forcing(exp_name, ywind_s, taux_s, yheat_s, theat_s, Ln=100e3, Ys=250e3, Qaabw=25,
                   add_local_Qaabw=False, save_new_forcing=False, save_comp=True, load_old_forcing=True):
    
    from scipy.interpolate import make_interp_spline, BSpline
    
    fz = 14
    
    if save_new_forcing:
        response = input("Warning: save_new_forcing=True. Overwrite existing forcing? [y/n]")
        if response=='y':
            pass
        else:
            print("Setting save_new_forcing=False...")
            save_new_forcing=False

    #debug_here()
    # load params and forcing
    print("loading params...")
    params, forcing = hf.get_exp_params_py(exp_name, add_hFac=False, add_forcing=load_old_forcing)
    taux = forcing['zonalWind'].mean(axis=1)
    
#     Ln = 100e3 # width of sponge layer. This quantity is not saved anywhere as far as I know
#     Ys = 250e3 # shelf width 
    Ly = np.round(params['Ly'])
    yy = params['yy']
    xx = params['xx']
    
    # define points for spline fit
    # NOTE: This is cumbersome and rather clumsy but I can't think of a better way to adjust these curves.
    # Ideally we would want a simple functional form for the wind profile but in lieu of that, we manually 
    # tinker with the spline fit.
#     tau_east = 0.3
#     tau_west = 0.1
#     taux_s = np.array([-tau_west*0.95, -tau_west*0.95, 0, tau_east*0.68, tau_east*0.68, 0, 0])
#     ywind_s = np.array([0, 250e3, tau0_y*1e3, 1250e3, 2050e3, Ly-Ln, Ly])

#     taux_s = np.array([-tau_west*0.99, -tau_west*0.99,  0, 
#                        tau_east*0.71, tau_east*0.4, 0, 0])
#     ywind_s = np.array([0, 100e3, tau0_y*1e3, 
#                          2010e3, 2190e3, Ly-Ln, Ly])
    
    #debug_here()
    spl = make_interp_spline(ywind_s, taux_s, k=3)  # type: BSpline
    taux_new = spl(yy)
    taux_new[yy>(Ly-Ln)] = 0
    
    # compute curl
    #debug_here()
    dtaux_dy = np.gradient(taux, yy, edge_order=2) 
    dtaux_new_dy = np.gradient(taux_new, yy, edge_order=2) 
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 8))
    plt.sca(axes[0])
    plt.plot(taux, yy/1e3, label='old')
    plt.plot(taux_new, yy/1e3, label='new')
    plt.plot(taux_s, ywind_s/1e3, 'or')
    plt.axvline(0, linestyle='--', linewidth=2, color='k')
    plt.axhline(Ys/1e3, linestyle='-', linewidth=2, color='r')
    plt.xlim(-0.3, 0.3)
    plt.yticks(np.arange(0, 2500, 250))
    plt.xlabel("Wind stress (N/m$^2$)", fontsize=fz)
    plt.grid(True)
    plt.ylim(0, Ly/1e3)
    plt.legend(loc=0, fontsize=fz-2)
    
    plt.sca(axes[1])
    plt.plot(-dtaux_dy*1e7, yy/1e3, ':', linewidth=2, label='old')
    plt.plot(-dtaux_new_dy*1e7, yy/1e3, ':', linewidth=2, label='new')
    plt.xlabel("Wind stress curl (N/m$^3$)", fontsize=fz)
    plt.xlim(-6, 6)
    plt.ylim(0, Ly/1e3)
    plt.grid(True)
    plt.axvline(0, linestyle='--', linewidth=2, color='k')
    plt.axhline(Ys/1e3, linestyle='-', linewidth=2, color='r')
    plt.legend(loc=0, fontsize=fz-2)
    plt.subplots_adjust(hspace=0.3)
    
    fname = os.path.join(params['input_path'], 'old_new_wind.pdf')
    if save_comp:
        if np.all(taux==taux_new):
            response = input("Warning: heat profile save as before. Overwrite existing plots? [y/n]")
            if response=='y':
                plt.savefig(fname, bbox_inches='tight')
            else:
                print("Comparison plot was not saved.")    
    
    #-------------- modify heat profile ----------------#
    qflux_old = forcing['surfQ']
    q_xmean_no_aabw = qflux_old[:, 0]
#     yheat_s = np.array([0, Ys/2, Q0_y*1e3, 1300e3, 2100e3, Ly-Ln, Ly])
#     theat_s = np.array([Qaabw*0.99, Qaabw*0.99, 0, Qacc*0.61, Qacc*0.61, 0, 0])  
    
    spl = make_interp_spline(yheat_s, theat_s, k=3)  # type: BSpline
    qflux_new = spl(yy)
    qflux_new[yy>(Ly-Ln)] = 0 
    
    
    col_new = 'red'
    plt.figure(figsize=(8,8))
    plt.plot(q_xmean_no_aabw, yy/1e3, label='old')
    plt.plot(qflux_new, yy/1e3, label='new', color=col_new)
    plt.plot(theat_s, yheat_s/1e3, 'or')
    plt.axvline(0, linestyle='--', linewidth=2, color='k')
    plt.axhline(Ys/1e3, linestyle='-', linewidth=2, color='r')
    plt.xlim(-10, 20)
    plt.xlabel("Heat flux out of ocean (W/m$^2$)", fontsize=fz)
    plt.grid(True)
    plt.ylim(0, Ly/1e3)
    plt.legend(loc=0, fontsize=fz-2)
    ax1 = plt.gca()
    
    # expand to 2D
    taux_new_2D = np.tile(taux_new[np.newaxis, :], (len(xx), 1))
    qflux_new_2D = np.tile(qflux_new[np.newaxis, :], (len(xx), 1))
    XX, YY = np.meshgrid(xx, yy)
    
    plt.figure(figsize=(12,5))
    vmin, vmax = -20, 20
    clvls = np.arange(vmin, vmax+0.1, 1)
    
    plt.contourf(XX/1e3, YY/1e3, qflux_new_2D.T, clvls, vmin=vmin, vmax=vmax, cmap=plt.cm.rainbow, extend='both')
    plt.axhline(Ys/1e3, linestyle='--', linewidth=1, color='k')
    plt.xticks(np.arange(-2000, 2001, 500))
    ax1.tick_params(axis='both', which='major', labelsize=fz)
    plt.ylabel("Y (km)", fontsize=fz)
    plt.xlabel("X (km)", fontsize=fz)
    plt.title("Ocean surface heat loss (W/m2)", fontsize=fz)
    plt.colorbar(extend='both')
    
    
    if add_local_Qaabw:
        xc = -500e3 
        Lc = 250e3

        qflux_Ys = Qaabw #qflux_new_2D[yi.T].max()
        qflux_abbw = qflux_Ys*np.exp(-((XX.T-xc)/Lc)**2 - ((YY.T+5e3)/(Ys))**2)
        qflux_new_2D = qflux_new_2D+qflux_abbw
        cs = plt.contour(XX/1e3, YY/1e3, qflux_abbw.T, np.arange(50, 201, 50), colors='k', vmin=vmin, vmax=vmax, extend='both')
        #plt.gca().clabel(cs, inline=1, fontsize=fz-4, fmt='%.0f')
        xidx = np.argmin(np.abs(xc-params['xx']))
        qflux_aabw = qflux_new_2D[xidx, :]
        #ax1.set_xlim(-10, Qaabw+10)
        
    
    print("old mean qflux: %.2f" %qflux_old.mean())
    print("new mean qflux: %.2f" %qflux_new_2D.mean())
    
    if save_comp:
        if np.all(q_xmean_no_aabw==qflux_new):
            response = input("Warning: heat profile same as before. Overwrite existing plots? [y/n]")
            if response=='y':
                fname = os.path.join(params['input_path'], 'old_new_heat.pdf')
                plt.savefig(fname, bbox_inches='tight')
            else:
                print("Comparison plot was not saved.")
    
    if save_new_forcing:
        print("saving new forcing")
#         debug_here()
        backup_file(os.path.join(params['input_path'], 'zonalWindFile.bin'), cleanup=True)
        backup_file(os.path.join(params['input_path'], 'surfQfile.bin'), cleanup=True)

        new_taux_fpath = os.path.join(params['input_path'], 'zonalWindFile.bin')
        #debug_here()
        with open(new_taux_fpath, 'wb') as f:
            taux_new_2D.T.astype('>f8').tofile(f)
            
        new_qflux_fpath = os.path.join(params['input_path'], 'surfQfile.bin')
        with open(new_qflux_fpath, 'wb') as f:
            qflux_new_2D.T.astype('>f8').tofile(f)
    
    
def estimateRunTime(nyrs, wall_mins_per_cpu_simyr=-999, ncpus=128, is_hires=False, cluster='HPC'):
    
    # TODO: this is a bad, bad script. To many variables hard coded.
    
    if wall_mins_per_cpu_simyr<=0:
    
        if cluster=='HPC':
            #wall_hrs_per_cpu_nyrs = 55/(128*51)
            if ncpus<=128:
                wall_mins_per_cpu_simyr = 70/128 # maybe no longer accurate
            else:
                if is_hires:
                    wall_mins_per_cpu_simyr = 457/256
                else:
                    wall_mins_per_cpu_simyr = 50/256

        else:
            wall_mins_per_cpu_simyr = 70/(128*1)
    
    est_hrs = wall_mins_per_cpu_simyr*ncpus*nyrs/60
    
    return est_hrs

# def getRunTime(nyrs, wall_mins_per_sim_yr, config_mods={}):
    
#     alpha = compute_alpha(wall_mins_per_sim_yr, config_mods={})
    
#     nTimeSteps = secsInYr/config['dt_secs']
#     nGridPerCpu = (config['Nx']/config['nPx'] +  config['ONx'])*(config['Ny']/config['nPy'] + config['ONy'])*config['Nr']
#     nStepsPerCpu = nTimeSteps*nGridPerCpu
    
#     est_hrs = alpha*nStepsPerCpu*nyrs
    
#     return est_hrs
    
# def compute_alpha(wall_time_mins, sim_time_yr=1, nPx=8, nPy=16):
#     # would nice to not have these hardcoded
#     Nx = 384
#     Ny = 240
#     Nr = 70
#     ONx = 3
#     ONy = 3
#     dt = 617
#     sim_time_secs = sim_time_yr*3600*24*365;
    
#     wall_time_hrs = wall_time_mins/60
#     nTimeSteps = sim_time_secs/dt
#     nGridPerCpu = (Nx/nPx + ONx)*(Ny/nPy + ONy)*Nr
#     nStepsPerCpu = nTimeSteps*nGridPerCpu
    
#     alpha = wall_time_hrs/nStepsPerCpu

#     return alpha

def compute_alpha(wall_mins_per_sim_yr, config_mods={}):
    
    # wall_time_per_sim_yr: time in minutes for model to output a simulation year 
    
    config = {}
    
    config['Nx'] = 384
    config['Ny'] = 240
    config['Nr'] = 70
    config['ONx'] = 3
    config['ONy'] = 3
    config['dt_secs'] = 617
    config['nPx'] = 16
    config['nPy'] = 16
    
    
#     config = {}
    
#     config['Nx'] = 160
#     config['Ny'] = 100
#     config['Nr'] = 70
#     config['ONx'] = 3
#     config['ONy'] = 3
#     config['dt_secs'] = 1629
#     config['nPx'] = 16
#     config['nPy'] = 10

#     config = {}
    
#     config['Nx'] = 80
#     config['Ny'] = 50
#     config['Nr'] = 70
#     config['ONx'] = 3
#     config['ONy'] = 3
#     config['dt_secs'] = 3258
#     config['nPx'] = 8
#     config['nPy'] = 10
    
    # update config
    for key in config_mods:
        if key in config:
            config[key] = config_mods[key]

#     Nx = 384
#     Ny = 240
#     Nr = 70
#     ONx = 3
#     ONy = 3
#     dt = 617
    #sim_time_secs = sim_time_yr*secsInYr;
    
    wall_time_hrs = wall_mins_per_sim_yr/60
    
    nTimeSteps = secsInYr/config['dt_secs']
    nGridPerCpu = (config['Nx']/config['nPx'] +  config['ONx'])*(config['Ny']/config['nPy'] + config['ONy'])*config['Nr']
    nStepsPerCpu = nTimeSteps*nGridPerCpu
    
    alpha = wall_time_hrs/nStepsPerCpu

    return alpha

def copyDataFiles(exp_name1, exp_name2, pickupNum, dryRun=True, overwrite=False):  
    
    if dryRun:
        print("WARNING: This is only a dry run. set dryRun=False to copy files")
        print("---------------------------------------")
    
    params1, _ = hf.get_exp_params_py(exp_name1)
    params2, _ = hf.get_exp_params_py(exp_name2)
    
    pickupStr = "%010d" %pickupNum
    
    results_files1 = os.listdir(params1['results_path'])
    #data_files = [fname for fname in results_files1 if pickupStr in fname]
    results_files2 = os.listdir(params2['results_path'])

    vbles2copy = ['THETA', 'THETA_inst', 'UVEL', 'UVEL_inst', 'VVEL', 'VVEL_inst', 'WVEL', 
                 'LaHs1TH', 'LaHw1TH', 'LaUH1TH', 'LaVH1TH', 'PHIHYD']
    ext_list = ['data', 'meta']     
    
    for vname in vbles2copy:
        for ext in ext_list:
            fname = ".".join([vname, pickupStr, ext])
            old_fpath = os.path.join(params1['results_path'], fname)
            new_fpath = os.path.join(params2['results_path'], fname)

            if fname in results_files1:
                if fname in results_files2 and overwrite==False:
                    print("%s already exists in target dir. Skipping..." %(fname))
                    continue
                else:
                    print("copying %s to target dir..." %(fname))
                    if dryRun==False:
                        shutil.copyfile(old_fpath, new_fpath) 
            else:
                print("Warning: %s not found in old path. Skipping..."% fname)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    