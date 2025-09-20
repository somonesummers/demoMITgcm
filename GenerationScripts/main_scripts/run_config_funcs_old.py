import numpy as np
import os
import re
import MITgcmutils # see https://mitgcm.readthedocs.io/en/latest/utilities/utilities.html
import pickle
import shutil
import run_analysis_funcs as raf
import warnings
import copy
            
def write_data(run_config, data_params, group_name='data', prec=8, lf='\r\n'):
    
    """
    TODO:
    - add doc file
    - replace run_config with input file path
    """
    group_name2 = group_name.upper() + "_"
    file_header =  group_name.upper()  
    fname = "data" + '.' + group_name
        
    if group_name=='data':
        section_titles = ['Continuous equation parameters', 'Elliptic solver parameters',
                          'Time stepping parameters', 'Gridding parameters', 
                          'Input datasets']
        group_name2 = ''
        file_header = 'Model'
        fname = 'data'
        
    elif group_name=='rbcs':
        section_titles = ['Relaxtion parameters', 'ptracer parameters']

        
    elif group_name=='obcs':
        section_titles = ['Open boundary types and locations', 'Orlanski parameters',
                          'Sponge-layer parameters', 'Stevens parameters']
        
    elif group_name=='gmredi':
        section_titles = ['Parameters for gmredi package']
        group_name2 = 'GM_'
     
    elif group_name=='layers':
        section_titles = ['Parameters for layers package']
        
    elif group_name=='diagnostics':
        section_titles = ['Diagnostics package choices', 'Statistics package choices']

    else:
        print("Error: unsupported group type")
        return
    

    input_fpath = os.path.join(run_config['run_dir'], 'input')
    fpath = os.path.join(input_fpath, fname)
                         
    with open(fpath, 'w') as f:
        
        header_main = "| %s parameters |" %file_header
        header_pad = "="*len(header_main)
        header = """# %s %s# %s %s# %s %s%s"""%(header_pad, lf, header_main, lf, header_pad, lf, lf)
        f.write(header)
        
        for ii, params in enumerate(data_params):
            
            f.write("# %s%s" %(section_titles[ii], lf))
            
            if group_name=='diagnostics':
                if ii==0:
                    f.write(" &DIAGNOSTICS_LIST"+lf)
                else:
                    f.write(" &DIAG_STATIS_PARMS"+lf)  
            else:
                f.write(" &%sPARM%02d" %(group_name2, ii+1)+lf)
            
            for key in params:
                val = params[key]
                
                if type(val)==bool:
                    val_str = ".%s."%str(val).upper()
                    
                elif type(val)==float or type(val)==np.float64:
                    val_str = ('{:.%ie}'%prec).format(val)
                    
                elif type(val)==np.ndarray:
                    
                    val_str = array2str(val, prec=8, maxCol=10)
                    
                elif type(val)==str:
                    val_str = "'%s'"%val
                    
                else:
                    val_str = str(val)
                
                f.write(' %s=%s,%s'%(key, val_str, lf))
                
            f.write(' &%s' %lf)
            f.write(lf)
            
    
def array2str(arr, numtype='float', prec=8, maxCol=10, lf='\r\n'):
    
    """
    Formats a 1d array of numbers as an input parameter for MITgcm.
    
    This code mimics the functionality of list2str.m written by Andrew Stewart.
    
    
    """

    arr_str_list = []
    ii=1
    for item in arr: 
        if type(item.flat[0])==np.float64:
            arr_str_list.append(('{:.%ie}'%prec).format(item))
        else:
            arr_str_list.append(str(item))
            
        if ii<len(arr):
            arr_str_list.append(', ')
            
            if ii%maxCol==0 and ii<len(arr):
                arr_str_list.append('%s      ' %lf)

        ii+=1
     
    arr_str = "".join(arr_str_list)
   
    return arr_str

    
def createLAYERSSIZEh(run_config, Nlayers, layers_maxNum=1):


    layerssizehtext = ['C ======================================================================\n',
                     'C * Compiled-in size options for the LAYERS package *\n',
                     'C\n',
                     'C  - Just as you have to define Nr in SIZE.h, you must define the number\n',
                     'C    of vertical layers for isopycnal averaging so that the proper array\n',
                     'C    sizes can be declared in the LAYERS.h header file.\n',
                     'C\n',
                     'C  - Variables -\n',
                     'C        NLayers      :: the number if isopycnal layers (must match data.layers)\n',
                     'C        FineGridFact :: how many fine-grid cells per dF cell\n',
                     'C        FineGridMax  :: the number of points in the finer vertical grid\n',
                     'C                         used for interpolation\n',
                     '      INTEGER    Nlayers, FineGridFact, FineGridMax, layers_maxNum\n',
                     '      PARAMETER( Nlayers = %s )\n'%Nlayers,
                     '      PARAMETER( FineGridFact = 10 )\n',
                     '      PARAMETER( FineGridMax = Nr * FineGridFact )\n',
                     '      PARAMETER( layers_maxNum = %s )\n'%layers_maxNum,
                     ' \n']

    code_fpath = os.path.join(run_config['run_dir'], 'code')
    fpath = os.path.join(code_fpath, 'LAYERS_SIZE.h')

    with open(fpath, 'w') as f:

        text_str = "".join(layerssizehtext)
        f.writelines(text_str)

            

def createDIAGSIZEh(run_config, numdiags, numLevels, lf='\r\n'):
    
    diagsizehtext = ['C     Diagnostics Array Dimension ',
                     'C     --------------------------- ',
                     'C     ndiagMax   :: maximum total number of available diagnostics ',
                     'C     numlists   :: maximum number of diagnostics list (in data.diagnostics) ',
                     'C     numperlist :: maximum number of active diagnostics per list (data.diagnostics) ',
                     'C     numLevels  :: maximum number of levels to write    (data.diagnostics) ',
                     'C     numDiags   :: maximum size of the storage array for active 2D/3D diagnostics ',
                     'C     nRegions   :: maximum number of regions (statistics-diagnostics) ',
                     'C     sizRegMsk  :: maximum size of the regional-mask (statistics-diagnostics) ',
                     'C     nStats     :: maximum number of statistics (e.g.: aver,min,max...) ',
                     'C     diagSt_size:: maximum size of the storage array for statistics-diagnostics ',
                     'C Note : may need to increase "numDiags" when using several 2D/3D diagnostics, ',
                     'C  and "diagSt_size" (statistics-diags) since values here are deliberately small. ',
                     '      INTEGER    ndiagMax ',
                     '      INTEGER    numlists, numperlist, numLevels ',
                     '      INTEGER    numDiags ',
                     '      INTEGER    nRegions, sizRegMsk, nStats ',
                     '      INTEGER    diagSt_size ',
                     '      PARAMETER( ndiagMax = 500 ) ',
                     '      PARAMETER( numlists = %s, numperlist = 1, numLevels=%s ) '%(numdiags, numLevels),
                     '      PARAMETER( numDiags = %s*%s ) '%(numdiags, numLevels),
                     '      PARAMETER( nRegions = 0 , sizRegMsk = 1 , nStats = 4 ) ',
                     '      PARAMETER( diagSt_size = 9*Nr ) ',
                     ' \n',
                     ' \n',
                     'CEH3 ;;; Local Variables: *** ',
                     'CEH3 ;;; mode:fortran *** ',
                     ' %s'%lf]
    
    input_fpath = os.path.join(run_config['run_dir'], 'code')
    fpath = os.path.join(input_fpath, 'DIAGNOSTICS_SIZE.h')

    with open(fpath, 'w') as f:

        text_str = lf.join(diagsizehtext)
        f.writelines(text_str)

    
def create_eedata(run_config, nTx, nTy, lf='\r\n'):
    # Main function of eedata file is to set the number of threads per
    # process in the x- and y- directions. Typically nTx=nSx and nTy=nSy,
    # so that each process handles nSx*nSy tiles, each with its own thread.
    
    listterm = '&'
    
    eedatatext = ['# "eedata" file ',
                  '# Lines beginning "#" are comments ',
                  '# nTx - No. threads per process in X ',
                  '# nTy - No. threads per process in Y ',
                  ' &EEPARMS ',
                  ' nTx=%s, '%(nTx),
                  ' nTy=%s, '%(nTy),
                  ' %s '%listterm,
                  '# Note: Some systems use & as the namelist terminator (as shown here). ',
                  '#       Other systems use a / character.']
    
    input_fpath = os.path.join(run_config['run_dir'], 'input')
    fpath = os.path.join(input_fpath, 'eedata')

    with open(fpath, 'w') as f:

        text_str = lf.join(eedatatext)
        f.writelines(text_str)


def createSIZEh(run_config, grid_params):
    
    # update to control line break 
    
    sizehtext = ['CBOP \r\n',
                 'C    !ROUTINE: SIZE.h \r\n',
                 'C    !INTERFACE: \r\n',
                 'C    include SIZE.h \r\n',
                 'C    !DESCRIPTION: \\bv \r\n',
                 'C     *==========================================================* \r\n',
                 'C     | SIZE.h Declare size of underlying computational grid.      \r\n',
                 'C     *==========================================================* \r\n',
                 'C     | The design here support a three-dimensional model grid     \r\n',
                 'C     | with indices I,J and K. The three-dimensional domain       \r\n',
                 'C     | is comprised of nPx*nSx blocks of size sNx along one axis  \r\n',
                 'C     | nPy*nSy blocks of size sNy along another axis and one      \r\n',
                 'C     | block of size Nz along the final axis.                     \r\n',
                 'C     | Blocks have overlap regions of size OLx and OLy along the  \r\n',
                 'C     | dimensions that are subdivided.                            \r\n',
                 'C     *==========================================================* \r\n',
                 'C     \\ev \r\n',
                 'CEOP \r\n',
                 'C     Voodoo numbers controlling data layout. \r\n',
                 'C     sNx :: No. X points in sub-grid. \r\n',
                'C     sNy :: No. Y points in sub-grid. \r\n',
                'C     OLx :: Overlap extent in X. \r\n',
                'C     OLy :: Overlat extent in Y. \r\n',
                'C     nSx :: No. sub-grids in X. \r\n',
                'C     nSy :: No. sub-grids in Y. \r\n',
                'C     nPx :: No. of processes to use in X. \r\n',
                'C     nPy :: No. of processes to use in Y. \r\n',
                'C     Nx  :: No. points in X for the total domain. \r\n',
                'C     Ny  :: No. points in Y for the total domain. \r\n',
                'C     Nr  :: No. points in Z for full process domain. \r\n',
                '      INTEGER sNx \r\n',
                '      INTEGER sNy \r\n',
                '      INTEGER OLx \r\n',
                '      INTEGER OLy \r\n',
                '      INTEGER nSx \r\n',
                '      INTEGER nSy \r\n',
                '      INTEGER nPx \r\n',
                '      INTEGER nPy \r\n',
                '      INTEGER Nx \r\n',
                '      INTEGER Ny \r\n',
                '      INTEGER Nr \r\n',
                '      PARAMETER ( \r\n',
                '     &           sNx =  %s, \r\n'%grid_params['sNx'],
                '     &           sNy =  %s, \r\n'%grid_params['sNy'],
                '     &           OLx =  %s, \r\n'%grid_params['OLx'],
                '     &           OLy =  %s, \r\n'%grid_params['OLy'],
                '     &           nSx =  %s, \r\n'%grid_params['nSx'],
                '     &           nSy =  %s, \r\n'%grid_params['nSy'],
                '     &           nPx =  %s, \r\n'%grid_params['nPx'],
                '     &           nPy =  %s, \r\n'%grid_params['nPy'],
                '     &           Nx  = sNx*nSx*nPx, \r\n',
                '     &           Ny  = sNy*nSy*nPy, \r\n',
                '     &           Nr  =  %s) \r\n'%grid_params['Nr'],
                ' \r\n',
                'C     MAX_OLX :: Set to the maximum overlap region size of any array \r\n',
                'C     MAX_OLY    that will be exchanged. Controls the sizing of exch \r\n',
                'C                routine buffers. \r\n',
                '      INTEGER MAX_OLX \r\n',
                '      INTEGER MAX_OLY \r\n',
                '      PARAMETER ( MAX_OLX = OLx, \r\n',
                '     &            MAX_OLY = OLy ) \r\n',
                ' \r\n']
    
    code_fpath = os.path.join(run_config['run_dir'], 'code')
    fpath = os.path.join(code_fpath, 'SIZE.h')

    with open(fpath, 'w') as f:

        text_str = "".join(sizehtext)
        f.writelines(text_str)    
    
def createSBATCHfile_Sherlock(run_config, cluster_params, walltime_hrs, email, 
                              min_walltime_hrs=2/60, max_walltime_hrs=48, mem_GB=1, 
                              part='serc', use_mpirun=False, specify_nodes=True):

    """
    function to generate sbatch code for slurm submission
    """

    if walltime_hrs<min_walltime_hrs:
        print("rounding walltime up to %s minutes" %(min_walltime_hrs*60))
        walltime_hrs = min_walltime_hrs # enforce minimum  
        
    if walltime_hrs>max_walltime_hrs:
        print("Requested walltime exceeds maximum. Capping walltime to %s hours" %max_walltime_hrs)
        walltime_hrs = max_walltime_hrs # enforce maximum 
        
    print('--------------------\n')
    # create walltime str (hh:mm:ss) 
    walltime_str = "%02d:%02d:%02d" %(walltime_hrs//1, (walltime_hrs%1)*60, (((walltime_hrs%1)*60)%1)*60)
    
    ncpus = run_config['ncpus_xy'][0]*run_config['ncpus_xy'][1]
    results_fpath = os.path.join(run_config['run_dir'], 'results')
    
    
    if use_mpirun:
        run_cmd = 'mpirun -np %s ./mitgcmuv'%(ncpus)
    else:
        # default to srun
        run_cmd = 'srun ./mitgcmuv'
       
    # compute the minimum number of nodes necessary (assuming we get all cpus on each node)
    nodes = np.ceil(ncpus/cluster_params['cpus_per_node']) 
    if specify_nodes:
        node_reqest_cmd = '#SBATCH --nodes=%i   #total number of nodes requested\n'%(nodes)
    else:
        # otherwise, let the scheduler be opportunistic
        node_reqest_cmd = '##SBATCH --nodes=%i   #total number of nodes requested\n'%(nodes)

    ### Set variables for PBS
    cmd_list = ['#!/bin/bash \n', '#SBATCH -J %s # job name \n' %run_config['run_name'],
                  '#SBATCH -o output_%j.txt # output and error file name (%j expands to jobID)\n',
                  '#SBATCH -n %i   #total number of mpi tasks requested\n' %ncpus,
                  node_reqest_cmd,
                  '#SBATCH --mem-per-cpu=%sG\n' %mem_GB,
                  '#SBATCH -p %s\n' %part,
                  '#SBATCH -t %s # run time (hh:mm:ss)\n'%walltime_str,
                  '#SBATCH --mail-user=%s\n'%email,
                  '#SBATCH --mail-type=begin  # email me when the job starts\n',
                  '#SBATCH --mail-type=end  # email me when the job finishes\n',
                  '\n',
                  'module load netcdf\n',
                  'module load openmpi/%s\n' %cluster_params['open_mpi_ver'],
                  '#module load openmpi\n',
                   run_cmd]


    ### Open output script and write header text
    fpath = os.path.join(results_fpath, 'run_mitgcm')
    
    with open(fpath, 'w') as f:
        text_str = "".join(cmd_list)
        f.writelines(text_str)
        
    print(text_str)
    
    
    
#################################################################
def getLastPickupIter(results_path):

    
#     if len(results_path)==0:
#         results_path = os.path.join(exp_dir, exp_name, 'results')
    
    # get most recent pickup file
    results_files = os.listdir(results_path)
    r = re.compile(r"^(pickup)\.\d+\.data")        
    pickup_files = sorted(list(filter(r.search, results_files)))

    # get latest iter number
    r2 = re.compile(r"\d+")
    
    try:
        iterN_str = r2.search(pickup_files[-1]).group()
    except IndexError:
        print(pickup_files)
        print("Warning: No pickup file found.")
        return 0
    
    return iterN_str

         
        
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
    
    dir_list = os.listdir(fname_dir)
    #print(dir_list)
    ver = 1
    
    new_fname = fname_root + '_old_v%03d'%ver + ext
    #print(new_fname)
    while new_fname in dir_list:
        ver+=1
        new_fname = fname_root + '_old_v%03d'%ver + ext
  
    #debug_here()
    new_fpath = os.path.join(fname_dir, new_fname)
    shutil.copyfile(fpath, new_fpath) 
    
    
def parseDataDiagnosticFile(params):
    
    """
    function to scrap data.diagnostics and extract dump frequencies
    """
    
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

    return diagnostic_dict
    
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
            
            
    
def extendRun(exp_path, newEndYr, ncpus, dump_config={}, new_diags={}, overWriteFiles=False, cleanup=True, 
              printNewRunScript=False, newRunName='', results_path='', overwrite_diag=False,
              wall_mins_per_sim_yr=60, max_walltime_hrs=48):
    
    
    secsInYr = 3600*24*365
    

    if overWriteFiles==False:
        print('------------------------------------------------------------------------------------')
        print("Warning: This is a dry run. New params and data are copied to temporary files.")
        print("To overwrite param and data files. set overWriteFiles=True")
        print('------------------------------------------------------------------------------------')
        
        
    if len(dump_config)>0 or len(new_diags)>0:
        warnings.warn("Code to update dump frequencies is still under construction.\n Ignoring dump_config and new_diags arguments")

    newEndTime = newEndYr*secsInYr 
    
    exp_name = exp_path.split('/')[-1]

    
    params,_ = raf.get_exp_params_py(exp_path)
    
    if len(newRunName)==0:
        newRunName = exp_name
    
        
    # backup old datafiles and save updated ones.
    backup_file(os.path.join(params['input_path'], 'data'), cleanup=cleanup)
    backup_file(os.path.join(params['input_path'], 'params.p'), cleanup=cleanup)
    backup_file(os.path.join(params['input_path'], 'data.diagnostics'), cleanup=cleanup)
    
    
    # get latest pickup file (TODO: allow restart from any previous pickup)
    results_path = params['results_path']
    lastIter = getLastPickupIter(results_path) # this method checks the last pickup file that was created
    newStartIter = int(lastIter)
    oldEndTime = int(newStartIter)*params['deltaT']
    
    #----------------update data file----------------------#
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
                        
                        
    
    #----------------update sbatch file----------------------#
    # update wall clock time
    backup_file(os.path.join(results_path, 'run_mitgcm'), cleanup=cleanup)
    pickupYear = int(lastIter)*params['deltaT']/(3600*24*365)
    nyrs = newEndYr - pickupYear
    walltime_hrs = wall_mins_per_sim_yr*nyrs/60
    
    # pad walltime by 5%
    walltime_hrs = np.ceil(walltime_hrs*1.05)
    print("Estimated runtime: %.1f hrs" %walltime_hrs)
    
    walltime_hrs_str = getWalltimeStr(walltime_hrs, max_walltime_hrs=max_walltime_hrs)

    oldLines = ['#SBATCH -t', 'cd /central', '#SBATCH -J'] # lines to be updated (only need starting string)
    newLines = ['#SBATCH -t %s  # run time (hh:mm:ss)\n' %walltime_hrs_str, '#cd %s\n' %results_path, 
               '#SBATCH -J %s\n' %newRunName] 
 
    runscript_fpath = os.path.join(params['results_path'], 'run_mitgcm')
    
    updateLineInText(runscript_fpath, oldLines, newLines, updateFileName=overWriteFiles, printNewFile=printNewRunScript)
    
       
    
    if overWriteFiles:
        # overwrites data and param files
        
        print("Ready to submit...")
        # update data, data.diagnostics, and param files 
        os.rename(os.path.join(params['input_path'], 'data_tmp'), os.path.join(params['input_path'], 'data'))
        
        
        #----------------update params file----------------------#
#         if update_params_file:
        with open(os.path.join(params['input_path'], 'params_tmp.p'), 'wb') as f_new:  
            with open(os.path.join(params['input_path'], 'params.p'), 'rb') as f_old:
                new_params = copy.deepcopy(params)
                new_params['mitgcm_params']['data'][2]['nIter0'] = newStartIter
                new_params['mitgcm_params']['data'][2]['endTime'] = newEndTime

                pickle.dump(new_params, f_new)
        
        
        
        #os.rename(os.path.join(params['input_path'], 'params_tmp.p'), os.path.join(params['input_path'], 'params.p'))
        if len(dump_config)==2:
            os.rename(os.path.join(params['input_path'], 'data.diagnostics.tmp'), 
                      os.path.join(params['input_path'], 'data.diagnostics'))
            
            
        print("Over-writing data and param files!")
            
def updateLineInText(fpath, oldLines, newLines, updateFileName=False, printNewFile=False):
    
    # TODO: double check this and update 
    
    fname = fpath.split('/')[-1]
    fname_dir = "/".join(fpath.split('/')[:-1])
    line_updated = False
    tmp_fpath = os.path.join(fname_dir, 'tmp_file.txt')
    oldLines_tup = tuple(oldLines)
    
    with open(tmp_fpath, 'w') as new_file:  
        with open(fpath, 'r') as old_file:
            for line in old_file:
                if line.startswith(oldLines_tup)==False:
                    new_file.write(line)
                else:
                    # find the specific parameter 
                    xx = np.array([line.startswith(pp) for pp in oldLines])
                    idx = np.arange(len(oldLines))[xx][0]
                    new_file.write(newLines[idx])
                    line_updated=True
                    
#     if line_updated:
#         print("%s updated!" %fname)
#     else:
#         print("WARNING: no lines in %s were updated" %fname)
      
    
    if printNewFile:
        with open(tmp_fpath, 'r') as f:
            print('---------------------------------')
            print(f.read())
        
    if updateFileName:
        os.rename(tmp_fpath, fpath)

    return


def getSavedIters(vname, results_path):
    

    results_files = os.listdir(results_path)
    vble_files = [f for f in results_files if f.startswith(vname) and f.endswith('.data')]
    
    if len(vble_files)==0:
        warnings.warn("No data files found!!")
    
    vble_files.sort()
    time_stmps = np.array([os.path.getmtime(os.path.join(results_path, f)) for f in vble_files])
    #time_stmps = time_stmps-time_stmps[0]

    # get iters 
    iters = [int(f.split('.')[1]) for f in vble_files]

    # get model time 
    iters = np.array([int(f.split('.')[1]) for f in vble_files])
    
    return iters


def rdmdsWrapper(fpath,dumpIter):
    
    """
    Convenience wrapper for rdmds.m that tries to read the specified file
    'fname' and iteration 'dumpIter'. If that name/iteration can't be
    found, we try dumpIter+1 and dumpIter-1 in case of rounding errors
    between Matlab and MITgcm.
    """
    
    A = MITgcmutils.mds.rdmds(fpath, dumpIter)        
    if len(A)==0:
        A = MITgcmutils.mds.rdmds(fpath, dumpIter-1)
    if len(A)==0:
        A = MITgcmutils.mds.rdmds(fpath, dumpIter+1)
        
    return A

def readIters(vname, tr, params, mask_zeros=True, returnTimeAvg=True, printStatus=True):
    
    if printStatus:
        print("Loading %s for years %s-%s..." %(vname, tr[0], tr[1]))
    
    dumpIters = params['dumpIters']
    if vname.startswith("La"):
        zlen = len(params['layers_bounds'])-1
    else:
        zlen = params['Nr']
    
    #print(dumpIters)
    #vble_t = np.zeros((params['Nx'], params['Ny'], zlen))
    #avg = np.zeros((params['Nx'], params['Ny'], zlen))
    #navg = 0
    secsInYear = 86400*365
    dt_sec = params['deltaT']
    ii = 0
    tvec = []
    tyears = params['dumpIters']*dt_sec/secsInYear
    
    # loop through output iterations
    for nn in range(len(dumpIters)):
        
        dump_time = params['dumpIters'][nn]*dt_sec
        dump_time_in_range = dump_time>=(tr[0]*secsInYear-2*dt_sec) and dump_time<=(tr[-1]*secsInYear+2*dt_sec)
        
        if dump_time_in_range:
            # Read the next iteration and check the data were found  
            try:
                vble_tmp = rdmdsWrapper(os.path.join(params['results_path'], vname), dumpIters[nn])
            except OSError:
                #break # continue?
                continue
                
            if len(vble_tmp)==0:
                print('Ran out of data at ITER= %s/%s t= %s days.' %(dumpIters[nn], params['nDumps'], tyears))
                break
            else:
                vble_tmp = vble_tmp.T
                if ii==0:
                    vble_t = vble_tmp[..., np.newaxis]
                else:
                    vble_t = np.concatenate([vble_t, vble_tmp[..., np.newaxis]], axis=3)
                    
                #avg += avg_tmp.T
                ii += 1
                tvec.append(dump_time)
                
    if ii==0:
        debug_here()
        raise ValueError("No output files found.")

    vble_t = np.ma.masked_invalid(vble_t)
    if mask_zeros:
        vble_t = np.ma.masked_where(vble_t==0, vble_t)

    
    # calculate average
    if returnTimeAvg:
        return vble_t.mean(axis=3)
    else:
        return vble_t, np.array(tvec)
    
    


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
       
    
                         