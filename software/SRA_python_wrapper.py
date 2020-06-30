#!/usr/bin/env python
# coding: utf-8

# In[12]:


import os
import subprocess
from multiprocessing import Pool
import warnings

print('see CD8T_RA download_dataset.py for an example of usage of SRA wrapper functions')

def download_sra(download_path, SRAaccession_name):
    '''Wrapper to download SRA data. Download_path is path to directory that is being downloaded to (not included /sra/)
    Requirements: 1) SRA toolkit must point to download_path 
    (see https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std or Google Docs Software notes)
    2) SRA_accession_name is a path/to/SRA_accession_list.txt'''
    
    warnings.warn("This function has not been fully tested")
    
    print('download SRA files')
    try:
        prf = 'prefetch --option-file ' + SRAaccession_name
        output = subprocess.check_output(prf, stderr=subprocess.STDOUT, universal_newlines=True, shell=True)
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)
        error_ = 'Did you change the sra toolkit path to the download_path? '
        error_ += 'See https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std'
        print(error_)

def sra_to_fastq(sra_file, local_data_path, sra_path, n_files, kwargs = {}):
    '''SRA file is the file name without a .sra extension. Output path is directory in which to start .fastq output
    kwargs is a dictionary of flags and commands'''
    print('File # {}: '.format(n_files) + sra_file)
    
    # generate fastq files
    if os.path.isfile(sra_path + 'sra/' + sra_file + '.sra'): # only spend seven minutes if the file exists
        try: 
            fqd = 'fastq-dump ' + sra_file + ' '
            for flag, option in kwargs.items(): # add flag options
                if option != '':
                    fqd += '--' + flag + ' ' + option + ' '
                else:
                    fqd += '--' + flag + ' '
            
            fqd = fqd[:-1] # get rid of space
            output = subprocess.check_output(fqd, stderr=subprocess.STDOUT, shell=True, 
                                             universal_newlines=True)
        except:
            str_ = sra_file + ' failed to generate fastq'
            with open(local_data_path + 'interim/SRA_to_fastq_errors.txt', 'a') as f:
                f.write(str_ + '\n')
            pass 
    else:
        str_ = sra_file + '.sra does not exist '
        with open(local_data_path + 'interim/SRA_to_fastq_errors.txt', 'a') as f:
            f.write(str_ + '\n')
        pass 

def convert_sra(sra_path, local_data_path, cr = 10, kwargs = {}):
    '''Wrapper to convert sra to fastq. Download path is path/to/sra_directory. 
    Will parallelize sra conversion to cr cores.
    sra_path is a path/to/sra_directory/parent/, 
    kwargs in the format {'flag1': 'option1', 'flag2': 'option2'}
    specify output path with kwargs = {'outdir': path/to/output_directory/} 
    This function has been tested and will write all errors to interim/SRA_to_fastq_errors.txt
    kwargs is any additional flags to add. '''
    # generate a file for any errors that occur
    with open(local_data_path + 'interim/SRA_to_fastq_errors.txt', 'w+') as f: 
        pass 
    
    print('Convert to fastq format')
    sra_files = os.listdir(sra_path + 'sra') # get downloaded sra files
#     # TEST (normal, file doesn't exist, one-sided from split-file)
#     sra_files = ['SRR4785592.sra', 'SRR2NKBOO.sra', 'SRR4785803.sra'] 
    try:
        sra_files = [file.split('.sra')[0] for file in sra_files] # remove extension
    except ValueError:
        print(sra_path + ' contains files that do not have a .sra extension')
        
    if 'outdir' in kwargs.keys(): #generate output directory if it does not already exist
        output_path = kwargs['outdir']
        if not os.path.isdir(output_path):
            print(output_path + ' does not exist, generating directory')
            os.mkdir(output_path)
        else:
            if len(os.listdir(output_path)) > 0:
                raise ValueError('Specified output directory must be empty')
    
    # parallel job
    i_2, i_3, i_4 = len(sra_files)*[local_data_path], len(sra_files)*[sra_path], range(1, len(sra_files) + 1)
    i_5 = len(sra_files)*[kwargs] 
    pool = Pool(processes=cr) 
    p = pool.starmap(sra_to_fastq, list(zip(sra_files, i_2, i_3, i_4, i_5)))
    pool.close()

    




# In[2]:


# def sra_to_fastq(sra_file, output_path, local_data_path, sra_path, kwargs):
#     '''SRA file is the file name without a .sra extension. Output path is directory in which to start .fastq output
#     kwargs is a dictionary of flags and commands'''
#     print(sra_file)
    
#     # generate fastq files
#     if os.path.isfile(sra_path + 'sra/' + sra_file + '.sra'): # only spend seven minutes if the file exists
#         try: 
#             fqd = 'fastq-dump '  #do conversion from sra to fastq
#             for flag, option in kwargs.items():
#                 if option != '':
#                     fqd += '--' + flag + ' ' + option + ' '
#                 else:
#                     fqd += '--' + flag + ' '
#             fqd += sra_file
#             output = subprocess.check_output(fqd, stderr=subprocess.STDOUT, shell=True, 
#                                              universal_newlines=True)
#             # time out 420: shouldn't take more than seven minutes to run 
#         except:
#             str_ = sra_file + ' failed to generate fastq'
#             with open(local_data_path + 'interim/SRA_to_fastq_errors.txt', 'a') as f:
#                 f.write(str_ + '\n')
#             pass # no files to move if this part failed
        
#         else: # move files if fastq part worked
#             try: # move read1 file to output path
#                 cmd = 'scp ' + sra_file + '_1.fastq ' + output_path
#                 output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, 
#                                                  universal_newlines=True)
#             except: 
#                 str_ = sra_file + ' does not generate a read 1 file'
#                 with open(local_data_path + 'interim/SRA_to_fastq_errors.txt', 'a') as f:
#                     f.write(str_ + '\n')
#             else:
#                 cmd = 'rm ' + sra_file + '_1.fastq '
#                 output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, 
#                                                  universal_newlines=True)

#             try: # move read2 file to output path
#                 cmd = 'scp ' + sra_file + '_2.fastq ' + output_path
#                 output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, 
#                                                  universal_newlines=True)
#             except: 
#                 str_ = sra_file + ' does not generate a read 2 file'
#                 with open(local_data_path + 'interim/SRA_to_fastq_errors.txt', 'a') as f:
#                     f.write(str_ + '\n')
#             else:
#                 cmd = 'rm ' + sra_file + '_2.fastq '
#                 output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, 
#                                                  universal_newlines=True)
#     else:
#         str_ = sra_file + '.sra does not exist '
#         with open(local_data_path + 'interim/SRA_to_fastq_errors.txt', 'a') as f:
#             f.write(str_ + '\n')
#         pass # no files to move if this part failed

