import numba
import time
import os
import sys
import subprocess

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)


import glob_var
import structures


def parse_task_mapping_file(task_mapping_file):
    mapping_dict = {}
    with open(task_mapping_file, 'r') as rf:
        for line in rf:
            stripped_line = line.rstrip()
            splitted_line = stripped_line.split('\t')
            task_id = splitted_line[0]
            file_index = splitted_line[1]
            mapping_dict[task_id] = file_index
    return mapping_dict


def get_env_variables():
    working_home_dir = os.getcwd()
    local_scratch = os.getenv('TMPDIR')
    task_id = os.environ["SGE_TASK_ID"]
    job_id = os.environ["JOB_ID"]

    env_variables_dict = {"working_dir": working_home_dir,
                          "local_scratch": local_scratch,
                          "task_id": task_id,
                          "job_id": job_id}
    return env_variables_dict


def print_qstat_proc(env_variables_dict, path_to_qstat):
    #subprocess.call(qstat_command, shell=False)  # qstat is an executable itself, you don't need shell to run it
    # also subprocess doesn't want to find qstat if I provide the arguments as a string
    # for more detail, see either the comment of jfs here https://stackoverflow.com/questions/18962785/oserror-errno-2-no-such-file-or-directory-while-using-python-subprocess-in-dj
    # or this post https://stackoverflow.com/questions/4795190/passing-variables-to-a-subprocess-call

    subprocess.call([path_to_qstat, '-j', env_variables_dict["job_id"]], shell=False)  # qstat is an executable itself, you don't need shell to run it


