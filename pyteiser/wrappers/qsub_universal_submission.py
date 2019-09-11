import argparse
import random
import os


# The overall submission scheme:
# This (qsub_universal_submission.py) is the master submission script
# You call it with the arguments, specifying:
# 1) the Python script you want to submit
# 2) list of input files you want to run it on (more specifically, input/output folders, prefixes and list of file numbers)
# 3) list of submission arguments, like time and memory constraints
# 4) list of arguments that need to be passed to Python script you want to run
#
# This script then creates a qsub command of a standard form of
# qsub script_to_call [qsub options] [array of jobs t] [script parameters]
# and then submits it.
# Each individual script that is called through qsub then has to:
# 1) move to the working directory of the original script and one level up to import the rest of pyteiser
# 2) read the environmental variables to find out what task number was assigned to it
# 3) read the task number-to-input file mapping file to find out which file should it work on
# 4) read its own parameters
# 5) run
#
# Theoretically, I could submit each individual script with all the necessary qsub parameters avoiding such
# complicated scheme with an external wrapper
# The reason for doing things is that stupid SGE doesn't allow listing the indices of array jobs with comma
# (even though the documentation says it does). Therefore, if I have, say, files number 13,27,57 and I want to run
# my script on them, I have to either submit 3 jobs (bad) or create a separate file with mapping
# To avoid creating mapping files for each script individually, I use a system where a wrapper script does
# it automatically given the list of files that need to be processed


# maybe http://45.76.113.195/?a/44467174 ?

def handler():
    parser = argparse.ArgumentParser()

    #script to submit itself
    parser.add_argument("--script_to_sumbit", help="", type=str)
    parser.add_argument("--do_print_command", help="if the script should print the command it's submitting", type=str)

    # list of input file indices (to create mapping to task numbers)
    parser.add_argument("--input_indices_list_file", help="input: list of indices of files to process", type=str)
    parser.add_argument("--mapping_task_ids_folder", help="output: mapping of file indices to task ids", type=str)

    # qsub submission parameters
    parser.add_argument("--python_binary", help="S parameter of qsub", type=str)
    parser.add_argument("--time_required", help="-l h_rt parameter of qsub", type=str)
    parser.add_argument("--mem_free", help="-l mem_free parameter of qsub", type=str)
    parser.add_argument("--include_mem_scratch_parameter", help="should the scratch parameter be included", type=str)
    parser.add_argument("--mem_scratch", help="-l scratch parameter of qsub", type=str)
    parser.add_argument("--stderr_file", help="-e parameter of qsub", type=str)
    parser.add_argument("--stdout_file", help="-o parameter of qsub", type=str)
    parser.add_argument("--include_queue_parameter", help="should the desired queue be passed as a parameter", type=str)
    parser.add_argument("--queue", help="-q parameter of qsub", type=str)
    parser.add_argument("--restart", help="-r parameter of qsub", type=str)

# script-specific arguments go into the -ac argument

    parser.set_defaults(
        python_binary='/wynton/home/goodarzi/khorms/miniconda3/bin/python',
        time_required='50:00:00',
        mem_free='1G',
        include_mem_scratch_parameter='y',
        mem_scratch='1G',
        stderr_file='/wynton/scratch/khorms/logs/test_stderr.txt',
        stdout_file='/wynton/scratch/khorms/logs/test_stdout.txt',
        include_queue_parameter='y',
        queue='long.q',
        restart='yes',
        do_print_command='no'
    )

    args, unknown = parser.parse_known_args()

    return args, unknown


def create_mapping(args):
    with open(args.input_indices_list_file, 'r') as rf:
        full_string = rf.read()
        full_string = full_string.rstrip()
    indices_list_str = full_string.split(', ')
    indices_list = [int(x) for x in indices_list_str]

    random_int_list = [str(random.randint(0, 9)) for p in range(7)]
    random_int_string = "".join(random_int_list)
    unique_masking_file = "%s/mapping_%s.txt" % (args.mapping_task_ids_folder, random_int_string)
    with open(unique_masking_file, 'w') as wf:
        for task_id, index in enumerate(indices_list):
            current_string = "%d\t%d\n" % (task_id+1, index)
            wf.write(current_string)
    number_of_tasks = len(indices_list)

    return unique_masking_file, number_of_tasks



def unknown_args_to_string(unknown_args):
    return " ".join(unknown_args)


def construct_command(unique_masking_file, number_of_tasks,
                      args, unknown_args):
    # syntax
    # qsub [ options ] [ command | -- [ command_args ]]
    # create a template
    command_template = ''
    # add qsub options
    command_template += "qsub -S {} -l h_rt={} -l mem_free={} ".format(args.python_binary, args.time_required, args.mem_free)

    # some servers do not ask for scratch parameter
    if args.include_mem_scratch_parameter == 'y' or args.include_mem_scratch_parameter == 'yes':
        command_template += "-l scratch={} ".format(args.mem_scratch)

    command_template += "-e {} -o {} ".format(args.stderr_file, args.stdout_file)

    if args.include_queue_parameter == 'y' or args.include_queue_parameter == 'yes':
        command_template += "-q {} ".format(args.queue)

    command_template += "-r {} ".format(args.restart)

    # add the array jobs option
    command_template += "-t {}-{} ".format(1, number_of_tasks)

    # add the script to run
    command_template += "{} ".format(args.script_to_sumbit)

    # add the parameter defining the mapping task ids file
    command_template += "--task_mapping_file {} ".format(unique_masking_file)

    # add the parameters for the script itself
    unknown_args_string = unknown_args_to_string(unknown_args)
    command_template += "{} ".format(unknown_args_string)

    return command_template


def get_submission_working_directory(script_to_sumbit):
    return os.path.dirname(script_to_sumbit)


def submit_job(wording_dir, command, do_print):
    if do_print == 'y' or do_print == 'yes':
        print("The command is: ")
        print(command)
    os.chdir(wording_dir)
    os.system(command)


def main():
    args, unknown_args = handler()

    # use the list of input files indices to consrtuct a mapping from task id to input file index
    unique_masking_file, number_of_tasks = create_mapping(args)

    # assemble qsub command
    command = construct_command(unique_masking_file, number_of_tasks,
                                args, unknown_args)

    # get the directory where the script to be submitted is located
    wording_dir = get_submission_working_directory(args.script_to_sumbit)

    # submit the actual job
    submit_job(wording_dir, command, args.do_print_command)


if __name__ == "__main__":
    main()