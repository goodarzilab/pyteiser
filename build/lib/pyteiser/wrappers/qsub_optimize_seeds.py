import argparse
import os
import math

from .. import IO

def handler():
    parser = argparse.ArgumentParser()

    # script to submit itself
    parser.add_argument("--script_to_sumbit", help="", type=str)
    parser.add_argument("--do_print_command", help="if the script should print the command it's submitting", type=str)

    # a single input file for optimize_seeds script to determine the number of chunks
    # and the size of chunk
    parser.add_argument("--unique_seeds_filename", help="best representatives of each family", type=str)
    parser.add_argument("--size_of_chunks", help="how many seeds should 1 process take on", type=float)

    # qsub submission parameters
    parser.add_argument("--python_binary", help="S parameter of qsub", type=str)
    parser.add_argument("--time_l_keyword", help="how to define time requirements in the -l parameters list", type=str)
    parser.add_argument("--time_required", help="-l h_rt parameter of qsub", type=str)
    parser.add_argument("--include_mem_free_parameter", help="should the mem_free parameter be included", type=str)
    parser.add_argument("--mem_free", help="-l mem_free parameter of qsub", type=str)
    parser.add_argument("--include_mem_scratch_parameter", help="should the scratch parameter be included", type=str)
    parser.add_argument("--mem_scratch", help="-l scratch parameter of qsub", type=str)
    parser.add_argument("--include_mem_parameter", help="should the mem parameter be included", type=str)
    parser.add_argument("--mem", help="-l mem parameter of qsub", type=str)
    parser.add_argument("--stderr_file", help="-e parameter of qsub", type=str)
    parser.add_argument("--stdout_file", help="-o parameter of qsub", type=str)
    parser.add_argument("--include_queue_parameter", help="should the desired queue be passed as a parameter", type=str)
    parser.add_argument("--queue", help="-q parameter of qsub", type=str)
    parser.add_argument("--restart", help="-r parameter of qsub", type=str)

# script-specific arguments go into the -ac argument

    parser.set_defaults(
        python_binary='/wynton/home/goodarzi/khorms/miniconda3/bin/python',
        time_l_keyword='h_rt',
        time_required='03:00:00',
        include_mem_free_parameter='y',
        mem_free='1G',
        include_mem_scratch_parameter='y',
        mem_scratch='1G',
        include_mem_parameter='n',
        mem='1G',
        stderr_file='/wynton/scratch/khorms/logs/test_stderr.txt',
        stdout_file='/wynton/scratch/khorms/logs/test_stdout.txt',
        include_queue_parameter='y',
        queue='long.q',
        restart='yes',
        do_print_command='no',

        size_of_chunks=10,
    )

    args, unknown = parser.parse_known_args()

    return args, unknown


def determine_number_of_chunks(args):
    seeds_initial = IO.read_motif_file(args.unique_seeds_filename)
    seeds_number = len(seeds_initial)
    number_of_chunks = math.ceil(seeds_number / args.size_of_chunks)
    return number_of_chunks



def unknown_args_to_string(unknown_args):
    return " ".join(unknown_args)


def construct_command(args, unknown_args, number_of_tasks):
    # syntax
    # qsub [ options ] [ command | -- [ command_args ]]
    # create a template
    command_template = ''
    # add qsub options
    command_template += "qsub -S {} -l {}={} ".format(args.python_binary, args.time_l_keyword,
                                                                     args.time_required)

    # some servers do not ask for scratch parameter
    if args.include_mem_free_parameter == 'y' or args.include_mem_free_parameter == 'yes':
        command_template += "-l mem_free={} ".format(args.mem_free)
    if args.include_mem_scratch_parameter == 'y' or args.include_mem_scratch_parameter == 'yes':
        command_template += "-l scratch={} ".format(args.mem_scratch)
    if args.include_mem_parameter == 'y' or args.include_mem_parameter == 'yes':
        command_template += "-l mem={} ".format(args.mem)
    if args.include_queue_parameter == 'y' or args.include_queue_parameter == 'yes':
        command_template += "-q {} ".format(args.queue)

    command_template += "-e {} -o {} ".format(args.stderr_file, args.stdout_file)

    command_template += "-r {} ".format(args.restart)

    # add the array jobs option
    command_template += "-t {}-{} ".format(1, number_of_tasks)

    # add the script to run
    command_template += "{} ".format(args.script_to_sumbit)

    # add unique_seeds_filename parameter
    command_template += "--unique_seeds_filename {} ".format(args.unique_seeds_filename)
    # add the parameter defining the size of each chunk
    command_template += "--size_of_chunks {} ".format(args.size_of_chunks)

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

    number_of_tasks = determine_number_of_chunks(args)

    # assemble qsub command
    command = construct_command(args, unknown_args, number_of_tasks)

    # get the directory where the script to be submitted is located
    wording_dir = get_submission_working_directory(args.script_to_sumbit)

    # submit the actual job
    submit_job(wording_dir, command, args.do_print_command)


if __name__ == "__main__":
    main()