##
# \file run_reconstruction_pipeline.py
# \brief      Script to execute entire reconstruction pipeline
#
# \author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
# \date       October 2017
#

import numpy as np
import os
import re

import niftymic.validation.simulate_stacks_from_reconstruction as \
    simulate_stacks_from_reconstruction
import niftymic.validation.evaluate_simulated_stack_similarity as \
    evaluate_simulated_stack_similarity
import niftymic.validation.show_evaluated_simulated_stack_similarity as \
    show_evaluated_simulated_stack_similarity
import niftymic.validation.export_side_by_side_simulated_vs_original_slice_comparison as \
    export_side_by_side_simulated_vs_original_slice_comparison
import pysitk.python_helper as ph
from niftymic.definitions import DIR_TEMPLATES
from niftymic.utilities.input_arparser import InputArgparser


def main():
    pid = "49"
    method     = "manual" #auto, manual
    input_dir  = "{0:}/input".format(pid)
    mask_dir   = "{0:}/mask_{1:}/mask".format(pid, method)
    output_dir = "{0:}/input_preprocess".format(pid)

    # get input stack names
    files = os.listdir(input_dir)
    input_files = []
    mask_files  = []
    for file in files:
        if (("nii.gz" in file)):
            input_files.append("{0:}/{1:}".format(input_dir, file))
            mask_name = "{0:}/{1:}".format(mask_dir, file)
            assert(os.path.isfile(mask_name))
            mask_files.append(mask_name)

    cmd_args = []
    cmd_args.append("--filenames %s" % (" ").join(input_files))
    cmd_args.append("--filenames-masks %s" % (" ").join(mask_files))
    cmd_args.append("--dir-output %s" % output_dir)
    cmd_args.append("--prefix-output ''")
    cmd = "niftymic_correct_bias_field %s" % (" ").join(cmd_args)
    time_start_bias = ph.start_timing()
    exit_code = ph.execute_command(cmd)
    elapsed_time_bias = ph.stop_timing(time_start_bias)
    print("Computational Time for Bias Field Correction: %s" %
          elapsed_time_bias)
    if exit_code != 0:
        raise RuntimeError("Bias field correction failed")
    return 0


if __name__ == '__main__':
    main()
