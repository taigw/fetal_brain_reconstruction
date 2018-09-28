"""
\file run_reconstruction_pipeline.py
\brief      Script to execute entire reconstruction pipeline
\author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
\date       December 2017
"""

import os
import sys
import re
import numpy as np

import pysitk.python_helper as ph
from niftymic.definitions import DIR_TEMPLATES
import niftymic.utilities.template_stack_estimator as tse
import niftymic.validation.evaluate_image_similarity as evaluate_image_similarity
import niftymic.validation.evaluate_slice_residual_similarity as esrs
import niftymic.application.multiply_stack_with_mask as mswm

import src.data_reader as dr
import src.utilities as utils
import src.analyse_image_similarities
import src.analyse_slice_residual_similarities
import src.statistical_evaluation as se

import pysitk.simple_itk_helper as sitkh
import pysitk.data_anonymizer as da

INDIVIDUAL_CASE_IDS = [

    "05_UZL11-Study16",
    "06_UZL11-Study1",  # T2_Brain_1101/T2_spine_1901 can be opened; nice recon

    "10_UZL8-Study13",
    "11_UZL8-Study1",

    "14_UZL10-Study15",
    "15_UZL10-Study1",  # 5/5
    "16_UZL00056-Study1",  # template-alignment incorrect
    "17_UZL00056-Study2",
    "18_UZL7-Study12",
    "19_UZL7-Study1",  # is there a coronal stack?
    "20_UZL00057-Study1",
    "21_UZL00057-Study2",

    "25_UZL00059-Study1",  # 3/5
    "26_UZL00059-Study2",

    "34_UZL00065-Study1",
    "35_UZL00065-Study2",  # 3/5
    "36_UZL9-Study14",
    "37_UZL9-Study1",
    "38_UZL00066-Study1",
    "39_UZL00066-Study2",
    "43_UZL62-Study3",
    "44_UZL62-Study1",
    # "45_UZL62-Study2",  # FU scan  # template-alignment incorrect

    # "49_UZL4-Study1",  # DISCARD (BTFE sequence)
    "50_UZL4-Study8",  # recon very grainy - discard associated stack?
    "51_UZL00072-Study1",  # 2/5; but detection much worse
    "52_UZL00072-Study2",
    "53_UZL17-Study1",  # 3/5
    "54_UZL17-Study2",  # 3/5
    "55_UZL17-Study21",  # 5/5

    "57_UZL00074-Study1",
    "58_UZL00074-Study2",
    "66_UZL00082-Study1",  # poor recon # template-alignment incorrect
    "67_UZL00082-Study2",  # not much information; 2 stacks cropped; poor recon
]

RECON_FAILED_CASE_IDS = [
    "05_UZL11-Study16",  # 1/5; difficult one; almost no image valid; Cor image? # Pre
    "11_UZL8-Study1",  # 1/5 + template alignment incorrect # Post
    "16_UZL00056-Study1",  # template alignment incorrect # Pre
    "52_UZL00072-Study2",  # template-alignment incorrect; recon 4/5 # Post
    "66_UZL00082-Study1",  # template-alignment incorrect; recon 3/5 # Pre
    "67_UZL00082-Study2",  # recon 2/5, not much information; 2 stacks cropped # Pre
]


SEG_MODES = [
    "seg_manual",
    "seg_auto",  # used to register to template
    "detect",
]
RECON_MODES = [
    "seg_manual",
    "seg_auto",
    "detect",
    "irtk",
]


def main():

    time_start = ph.start_timing()

    flag_individual_cases_only = 1

    flag_batch_script = 0
    batch_ctr = [32]

    flag_correct_bias_field = 0
    # flag_correct_intensities = 0

    flag_collect_segmentations = 0
    flag_select_images_segmentations = 0

    flag_reconstruct_volume_subject_space = 0
    flag_reconstruct_volume_subject_space_irtk = 0
    flag_reconstruct_volume_subject_space_show_comparison = 0
    flag_register_to_template = 0
    flag_register_to_template_irtk = 0
    flag_show_srr_template_space = 0
    flag_reconstruct_volume_template_space = 0
    flag_collect_volumetric_reconstruction_results = 0
    flag_show_volumetric_reconstruction_results = 0

    flag_rsync_stuff = 0

    # Analysis
    flag_remove_failed_cases_for_analysis = 1
    flag_postop = 2  # 0... preop, 1...postop, 2... pre+postop

    flag_evaluate_image_similarities = 0
    flag_analyse_image_similarities = 1

    flag_evaluate_slice_residual_similarities = 0
    flag_analyse_slice_residual_similarities = 0

    flag_analyse_stacks = 0
    flag_analyse_qualitative_assessment = 0

    flag_collect_data_blinded_analysis = 0
    flag_anonymize_data_blinded_analysis = 0

    provide_comparison = 0
    intensity_correction = 1
    isotropic_resolution = 0.75
    alpha = 0.02
    outlier_rejection = 1
    threshold = 0.7
    threshold_first = 0.6

    # metric = "ANTSNeighborhoodCorrelation"
    # metric_radius = 5
    # multiresolution = 0

    prefix_srr = "srr_"
    prefix_srr_qa = "masked_"

    # ----------------------------------Set Up---------------------------------
    if flag_correct_bias_field:
        dir_batch = os.path.join(utils.DIR_BATCH_ROOT, "BiasFieldCorrection")
    elif flag_reconstruct_volume_subject_space:
        dir_batch = os.path.join(utils.DIR_BATCH_ROOT,
                                 "VolumetricReconstructionSubjectSpace")
    elif flag_register_to_template:
        dir_batch = os.path.join(utils.DIR_BATCH_ROOT,
                                 "VolumetricReconstructionRegisterToTemplate")
    elif flag_reconstruct_volume_template_space:
        dir_batch = os.path.join(utils.DIR_BATCH_ROOT,
                                 "VolumetricReconstructionTemplateSpace")
    else:
        dir_batch = os.path.join(utils.DIR_BATCH_ROOT, "foo")
    file_prefix_batch = os.path.join(
        dir_batch,
        "command")

    if flag_batch_script:
        verbose = 0
    else:
        verbose = 1

    data_reader = dr.ExcelSheetDataReader(utils.EXCEL_FILE)
    data_reader.read_data()
    cases = data_reader.get_data()

    if flag_analyse_qualitative_assessment:
        data_reader = dr.ExcelSheetQualitativeAssessmentReader(utils.QA_FILE)
        data_reader.read_data()
        qualitative_assessment = data_reader.get_data()

        statistical_evaluation = se.StatisticalEvaluation(
            qualitative_assessment)
        statistical_evaluation.run_tests(ref="seg_manual")
        ph.exit()

    cases_similarities = []
    cases_stacks = []

    if flag_individual_cases_only:
        N_cases = len(INDIVIDUAL_CASE_IDS)
    else:
        N_cases = len(cases.keys())

    i_case = 0
    for case_id in sorted(cases.keys()):
        if flag_individual_cases_only and case_id not in INDIVIDUAL_CASE_IDS:
            continue
        if not flag_analyse_image_similarities and \
                not flag_analyse_slice_residual_similarities:
            i_case += 1
            ph.print_title("%d/%d: %s" % (i_case, N_cases, case_id))

        if flag_rsync_stuff:
            dir_output = utils.get_directory_case_recon_seg_mode(
                case_id=case_id,
                recon_space="template_space",
                seg_mode="")

            dir_input = re.sub(
                "Volumes/spina/",
                "Volumes/medic-volumetric_res/SpinaBifida/",
                dir_output)
            cmd = "rsync -avuhn --exclude 'motion_correction' %sseg_manual %s" % (
                dir_input, dir_output)
            ph.print_execution(cmd)
            # ph.execute_command(cmd)

        # -------------------------Correct Bias Field--------------------------
        if flag_correct_bias_field:
            filenames = utils.get_filenames_preprocessing_bias_field(
                case_id)
            paths_to_filenames = [os.path.join(
                utils.get_directory_case_original(case_id), f) for f in filenames]
            dir_output = utils.get_directory_case_preprocessing(
                case_id, stage="01_N4ITK")

            # no image found matching the pattern
            if len(paths_to_filenames) == 0:
                continue

            cmd_args = []
            cmd_args.append("--filenames %s" % " ".join(paths_to_filenames))
            cmd_args.append("--dir-output %s" % dir_output)
            cmd_args.append("--prefix-output ''")
            cmd = "niftymic_correct_bias_field %s" % (" ").join(cmd_args)

            ph.execute_command(
                cmd,
                flag_print_to_file=flag_batch_script,
                path_to_file="%s%d.txt" %
                (file_prefix_batch, ph.add_one(batch_ctr)))

        # # Skip case in case segmentations have not been provided yet
        # if not ph.directory_exists(utils.get_directory_case_segmentation(
        #         case_id, utils.SEGMENTATION_INIT, SEG_MODES[0])):
        #     continue

        # ------------------------Collect Segmentations------------------------
        if flag_collect_segmentations:
            # Skip case in case segmentations have been collected already
            if ph.directory_exists(utils.get_directory_case_segmentation(
                    case_id, utils.SEGMENTATION_SELECTED, SEG_MODES[0])):
                ph.print_info("skipped")
                continue

            filenames = utils.get_segmented_image_filenames(
                case_id,
                subfolder=utils.SEGMENTATION_INIT)

            for i_seg_mode, seg_mode in enumerate(SEG_MODES):
                directory_selected = utils.get_directory_case_segmentation(
                    case_id, utils.SEGMENTATION_SELECTED, seg_mode)
                ph.create_directory(directory_selected)
                paths_to_filenames_init = [
                    os.path.join(utils.get_directory_case_segmentation(
                        case_id, utils.SEGMENTATION_INIT, seg_mode), f)
                    for f in filenames]
                paths_to_filenames_selected = [
                    os.path.join(directory_selected, f)
                    for f in filenames]
                for i in range(len(filenames)):
                    cmd = "cp -p %s %s" % (
                        paths_to_filenames_init[i], paths_to_filenames_selected[i])
                    # ph.print_execution(cmd)
                    ph.execute_command(cmd)

        if flag_select_images_segmentations:
            filenames = utils.get_segmented_image_filenames(
                case_id,
                subfolder=utils.SEGMENTATION_SELECTED)
            paths_to_filenames = [os.path.join(
                utils.get_directory_case_preprocessing(
                    case_id, stage="01_N4ITK"), f) for f in filenames]
            paths_to_filenames_masks = [os.path.join(
                utils.get_directory_case_segmentation(
                    case_id, utils.SEGMENTATION_SELECTED, "seg_manual"), f)
                for f in filenames]
            for i in range(len(filenames)):
                ph.show_niftis(
                    [paths_to_filenames[i]],
                    segmentation=paths_to_filenames_masks[i],
                    # viewer="fsleyes",
                )
                ph.pause()
                ph.killall_itksnap()

        # # -------------------------Correct Intensities-----------------------
        # if flag_correct_intensities:
        #     filenames = utils.get_segmented_image_filenames(case_id)
        #     paths_to_filenames_bias = [os.path.join(
        #         utils.get_directory_case_preprocessing(
        #             case_id, stage="01_N4ITK"), f) for f in filenames]
        #     print paths_to_filenames_bias

        # -----------------Reconstruct Volume in Subject Space-----------------
        if flag_reconstruct_volume_subject_space:

            filenames = utils.get_segmented_image_filenames(
                case_id,
                subfolder=utils.SEGMENTATION_SELECTED)
            # filenames = filenames[0:2]

            paths_to_filenames = [os.path.join(
                utils.get_directory_case_preprocessing(
                    case_id, stage="01_N4ITK"), f) for f in filenames]

            # Estimate target stack
            target_stack_index = utils.get_target_stack_index(
                case_id, utils.SEGMENTATION_SELECTED, "seg_auto", filenames)

            for i, seg_mode in enumerate(SEG_MODES):
                # Get mask filenames
                paths_to_filenames_masks = [os.path.join(
                    utils.get_directory_case_segmentation(case_id, utils.SEGMENTATION_SELECTED, seg_mode), f)
                    for f in filenames]

                if flag_reconstruct_volume_subject_space_irtk:
                    if seg_mode != "seg_manual":
                        continue
                    utils.export_irtk_call_to_workstation(
                        case_id=case_id,
                        filenames=filenames,
                        seg_mode=seg_mode,
                        isotropic_resolution=isotropic_resolution,
                        target_stack_index=target_stack_index,
                        kernel_mask_dilation=(15, 15, 4))

                else:
                    dir_output = utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="subject_space",
                        seg_mode=seg_mode)
                    # dir_output = "/tmp/foo"

                    cmd_args = []
                    cmd_args.append("--filenames %s" %
                                    " ".join(paths_to_filenames))
                    cmd_args.append("--filenames-masks %s" %
                                    " ".join(paths_to_filenames_masks))
                    cmd_args.append("--dir-output %s" % dir_output)
                    cmd_args.append("--use-masks-srr 0")
                    cmd_args.append("--isotropic-resolution %f" %
                                    isotropic_resolution)
                    cmd_args.append("--target-stack-index %d" %
                                    target_stack_index)
                    cmd_args.append("--intensity-correction %d" %
                                    intensity_correction)
                    cmd_args.append("--outlier-rejection %d" %
                                    outlier_rejection)
                    cmd_args.append("--threshold-first %f" % threshold_first)
                    cmd_args.append("--threshold %f" % threshold)
                    # cmd_args.append("--metric %s" % metric)
                    # cmd_args.append("--multiresolution %d" % multiresolution)
                    # cmd_args.append("--metric-radius %s" % metric_radius)
                    # if i > 0:
                    #     cmd_args.append("--reconstruction-space %s" % (
                    #         utils.get_path_to_recon(
                    #             utils.get_directory_case_recon_seg_mode(
                    #                 case_id, "seg_manual"))))
                    # cmd_args.append("--two-step-cycles 0")
                    cmd_args.append("--verbose %d" % verbose)
                    cmd_args.append("--provide-comparison %d" %
                                    provide_comparison)
                    # cmd_args.append("--iter-max 1")

                    cmd = "niftymic_reconstruct_volume %s" % (
                        " ").join(cmd_args)

                    ph.execute_command(
                        cmd,
                        flag_print_to_file=flag_batch_script,
                        path_to_file="%s%d.txt" %
                        (file_prefix_batch, ph.add_one(batch_ctr)))

        if flag_reconstruct_volume_subject_space_show_comparison:
            recon_paths = []
            for seg_mode in SEG_MODES:
                path_to_recon = utils.get_path_to_recon(
                    utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="subject_space",
                        seg_mode=seg_mode))
                recon_paths.append(path_to_recon)
            recon_path_irtk = os.path.join(
                utils.get_directory_case_recon_seg_mode(
                    case_id=case_id,
                    recon_space="subject_space",
                    seg_mode="IRTK"),
                "IRTK_SRR.nii.gz")
            show_modes = list(SEG_MODES)
            if ph.file_exists(recon_path_irtk):
                recon_paths.append(recon_path_irtk)
                show_modes.append("irtk")
            ph.show_niftis(recon_paths)
            ph.print_info("Sequence: %s" % (" -- ").join(show_modes))
            ph.pause()
            ph.killall_itksnap()

        # -------------------------Register to template------------------------
        if flag_register_to_template:
            for seg_mode in SEG_MODES:

                cmd_args = []
                # register seg_auto-recon to template space
                if seg_mode == "seg_auto":

                    path_to_recon = utils.get_path_to_recon(
                        utils.get_directory_case_recon_seg_mode(
                            case_id=case_id,
                            recon_space="subject_space",
                            seg_mode=seg_mode))

                    template_stack_estimator = \
                        tse.TemplateStackEstimator.from_mask(
                            ph.append_to_filename(path_to_recon, "_mask"))
                    path_to_reference = \
                        template_stack_estimator.get_path_to_template()

                    dir_input_motion_correction = os.path.join(
                        utils.get_directory_case_recon_seg_mode(
                            case_id=case_id,
                            recon_space="subject_space",
                            seg_mode=seg_mode),
                        "motion_correction")

                    dir_output = utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="template_space",
                        seg_mode=seg_mode)
                    # dir_output = "/home/mebner/tmp"
                    # # ------- DELETE -----
                    # dir_output = re.sub("data", "foo+1", dir_output)
                    # dir_output = re.sub(
                    #     "volumetric_reconstruction/20180126/template_space/seg_auto",
                    #     "", dir_output)
                    # # -------
                    # cmd_args.append("--use-fixed-mask 1")
                    cmd_args.append("--use-moving-mask 1")

                    # HACK
                    path_to_initial_transform = os.path.join(
                        utils.DIR_INPUT_ROOT_DATA, case_id,
                        "volumetric_reconstruction", "20180126",
                        "template_space", "seg_manual",
                        "registration_transform_sitk.txt")
                    cmd_args.append("--initial-transform %s" %
                                    path_to_initial_transform)
                    cmd_args.append("--use-flirt 0")
                    cmd_args.append("--use-regaladin 1")
                    cmd_args.append("--test-ap-flip 0")

                # register remaining recons to registered seg_auto-recon
                else:
                    path_to_reference = utils.get_path_to_recon(
                        utils.get_directory_case_recon_seg_mode(
                            case_id=case_id,
                            recon_space="template_space",
                            seg_mode="seg_auto"),
                        suffix="ResamplingToTemplateSpace",
                    )
                    path_to_initial_transform = os.path.join(
                        utils.get_directory_case_recon_seg_mode(
                            case_id=case_id,
                            recon_space="template_space",
                            seg_mode="seg_auto"),
                        "registration_transform_sitk.txt")

                    path_to_recon = utils.get_path_to_recon(
                        utils.get_directory_case_recon_seg_mode(
                            case_id=case_id,
                            recon_space="subject_space",
                            seg_mode=seg_mode))
                    dir_input_motion_correction = os.path.join(
                        utils.get_directory_case_recon_seg_mode(
                            case_id=case_id,
                            recon_space="subject_space",
                            seg_mode=seg_mode),
                        "motion_correction")
                    dir_output = utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="template_space",
                        seg_mode=seg_mode)

                    cmd_args.append("--use-fixed-mask 0")
                    cmd_args.append("--use-moving-mask 0")
                    cmd_args.append("--initial-transform %s" %
                                    path_to_initial_transform)
                    cmd_args.append("--use-flirt 0")
                    cmd_args.append("--use-regaladin 1")
                    cmd_args.append("--test-ap-flip 0")

                cmd_args.append("--moving %s" % path_to_recon)
                cmd_args.append("--fixed %s" % path_to_reference)
                cmd_args.append("--dir-input %s" % dir_input_motion_correction)
                cmd_args.append("--dir-output %s" % dir_output)
                cmd_args.append("--write-transform 1")
                cmd_args.append("--verbose %d" % verbose)
                cmd = "niftymic_register_image %s" % (" ").join(cmd_args)

                ph.execute_command(
                    cmd,
                    flag_print_to_file=flag_batch_script,
                    path_to_file="%s%d.txt" %
                    (file_prefix_batch, ph.add_one(batch_ctr)))

        if flag_register_to_template_irtk:
            dir_input = utils.get_directory_case_recon_seg_mode(
                case_id=case_id,
                recon_space="subject_space",
                seg_mode="IRTK")
            dir_output = utils.get_directory_case_recon_seg_mode(
                case_id=case_id,
                recon_space="template_space",
                seg_mode="IRTK")
            path_to_recon = os.path.join(dir_input, "IRTK_SRR.nii.gz")
            path_to_reference = utils.get_path_to_recon(
                utils.get_directory_case_recon_seg_mode(
                    case_id=case_id,
                    recon_space="template_space",
                    seg_mode="seg_manual"),
                suffix="ResamplingToTemplateSpace",
            )
            path_to_initial_transform = os.path.join(
                utils.get_directory_case_recon_seg_mode(
                    case_id=case_id,
                    recon_space="template_space",
                    seg_mode="seg_manual"),
                "registration_transform_sitk.txt")

            cmd_args = []
            cmd_args.append("--fixed %s" % path_to_reference)
            cmd_args.append("--moving %s" % path_to_recon)
            cmd_args.append("--initial-transform %s" %
                            path_to_initial_transform)
            cmd_args.append("--use-fixed-mask 0")
            cmd_args.append("--use-moving-mask 0")
            cmd_args.append("--use-flirt 0")
            cmd_args.append("--use-regaladin 1")
            cmd_args.append("--test-ap-flip 0")
            cmd_args.append("--dir-output %s" % dir_output)
            cmd_args.append("--verbose %d" % verbose)
            cmd = "niftymic_register_image %s" % (" ").join(cmd_args)
            ph.execute_command(cmd)

        if flag_show_srr_template_space:
            recon_paths = []
            show_modes = list(SEG_MODES)
            # show_modes.append("IRTK")
            for seg_mode in show_modes:
                dir_input = utils.get_directory_case_recon_seg_mode(
                    case_id=case_id,
                    recon_space="template_space",
                    seg_mode=seg_mode)
                # # ------- DELETE -----
                # dir_input = re.sub("data", "foo+1", dir_input)
                # dir_input = re.sub(
                #     "volumetric_reconstruction/20180126/template_space/seg_auto",
                #     "", dir_input)
                # # -------
                path_to_recon_space = utils.get_path_to_recon(
                    dir_input, suffix="ResamplingToTemplateSpace",
                )
                recon_paths.append(path_to_recon_space)
            ph.show_niftis(recon_paths)
            ph.print_info("Sequence: %s" % (" -- ").join(show_modes))
            ph.pause()
            ph.killall_itksnap()

        # -----------------Reconstruct Volume in Template Space----------------
        if flag_reconstruct_volume_template_space:
            for seg_mode in SEG_MODES:
                path_to_recon_space = utils.get_path_to_recon(
                    utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="template_space",
                        seg_mode=seg_mode),
                    suffix="ResamplingToTemplateSpace",
                )
                dir_input = os.path.join(
                    utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="template_space",
                        seg_mode=seg_mode),
                    "motion_correction")
                dir_output = utils.get_directory_case_recon_seg_mode(
                    case_id=case_id,
                    recon_space="template_space",
                    seg_mode=seg_mode)
                # dir_output = os.path.join("/tmp/spina/template_space/%s-%s" % (
                #     case_id, seg_mode))

                cmd_args = []
                cmd_args.append("--dir-input %s" % dir_input)
                cmd_args.append("--dir-output %s" % dir_output)
                cmd_args.append("--reconstruction-space %s" %
                                path_to_recon_space)
                cmd_args.append("--alpha %s" % alpha)
                cmd_args.append("--verbose %s" % verbose)
                cmd_args.append("--use-masks-srr 0")

                # cmd_args.append("--minimizer L-BFGS-B")
                # cmd_args.append("--alpha 0.006")
                # cmd_args.append("--reconstruction-type HuberL2")
                # cmd_args.append("--data-loss arctan")
                # cmd_args.append("--iterations 5")
                # cmd_args.append("--data-loss-scale 0.7")

                cmd = "niftymic_reconstruct_volume_from_slices %s" % \
                    (" ").join(cmd_args)
                ph.execute_command(
                    cmd,
                    flag_print_to_file=flag_batch_script,
                    path_to_file="%s%d.txt" %
                    (file_prefix_batch, ph.add_one(batch_ctr)))

        # ----------------Collect SRR results in Template Space----------------
        if flag_collect_volumetric_reconstruction_results:
            directory = utils.get_directory_case_recon_summary(case_id)
            ph.create_directory(directory)

            # clear potentially existing files
            cmd = "rm -f %s/*.nii.gz" % (directory)
            ph.execute_command(cmd)

            # Collect SRRs
            for seg_mode in SEG_MODES:
                path_to_recon_src = utils.get_path_to_recon(
                    utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="template_space",
                        seg_mode=seg_mode),
                )
                path_to_recon = os.path.join(
                    directory, "%s%s.nii.gz" % (prefix_srr, seg_mode))

                cmd = "cp -p %s %s" % (path_to_recon_src, path_to_recon)
                ph.execute_command(cmd)

            # Collect IRTK recon
            path_to_recon_src = os.path.join(
                utils.get_directory_case_recon_seg_mode(
                    case_id=case_id,
                    recon_space="template_space",
                    seg_mode="IRTK"),
                "IRTK_SRR_LinearResamplingToTemplateSpace.nii.gz"
            )

            path_to_recon = os.path.join(
                directory, "%s%s.nii.gz" % (prefix_srr, "irtk"))

            cmd = "cp -p %s %s" % (path_to_recon_src, path_to_recon)
            ph.execute_command(cmd)

            # Collect evaluation mask
            path_to_recon = utils.get_path_to_recon(
                utils.get_directory_case_recon_seg_mode(
                    case_id=case_id,
                    recon_space="subject_space",
                    seg_mode="seg_auto"))

            template_stack_estimator = \
                tse.TemplateStackEstimator.from_mask(
                    ph.append_to_filename(path_to_recon, "_mask"))
            path_to_template = \
                template_stack_estimator.get_path_to_template()
            path_to_template_mask_src = ph.append_to_filename(
                path_to_template, "_mask_dil")
            path_to_template_mask = "%s/" % directory

            cmd = "cp -p %s %s" % (path_to_template_mask_src,
                                   path_to_template_mask)
            ph.execute_command(cmd)

        if flag_show_volumetric_reconstruction_results:
            dir_output = utils.get_directory_case_recon_summary(case_id)
            paths_to_recons = []
            for seg_mode in RECON_MODES:
                path_to_recon = os.path.join(
                    dir_output, "%s%s.nii.gz" % (prefix_srr, seg_mode))
                paths_to_recons.append(path_to_recon)
            path_to_mask = "%s/STA*.nii.gz" % dir_output
            cmd = ph.show_niftis(paths_to_recons,
                                 segmentation=path_to_mask)
            sitkh.write_executable_file([cmd], dir_output=dir_output)
            ph.pause()
            ph.killall_itksnap()

        # ---------------------Evaluate Image Similarities---------------------
        if flag_evaluate_image_similarities:
            dir_input = utils.get_directory_case_recon_summary(case_id)
            dir_output = utils.get_directory_case_recon_similarities(case_id)
            paths_to_recons = []
            for seg_mode in ["seg_auto", "detect", "irtk"]:
                path_to_recon = os.path.join(
                    dir_input, "%s%s.nii.gz" % (prefix_srr, seg_mode))
                paths_to_recons.append(path_to_recon)
            path_to_reference = os.path.join(
                dir_input, "%s%s.nii.gz" % (prefix_srr, "seg_manual"))
            path_to_reference_mask = utils.get_path_to_mask(dir_input)

            cmd_args = []
            cmd_args.append("--filenames %s" % " ".join(paths_to_recons))
            cmd_args.append("--reference %s" % path_to_reference)
            cmd_args.append("--reference-mask %s" % path_to_reference_mask)
            # cmd_args.append("--verbose 1")
            cmd_args.append("--dir-output %s" % dir_output)

            exe = re.sub("pyc", "py", os.path.abspath(
                evaluate_image_similarity.__file__))
            cmd_args.insert(0, exe)

            # clear potentially existing files
            cmd = "rm -f %s/*.txt" % (dir_output)
            ph.execute_command(cmd)

            cmd = "python %s" % " ".join(cmd_args)
            ph.execute_command(cmd)

        # -----------------Evaluate Slice Residual Similarities----------------
        if flag_evaluate_slice_residual_similarities:

            path_to_reference_mask = utils.get_path_to_mask(
                utils.get_directory_case_recon_summary(case_id))

            dir_output_root = \
                utils.get_directory_case_slice_residual_similarities(case_id)

            # clear potentially existing files
            # cmd = "rm -f %s/*.txt" % (dir_output_root)
            # ph.execute_command(cmd)

            for seg_mode in SEG_MODES:
                dir_input = os.path.join(
                    utils.get_directory_case_recon_seg_mode(
                        case_id=case_id,
                        recon_space="template_space",
                        seg_mode=seg_mode,
                    ),
                    "motion_correction")
                path_to_reference = os.path.join(
                    utils.get_directory_case_recon_summary(case_id),
                    "%s%s.nii.gz" % (prefix_srr, seg_mode))
                dir_output = os.path.join(dir_output_root, seg_mode)

                cmd_args = []
                cmd_args.append("--dir-input %s" % dir_input)
                cmd_args.append("--reference %s" % path_to_reference)
                cmd_args.append("--reference-mask %s" % path_to_reference_mask)
                cmd_args.append("--use-reference-mask 1")
                cmd_args.append("--use-slice-masks 0")
                # cmd_args.append("--verbose 1")
                cmd_args.append("--dir-output %s" % dir_output)

                exe = re.sub("pyc", "py", os.path.abspath(esrs.__file__))
                cmd_args.insert(0, exe)

                cmd = "python %s" % " ".join(cmd_args)
                ph.execute_command(cmd)

        # Collect data for blinded analysis
        if flag_collect_data_blinded_analysis:
            if flag_remove_failed_cases_for_analysis and case_id in RECON_FAILED_CASE_IDS:
                continue

            dir_input = utils.get_directory_case_recon_summary(case_id)
            # pattern = "STA([0-9]+)[_]mask.nii.gz"
            pattern = "STA([0-9]+)[_]mask_dil.nii.gz"
            p = re.compile(pattern)
            gw = [p.match(f).group(1)
                  for f in os.listdir(dir_input) if p.match(f)][0]

            dir_output = os.path.join(
                utils.get_directory_blinded_analysis(case_id, "open"),
                case_id)

            exe = re.sub("pyc", "py", os.path.abspath(
                mswm.__file__))

            recons = []

            for seg_mode in RECON_MODES:
                path_to_recon = os.path.join(
                    dir_input, "%s%s.nii.gz" % (prefix_srr, seg_mode))

                cmd_args = []
                cmd_args.append("--filename %s" % path_to_recon)
                cmd_args.append("--gestational-age %s" % gw)
                cmd_args.append("--dir-output %s" % dir_output)
                cmd_args.append("--prefix-output %s" % prefix_srr_qa)
                cmd_args.append("--verbose 0")
                cmd_args.insert(0, exe)

                cmd = "python %s" % " ".join(cmd_args)
                # ph.execute_command(cmd)

                recon = "%s%s" % (
                    prefix_srr_qa, os.path.basename(path_to_recon))
                recons.append(recon)
            ph.write_show_niftis_exe(recons, dir_output)

        if flag_anonymize_data_blinded_analysis:
            dir_input = os.path.join(
                utils.get_directory_blinded_analysis(case_id, "open"), case_id)
            dir_output_dictionaries = utils.get_directory_anonymized_dictionares(
                case_id)
            dir_output_anonymized_images = os.path.join(
                utils.get_directory_blinded_analysis(case_id, "anonymized"), case_id)

            if not ph.directory_exists(dir_input):
                continue
            ph.create_directory(dir_output_dictionaries)
            ph.create_directory(dir_output_anonymized_images)

            data_anonymizer = da.DataAnonymizer()
            # Create random dictionary (only required once)
            # data_anonymizer.set_prefix_identifiers("%s_" % case_id)
            # data_anonymizer.read_nifti_filenames_from_directory(dir_input)
            # data_anonymizer.generate_identifiers()
            # data_anonymizer.generate_randomized_dictionary()
            # data_anonymizer.write_dictionary(
            #     dir_output_dictionaries, "dictionary_%s" % case_id)

            # Read dictionary
            data_anonymizer.read_dictionary(
                dir_output_dictionaries, "dictionary_%s" % case_id)

            # Anonymize files
            if 0:
                ph.clear_directory(dir_output_anonymized_images)
                data_anonymizer.anonymize_files(
                    dir_input, dir_output_anonymized_images)

                # Write executable script
                filenames = ["%s.nii.gz" % f
                             for f in sorted(data_anonymizer.get_identifiers())]
                ph.write_show_niftis_exe(
                    filenames, dir_output_anonymized_images)

            # Reveal anonymized files
            if 1:
                filenames = data_anonymizer.reveal_anonymized_files(
                    dir_output_anonymized_images)
                filenames = sorted(["%s" % f for f in filenames])
                ph.write_show_niftis_exe(
                    filenames, dir_output_anonymized_images)

            # Reveal additional, original files
            # data_anonymizer.reveal_original_files(dir_output)

            # relative_directory = re.sub(
            #     utils.get_directory_blinded_analysis(case_id, "anonymized"),
            #     ".",
            #     dir_output_anonymized_images)
            # paths_to_filenames = [os.path.join(
            #     relative_directory, f) for f in filenames]

        # ---------------------Analyse Image Similarities---------------------
        if flag_analyse_image_similarities or \
                flag_analyse_slice_residual_similarities or \
                flag_analyse_stacks:
            if flag_remove_failed_cases_for_analysis:
                if case_id in RECON_FAILED_CASE_IDS:
                    continue
            if cases[case_id]['postrep'] == flag_postop or flag_postop == 2:
                cases_similarities.append(case_id)
                cases_stacks.append(utils.get_segmented_image_filenames(
                    case_id,
                    # subfolder=utils.SEGMENTATION_INIT,
                    subfolder=utils.SEGMENTATION_SELECTED,
                ))

        dir_output_analysis = os.path.join(
            # "/Users/mebner/UCL/UCL/Publications",
            "/home/mebner/Dropbox/UCL/Publications",
            "2018_MICCAI/brain_reconstruction_paper")

    if flag_analyse_image_similarities:
        dir_inputs = []
        filename = "image_similarities_postop%d.txt" % flag_postop
        for case_id in cases_similarities:
            dir_inputs.append(
                utils.get_directory_case_recon_similarities(case_id))
        cmd_args = []
        cmd_args.append("--dir-inputs %s" % " ".join(dir_inputs))
        cmd_args.append("--dir-output %s" % dir_output_analysis)
        cmd_args.append("--filename %s" % filename)

        exe = re.sub("pyc", "py", os.path.abspath(
            src.analyse_image_similarities.__file__))
        cmd_args.insert(0, exe)

        cmd = "python %s" % " ".join(cmd_args)
        ph.execute_command(cmd)

    if flag_analyse_slice_residual_similarities:
        dir_inputs = []
        filename = "slice_residuals_postop%d.txt" % flag_postop
        for case_id in cases_similarities:
            dir_inputs.append(
                utils.get_directory_case_slice_residual_similarities(case_id))
        cmd_args = []
        cmd_args.append("--dir-inputs %s" % " ".join(dir_inputs))
        cmd_args.append("--subfolder %s" % " ".join(SEG_MODES))
        cmd_args.append("--dir-output %s" % dir_output_analysis)
        cmd_args.append("--filename %s" % filename)

        exe = re.sub("pyc", "py", os.path.abspath(
            src.analyse_slice_residual_similarities.__file__))
        cmd_args.insert(0, exe)

        cmd = "python %s" % " ".join(cmd_args)
        # print len(cases_similarities)
        # print cases_similarities
        ph.execute_command(cmd)

    if flag_analyse_stacks:
        cases_stacks_N = [len(s) for s in cases_stacks]
        ph.print_subtitle("%d cases -- Number of stacks" % len(cases_stacks))
        ph.print_info("min: %g" % np.min(cases_stacks_N))
        ph.print_info("mean: %g" % np.mean(cases_stacks_N))
        ph.print_info("median: %g" % np.median(cases_stacks_N))
        ph.print_info("max: %g" % np.max(cases_stacks_N))

    elapsed_time = ph.stop_timing(time_start)
    ph.print_title("Summary")
    print("Computational Time for Pipeline: %s" %
          (elapsed_time))

    return 0


if __name__ == '__main__':
    main()
