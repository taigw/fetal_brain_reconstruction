import numpy as np
import os
import pysitk.python_helper as ph

mask_names = ['manual', 'auto', 'detect']
srr_names  = ['baseline', # no outlier rejection and gaussian process regularization
                     'outlier',  # outlier rejection but not gaussian process regularization
                     'outlier_gpr'] # outlier rejection and gaussian process regularization
space_names = ['subject', 'template']


patients = ['57_2']
ages     = [25]
mask_indices = [2]
srr_indices  = [2]
space_indices = [1]

def main():
    for patient_i in range(len(patients)):
        for mask_i in mask_indices:
            for srr_i in srr_indices:
                patient = patients[patient_i]
                method_mask = mask_names[mask_i]
                method_srr  = srr_names[srr_i]
                dir_input = "{0:}/input_preprocess".format(patient)
                if(not os.path.isdir(dir_input)):
                    dir_input = "{0:}/input".format(patient)
                cmd_args = []
                cmd_args.append(" --dir-input {0:}".format(dir_input))
                cmd_args.append(" --dir-mask {0:}/mask_{1:}/mask".format(patient, method_mask))

                dir_output = "{0:}/mask_{1:}/reconstruct_{2:}".format(patient, method_mask, method_srr)
                cmd_args.append(" --dir-output {0:}".format(dir_output))
                
                patient_prefix = patient
                if("_" in patient):
                    patient_prefix = patient.split("_")[0]
                subject_space_arg = 1 if 0 in space_indices else 0
                template_space_arg = 1 if 1 in space_indices else 0
                cmd_args.append(" --gestational-age {0:}".format(ages[patient_i]))
                cmd_args.append(" --suffix-mask _mask")
                cmd_args.append(" --alpha 0.02")
                cmd_args.append(" --bias-field-correction 0")
                cmd_args.append(" --run-data-vs-simulated-data 0")
                cmd_args.append(" --run-recon-subject-space {0:}".format(subject_space_arg))
                cmd_args.append(" --run-recon-template-space {0:}".format(template_space_arg))
                
                outlier_rejct   = 0
                gaussian_smooth = 0
                if(srr_i == 1):
                    outlier_rejct = 1
                elif(srr_i == 2):
                    outlier_rejct = 1
                    gaussian_smooth = 1
                cmd_args.append(" --outlier-rejection {0:}".format(outlier_rejct))
                cmd_args.append(" --use-robust-registration {0:}".format(gaussian_smooth))
                cmd = "python reconstruction.py {0:}".format(' '.join(cmd_args))
                exit_code = ph.execute_command(cmd)

if __name__ == '__main__':
    main()
