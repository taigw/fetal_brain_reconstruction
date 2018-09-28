export PYTHONPATH=/Users/guotaiwang/Documents/workspace/fetalSR/ITK_NiftyMIC-build/Wrapping/Generators/Python:/Users/guotaiwang/Documents/workspace/fetalSR/ITK_NiftyMIC-build/lib
PATH=${PATH}:/Users/guotaiwang/CIVTK/install/niftyreg_install/bin:/Users/guotaiwang/Documents/workspace/fetalSR/ITK_NiftyMIC-build/bin:/Users/guotaiwang/CIVTK/install/c3d_install/bin
export THEANO_FLAGS=mode=FAST_RUN,device=cuda,floatX=float32
python reconstruction_multi_patients.py
