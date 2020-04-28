git clone git@xgitlab.cels.anl.gov:darshan/darshan.git
cd darshan
git checkout dev-detailed-hdf5-mod


module swap craype-mic-knl craype-haswell

module swap PrgEnv-intel PrgEnv-gnu

module load cray-hdf5-parallel/1.10.5.2

./configure CC=cc --with-mem-align=8 --prefix=/usr/projects/eap/users/gshipman/darshan-3.1.8/install --with-log-path-by-env=DARSHAN_OUTPUT_DIR,SLURM_SUBMIT_DIR,PWD --with-jobid-env=SLURM_JOBID --enable-hdf5-mod CFLAGS='-I /opt/cray/pe/hdf5-parallel/1.10.5.2/GNU/8.2/include' LDFLAGS='-L /opt/cray/pe/hdf5-parallel/1.10.5.2/GNU/8.2/lib -lhdf5'


make
make install
module use /usr/projects/eap/users/gshipman/darshan/install/share/craype-2.x/modulefiles
module load darshan/3.2.0-pre1