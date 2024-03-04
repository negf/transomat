#$ -pe qmpi 96
#$ -cwd


export uhf=0

export ppn=32
export ppn_conquest=32
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMP_THREADS=1
export MKL_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=sockets

rm *.sh.* -f

#source ~/prog/mpich_gnu.sh
source ~/prog/ucx_intel.sh
#source ~/prog/libfabric_1.20.1_intel.sh
source ~/prog/petsc_3.20.3_intelmpi_intel.sh
#source ~/prog/petsc_3.19_omp_intelmpi_intel.sh
#source ~/prog/mpich_ucx_intel.sh

export UCX_TLS=self,knem,posix,dc_mlx5
ucx_shm_net_info.sh > ucx_shm_net_info.out 2>ucx_shm_net_info.out


cat $TMPDIR/machines | awk '{print $1 ":'$ppn'"}' > machines
cat $PE_HOSTFILE > pe_hostfile.out

unset DISPLAY
rm -f cleanup.out
for i in `cat pe_hostfile.out | awk '{print $1}'`
do
  echo $i >> cleanup.out
  ssh $i "/home/marius/bin/cleanup_processes.sh" -x >> cleanup.out 2>>cleanup.out
done

cat $TMPDIR/machines | awk '{print $1 ":'$ppn'"}' > machines
cat $PE_HOSTFILE > pe_hostfile.out

cat $TMPDIR/machines | awk '{print $1 ":'$ppn_conquest'"}' > machines_conquest
cat $PE_HOSTFILE > pe_hostfile.out

which mpirun > mpirun.which
which ucx_info > ucx_info.which

#export I_MPI_OFI_PROVIDER="psm3"
export I_MPI_OFI_PROVIDER="mlx"
export FI_PROVIDER=$I_MPI_OFI_PROVIDER


export I_MPI_EXTRA_FILESYSTEM=1
export I_MPI_EXTRA_FILESYSTEM_FORCE=nfs
export ROMIO_HINTS=/home/marius/tmp/hints.dat
export I_MPI_ADJUST_ALLREDUCE=6
#export ROMIO_PRINT_HINTS=1
#export ROMIO_FSTYPE_FORCE="nfs:"
rm negf.stop -f
./do_scf.sh >doscf.out 2>doscf.err
#transomat_bin=~/prog/transomat/transomat_add_channelswf_12012023/bin/transomat
#mpirun -ppn $ppn  $transomat_bin  > trans.out 2>trans.err
