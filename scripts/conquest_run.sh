#$ -pe qmpi 64
#$ -cwd



export ppn=32
export ppn_conquest=32
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMP_THREADS=1
export MKL_THREADS=1
#export OMP_PROC_BIND=close
#export OMP_PLACES=sockets

rm *.sh.* -f

#source ~/prog/mpich_gnu.sh
source ~/prog/ucx_intel.sh
#source ~/prog/petsc_3.17.4_intel.sh
#source /home/marius/prog/petsc_3.17.4_mpich_intel.sh
#source ~/prog/mpich_ucx_intel.sh

export UCX_TLS=sysv,knem,dc_mlx5
#export UCX_TLS=posix,sysv,ud_mlx5



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

export I_MPI_OFI_PROVIDER="psm3"
#export I_MPI_OFI_PROVIDER="mlx"
export FI_PROVIDER=$I_MPI_OFI_PROVIDER


export I_MPI_EXTRA_FILESYSTEM=1
export I_MPI_EXTRA_FILESYSTEM_FORCE=nfs
export ROMIO_HINTS=/home/marius/tmp/hints.dat
#export ROMIO_PRINT_HINTS=1
#export ROMIO_FSTYPE_FORCE="nfs:"

mpirun -machinefile machines ~/prog/conquest/CONQUEST-release_merge/src_negf_uhf/Conquest > conquest.out 2>conquest.err
