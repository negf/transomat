ecc_dir=`pwd`
transomat_bin=~/prog/transomat/transomat_17012024_add_re_interror/bin/transomat
conquest_bin=~/prog/conquest/CONQUEST-release_merge/src_negf_uhf_fix0.5/Conquest
#conquest_bin=~/prog/conquest/CONQUEST-release_merge/src_negf_uhf/Conquest

rm  subdiv* -f
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
if [ ! -f pulay.negf ] ; then
echo "0" > pulay.negf
echo "-1" >> pulay.negf
echo "-1" >> pulay.negf
fi

i=`head -n 1 pulay.negf | awk '{print $1}'`

if [[ ( ! -f ./negf_"$i"/transp.petsc && ! -f do_conquest_step) ]]
then
  mkdir negf_"$i"
  scp * negf_"$i"
fi

if [ ! -f do_conquest_step ] ; then
  { time mpirun -machinefile ./machines  $transomat_bin transp_a.ini > trans_a.out 2>trans_a.err ;} 2>> transomat_timing.txt
fi
err=`echo $?`
if [ $err -ne 0 ] ; then
  echo "transomat failed " $err
  exit 1
fi

if [ $uhf -eq 1 ] ; then
  if [ ! -f do_conquest_step ] ; then
    { time mpirun -machinefile ./machines  $transomat_bin transp_b.ini > trans_b.out 2>trans_b.err ;} 2>> transomat_timing.txt
  fi
  err=`echo $?`
  if [ $err -ne 0 ] ; then
    echo "transomat failed " $err
    exit 1
  fi
fi


touch do_conquest_step

if [ ! -f negf.converged ]
then

  export OMP_NUM_THREADS=$OMP_THREADS
  export MKL_NUM_THREADS=$MKL_THREADS


  if [ ! -f ./conquest_"$i"/Conquest_input ]
  then
    mkdir conquest_"$i"
    scp * conquest_"$i"
  fi

  { time mpirun -machinefile ./machines_conquest  $conquest_bin > conquest.out 2>conquest.err ;} 2>> conquest_timing.txt

  err=`echo $?`
  if [ $err -ne 0 ] ; then
    echo "Conquest failed " $err
    exit 2
  fi
  rm -f do_conquest_step

  i=`tail -n 1 convergence.negf | awk '{print $1}'`

  plt21d_8 hartree_pot_diff.plt z -1 -1  > h_"$i".dat
  plt21d_8 hartree_pot_diff.plt z 0 0  > h_"$i"_ave.dat
  plt21d_8 td.plt z -1 -1  > d_"$i".dat
  plt21d_8 td_diff.plt z -1 -1  > dd_"$i".dat
  plt21d_8 td.plt z  0 0  > d_"$i"_ave.dat
  plt21d_8 td_diff.plt z 0 0  > dd_"$i"_ave.dat
  plt21d_8 density_u_0.dat z 0 0 > d_u_"$i"_ave.dat
  plt21d_8 rho_negf_u.dat z 0 0 > r_n_u_"$i"_ave.dat

  scp trans_a.out trans_a.out."$i"
  scp trans_a.err trans_a.err."$i"
  if [ $uhf -eq 1 ] ; then
    scp trans_b.out trans_b.out."$i"
    scp trans_b.err trans_b.err."$i"
    plt21d_8 density_d_0.dat z 0 0 > d_d_"$i"_ave.dat
  fi

  scp conquest.out conquest.out."$i"
  scp conquest.err conquest.err."$i"
  scp Conquest_out Conquest_out."$i"

else
 echo "NEGF converged"
fi
