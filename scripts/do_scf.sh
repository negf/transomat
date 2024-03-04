date > scf_start_end.txt
date +%s >> scf_start_end.txt
for i in `seq 1 2000` ; do
echo "----------" $i "-----------"
  if [ -f "negf.converged" ] || [ -f negf.stop ] ; then
    date >> scf_start_end.txt
    date +%s >> scf_start_end.txt
    rm -f do_conquest_step
    exit
  else
    { time ./scf.sh > scf.out 2>scf.err ;} 2>> scf_timing.txt
    err=`echo $?`
    if [ $err -ne 0 ] ; then
      echo "SCF failed " $err
      date >> scf_start_end.txt
      date +%s >> scf_start_end.txt
      exit $err
    fi
  fi
done

