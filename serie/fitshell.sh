for K in 100 200 300 400 500 600 700 800
do
  #rm -r $K
  mkdir $K
  cp *.f90 $K
  cp compile $K
  cd $K
  mkdir data
  sed -i -e 's/M=100/M='$K'/g'  main_trajs.f90
  ./compile
  ./van_der_waals
  cd ..
done
