for K in 100 200 300 400 500 600 700 800
do
  cd $K
  gprof van_der_waals > ../data/prof$K.txt
  cd ..
done
