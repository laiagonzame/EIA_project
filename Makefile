#Makefile creada per Cristina Roncero Barrero 
#per al projecte de l'assignatura d'Eines informatiques avançades

.PHONY : help
help: Makefile
	@sed -n 's/^##//p' $<

##compile : Compila i fa el "link"  del programa de dinamica.
.PHONY : compile
compile:  serie/pbc.o serie/integrator.o serie/ini.o serie/forces.o serie/write_vmd.o serie/main.o link

##compile_time : Compila i fa el "link" per mirar els temps. 
.PHONY : compile_time
compile_time : pbc_time integrator_time ini_time forces_time write_vmd_time main_time link_time 

##compile_mpi : Compila el programa paral·lelitzat 
.PHONY: compile_mpi
compile_mpi: paralel/forces-parallel.o paralel/initial.o paralel/integrator_para.o paralel/pbc_paralel.o paralel/write_vmd.o paralel/main_trajs.o link_mpi 

##estadistiques : Executa totes les estadistiques 
.PHONY: estadistiques
estadistiques: gr msd scala ener_temp

##time : execució del codi en serie (només despres de compile_time) i crida a gprof
.PHONY : time
time :    
	./van_der_waals
	gprof van_der_waals > time.out
	@echo "temps a time.out"


##login : login al cluster per compilar
.PHONY: login
login:
	@echo "caldra fer \"module load openmpi/1.4.2_intel-11.1.072 \""
	qrsh -q cerqt2.q -pe smp 1

#-----------------------------------------------------------------------------------------------#

#Compilacio dels diferents moduls, individual
#main.o : Compila el main (Sense links)
serie/main.o : ./serie/main_trajs.f90
	gfortran -c ./serie/main_trajs.f90 -o ./serie/main.o

#pbc.o : Compila el modul pbc.f90 
serie/pbc.o : ./serie/pbc.f90
	gfortran -c ./serie/pbc.f90 -o ./serie/pbc.o

#integrator.o : Compila el modul integrator.f90
serie/integrator.o : ./serie/integrator.f90
	gfortran -c ./serie/integrator.f90 -o ./serie/integrator.o

#ini.o : Compila el modul ini.f90
serie/ini.o : ./serie/ini.f90
	gfortran -c ./serie/ini.f90 -o ./serie/ini.o

#forces.o : Compila el modul forces.f90
serie/forces.o : ./serie/forces.f90
	gfortran -c ./serie/forces.f90 -o ./serie/forces.o

#write_vmd.o : Compila el modul write_vmd.f90
serie/write_vmd.o : ./serie/write_vmd.f90
	gfortran -c ./serie/write_vmd.f90 -o ./serie/write_vmd.o

#link : Link de tots els moduls i el main
.PHONY : link
link:
	gfortran ./serie/main.o ./serie/ini.o ./serie/pbc.o ./serie/forces.o ./serie/integrator.o ./serie/write_vmd.o -o van_der_waals 
	@echo "L'executable s'anomena van_der_waals"


#Compilacio individual amb la  per veure els temps
#main_time : Compilacio del main amb -pg
.PHONY : main_temps
main_time : 
	gfortran -pg -c ./serie/main_trajs.f90 -o ./serie/main.o

#pbc_time : Compila el modul pbc.f90 amb -pg
.PHONY : pbc_time
pbc_time: 
	gfortran -pg -c ./serie/pbc.f90 -o ./serie/pbc.o

#integrator_time: Compila el modul integrator.f90 amb -pg
.PHONY : integrator_time
integrator_time : 
	gfortran -pg -c ./serie/integrator.f90 -o ./serie/integrator.o

#ini.o : Compila el modul ini.f90 amb -pg
.PHONY : ini_time
ini_time : 
	gfortran -pg -c ./serie/ini.f90 -o ./serie/ini.o

#forces.o : Compila el modul forces.f90 amb -pg
.PHONY: forces_time
forces_time : 
	gfortran -pg -c ./serie/forces.f90 -o ./serie/forces.o

#write_vmd.o : Compila el modul write_vmd.f90 amb -pg
.PHONY:write_vmd_time
write_vmd_time :
	gfortran -pg -c ./serie/write_vmd.f90 -o ./serie/write_vmd.o

#link_time : Link de tots els moduls que permet emprear gprof
.PHONY : link_time
link_time:
	gfortran -pg ./serie/main.o ./serie/ini.o ./serie/pbc.o ./serie/forces.o ./serie/integrator.o ./serie/write_vmd.o -o van_der_waals 
	@echo "L'executable s'anomena van_der_waals"
	@echo "Cal ejecutar el codi amb \"make time\" per tal de veure els temps"



#Compilacio del codi paral·lel, dins del cluster
#forces-parallel.o : Compila el modul de forces paral·lel
paralel/forces-parallel.o : ./paralel/forces-parallel.f90
	mpif90 -g -c ./paralel/forces-parallel.f90 -o ./paralel/forces-parallel.o

#integrator_para.o : Compila el modul del d'integració en paral·lel
paralel/integrator_para.o : ./paralel/integrator_para.f90
	mpif90 -g -c ./paralel/integrator_para.f90 -o ./paralel/integrator_para.o

#pbc_paralel.o : Compila el modul de PBC paral·lelitzat
paralel/pbc_paralel.o : ./paralel/pbc_paralel.f90
	mpif90 -g -c ./paralel/pbc_paralel.f90 -o ./paralel/pbc_paralel.o

#initial.o : Compila el modul d'inicialització paral·lel
paralel/initial.o: ./paralel/initial.f90
	mpif90 -g -c ./paralel/initial.f90 -o ./paralel/initial.o

#write_vmd.o : Compila el modul write_vmd paral·lel
paralel/write_vmd.o : ./paralel/write_vmd.f90
	mpif90 -g -c ./paralel/write_vmd.f90 -o ./paralel/write_vmd.o

#main_trajs.o : Compila el main paral·lel
paralel/main_trajs.o : ./paralel/main_trajs.f90
	mpif90 -g -c  ./paralel/main_trajs.f90 -o ./paralel/main_trajs.o
#link de tots els moduls del codi paral·lelitzat. 
.PHONY: link_mpi
link_mpi:
	mpif90 -g  ./paralel/initial.o ./paralel/forces-parallel.o ./paralel/pbc_paralel.o ./paralel/integrator_para.o ./paralel/main_trajs.o  ./paralel/write_vmd.o -o van_der_waals_paralel
	@echo "executable file: van_der_waals_paralel"


#---------------------------------------------------------------------------------------#

#Compilació i execució dels programes que calculan les estadistiques

.PHONY:  Estadistiques
Estadistiques: gr msd 

.PHONY: gr
gr:
	ipython ./analysis/gr.py 

.PHONY: msd
msd:
	ipython ./analysis/msd.py 

.PHONY: ener_temp
ener_temp:
	ipython ./analysis/ener_temp.py 

.PHONY: scala
scala::
	ipython ./analysis/scala.py 



#--------------------------------------------------------------------------------------#

##clean : Esborra tots els *.o i *.mod generats
.PHONY : clean
clean:
	rm -f ./paralel/*.o
	rm -f ./serie/*.o
	rm -f *.mod
