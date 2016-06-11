.PHONY : help
help: Makefile
	@sed -n 's/^##//p' $<

##compile : Compila i fa el "link"  del programa de dinamica.
.PHONY : compile
compile: main.o pbc.o integrator.o ini.o forces.o write_vmd.o link

##compile_time : Compila i fa el "link" per mirar els temps. 
.PHONY : compile_time
compile_time : main_time pbc_time integrator_time ini_time forces_time write_vmd_time link_time

#main.o : Compila el main (Sense links)
main.o : ./serie/main_trajs.f90
	gfortran -c ./serie/main_trajs.f90 -o main.o

#pbc.o : Compila el modul pbc.f90 
pbc.o : ./serie/pbc.f90
	gfortran -c ./serie/pbc.f90 -o pbc.o

#integrator.o : Compila el modul integrator.f90
integrator.o : ./serie/integrator.f90
	gfortran -c ./serie/integrator.f90 -o integrator.o

#ini.o : Compila el modul ini.f90
ini.o : ./serie/ini.f90
	gfortran -c ./serie/ini.f90 -o ini.o

#forces.o : Compila el modul forces.f90
forces.o : ./serie/forces.f90
	gfortran -c ./serie/forces.f90 -o forces.o

#write_vmd.o : Compila el modul write_vmd.f90
write_vmd.o : ./serie/write_vmd.f90
	gfortran -c ./serie/write_vmd.f90 -o write_vmd.o

#link : Link de tots els moduls i el main
.PHONY : link
link:
	gfortran main.o ini.o pbc.o forces.o integrator.o write_vmd.o -o van_der_waals 
	@echo "L'executable s'anomena van_der_waals"


#Compilacio per veure els temps
#main_time : Compilacio del main amb -pg
.PHONY : main_temps
main_time : 
	gfortran -pg -c ./serie/main_trajs.f90 -o main.o

#pbc_time : Compila el modul pbc.f90 amb -pg
.PHONY : pbc_time
pbc_time: 
	gfortran -pg -c ./serie/pbc.f90 -o pbc.o

#integrator_time: Compila el modul integrator.f90 amb -pg
.PHONY : integrator_time
integrator_time : 
	gfortran -pg -c ./serie/integrator.f90 -o integrator.o

#ini.o : Compila el modul ini.f90 amb -pg
.PHONY : ini_time
ini_time : 
	gfortran -pg -c ./serie/ini.f90 -o ini.o

#forces.o : Compila el modul forces.f90 amb -pg
.PHONY: forces_time
forces_time : 
	gfortran -pg -c ./serie/forces.f90 -o forces.o

#write_vmd.o : Compila el modul write_vmd.f90 amb -pg
.PHONY:write_vmd_time
write_vmd_time :
	gfortran -pg -c ./serie/write_vmd.f90 -o write_vmd.o

#link_time : Link de tots els moduls que permet emprear gprof
.PHONY : link_time
link_time:
	gfortran -pg main.o ini.o pbc.o forces.o integrator.o write_vmd.o -o van_der_waals 
	@echo "L'executable s'anomena van_der_waals"
	@echo "Cal executar el codi i cridar \"make time\" per tal de veure els temps"

##time : execució del codi i gprof
.PHONY : time
time :    
	./van_der_waals
	gprof van_der_waals > time.out
	@echo "temps a time.out"
 

#Compilacio del programa que calcula les estadistiques
#compilacio_estadistiques: compilacio dels moduls d'estadistiques
#.PHONY : compilacio_estadistiques
#compilacio_estadistiques: stadistics.o temperature.o

#statistics.o :./serie/statistics.f90
#	gfortran -c ./serie/statistics.f90 -o statistics.o

#temperature.o :./serie/temperature.f90
#	gfortran -c ./serie/temperature.f90 -o temperature.0


##clean : Esborra tots els .o generats
.PHONY : clean
clean:
	rm -f *.o
