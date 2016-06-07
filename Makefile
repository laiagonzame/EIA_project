.PHONY : help
help: Makefile
	@sed -n 's/^##//p' $<

##compile : Compila i fa el "link"  del programa de dinamica.
.PHONY : compile
compile: main.o pbc.o integrator.o ini.o forces.o write_vmd.o link

##compile_time : Compila i fa el "link" per mirar els temps. 
.PHONY : compile_time
compile_time : main_time pbc_time integrator_time ini_time forces_time write_vmd_time link_time

##main.o : Compila el main (Sense links)
main.o : main_trajs.f90
	gfortran -c main_trajs.f90 -o main.o

##pbc.o : Compila el modul pbc.f90 
pbc.o : pbc.f90
	gfortran -c pbc.f90 -o pbc.o

##integrator.o : Compila el modul integrator.f90
integrator.o : integrator.f90
	gfortran -c integrator.f90 -o integrator.o

##ini.o : Compila el modul ini.f90
ini.o : ini.f90
	gfortran -c ini.f90 -o ini.o

##forces.o : Compila el modul forces.f90
forces.o : forces.f90
	gfortran -c forces.f90 -o forces.o

##write_vmd.o : Compila el modul write_vmd.f90
write_vmd.o : write_vmd.f90
	gfortran -c write_vmd.f90 -o write_vmd.o

##link : Link de  tots els moduls del programa de dinamica amb el main,genera l'executable amb el nom de van_der_walls
.PHONY : link
link:
	gfortran main.o ini.o pbc.o forces.o integrator.o write_vmd.o -o van_der_walls 
	@echo "L'executable s'anomena van_der_walls"


#Compilacio per veure els temps
##main_time : Compila el main per tal de poder estudiar els temps
.PHONY : main_temps
main_time : 
	gfortran -pg -c main_trajs.f90 -o main.o

##pbc_time : Compila el modul pbc.f90 per tal de poder estudiar el temps
.PHONY : pbc_time
pbc_time: 
	gfortran -pg -c pbc.f90 -o pbc.o

##integrator_time: Compila el modul integrator.f90 per tal de poder estudiar el temps
.PHONY : integrator_time
integrator_time : 
	gfortran -pg -c integrator.f90 -o integrator.o

##ini.o : Compila el modul ini.f90
.PHONY : ini_time
ini_time : 
	gfortran -pg -c ini.f90 -o ini.o

##forces.o : Compila el modul forces.f90
.PHONY: forces_time
forces_time : 
	gfortran -pg -c forces.f90 -o forces.o

##write_vmd.o : Compila el modul write_vmd.f90
.PHONY:write_vmd_time
write_vmd_time :
	gfortran -pg -c write_vmd.f90 -o write_vmd.o

##link_time : Link de  tots els moduls del programa de dinamica amb el main,genera l'executable amb el nom de van_der_walls
.PHONY : link_time
link_time:
	gfortran -pg main.o ini.o pbc.o forces.o integrator.o write_vmd.o -o van_der_walls 
	@echo "L'executable s'anomena van_der_walls"
	@echo "Cal executar el codi i cridar \"make time\" per tal de veure els temps"
.PHONY : time
time :
	gprof van_der_walls > time.out
	@echo "temps a time.out"
 
#Compilacio del programa que calcula les estadistiques

##clean : Esborra tots els .o generats
.PHONY : clean
clean:
	rm -f *.o
