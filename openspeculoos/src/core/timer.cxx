/*
* 
* This file is part of OpenSPECULOOS
* Copyright (c) 2008 Ecole Polytechnique Federale de Lausanne (EPFL)
*
* Contact: Nicolas Fietier, EPFL, Station 9, 1015 Lausanne, Switzerland
* E-mail:  nicolas.fietier@epfl.ch
*
* OpenSPECULOOS is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* OpenSPECULOOS is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with OpenSPECULOOS.  If not, see <http://www.gnu.org/licenses/>.
* 
*/

// timer.cxx

#include "core/timer.hxx"
#include "core/util.hxx"
#include "core/options.hxx"
#include <unistd.h>


Timer :: Timer ()
{

}

Timer :: Timer (string aname)
{
    // Constructor. Creates a timer with title 'aname'.

    nameTimer  = aname ;
    cpu_clktck = (float) sysconf(_SC_CLK_TCK) ;
    isActive   = false ;

    tot_etime             = 0 ;
    tot_cstime.tms_utime  = 0 ;
    tot_cstime.tms_stime  = 0 ;
    tot_cstime.tms_cutime = 0 ;
    tot_cstime.tms_cstime = 0 ;
}


Timer :: ~Timer ()
{
    // Destructor.

}


void Timer :: Print ()
{
    // Prints the measured times of the receiver on standard output.

# ifdef REQUIRE
    Require("inactive timer", ! isActive) ;
    InMethod("Timer::Print()") ;
# endif

    float tot_cpu, tot_sys, tot_ela ;

    tot_cpu = (float) (tot_cstime.tms_utime + tot_cstime.tms_cutime) / 
	cpu_clktck ;
    tot_sys = (float) (tot_cstime.tms_stime + tot_cstime.tms_cstime) / 
	cpu_clktck ;
    tot_ela = (float) tot_etime / cpu_clktck ;

    printf("\n ooo Timer '%s'\n",nameTimer.c_str()) ;
    printf("   CPU: %6.1f | SYSTEM: %6.1f | ELAPSED: %6.1f\n\n",
	   tot_cpu,tot_sys,tot_ela) ;
}


void Timer :: PrintIn (FILE* file)
{
    // Prints the measured times of the receiver on standard output.

# ifdef REQUIRE
    Require("inactive timer", ! isActive) ;
    InMethod("Timer::PrintIn(file)") ;
# endif

    float tot_cpu, tot_sys, tot_ela ;

    tot_cpu = (float) (tot_cstime.tms_utime + tot_cstime.tms_cutime) / 
	cpu_clktck ;
    tot_sys = (float) (tot_cstime.tms_stime + tot_cstime.tms_cstime) / 
	cpu_clktck ;
    tot_ela = (float) tot_etime / cpu_clktck ;

    fprintf(file,"CPU: %6.1f | SYSTEM: %6.1f | ELAPSED: %6.1f\n\n",
	    tot_cpu,tot_sys,tot_ela) ;
}

void Timer :: GetTimeFloat (float &tot_cpu, float &tot_sys,float &tot_ela)
{
    // Return the measured times of the receiver

# ifdef REQUIRE
    Require("inactive timer", ! isActive) ;
    InMethod("Timer::GetTIMEFLOAT (file)") ;
# endif

    //  float tot_cpu, tot_sys, tot_ela ;

    tot_cpu = (float) (tot_cstime.tms_utime + tot_cstime.tms_cutime) ;
    tot_sys = (float) (tot_cstime.tms_stime + tot_cstime.tms_cstime) ; 
    tot_ela = (float) tot_etime ;
}

float Timer :: GetCPU ()
{
    // Return the measured times of the receiver
    float tot_cpu;

# ifdef REQUIRE
    Require("inactive timer", ! isActive) ;
    InMethod("Timer::GetTIMEFLOAT (file)") ;
# endif

    //  float tot_cpu, tot_sys, tot_ela ;

    tot_cpu = (float) (tot_cstime.tms_utime + tot_cstime.tms_cutime)/ 
	cpu_clktck ; ;
    return tot_cpu;
}


void Timer :: Start ()
{
    // Starts the clock.

# ifdef REQUIRE
    Require("inactive timer", ! isActive) ;
    InMethod("Timer::Start()") ;
# endif

    start_etime = times(&start_cstime) ;
    isActive    = true ;
}


void Timer :: Stop ()
{
    // Stops the clock.

# ifdef REQUIRE
    Require("active timer", isActive) ;
    InMethod("Timer::Stop()") ;
# endif

    stop_etime = times(&stop_cstime) ;

    tot_etime             += stop_etime - start_etime ;
    tot_cstime.tms_utime  += stop_cstime.tms_utime - start_cstime.tms_utime ;
    tot_cstime.tms_stime  += stop_cstime.tms_stime - start_cstime.tms_stime ;
    tot_cstime.tms_cutime += stop_cstime.tms_cutime - start_cstime.tms_cutime ;
    tot_cstime.tms_cstime += stop_cstime.tms_cstime - start_cstime.tms_cstime ;

    isActive = false ;
}


void Timer :: PrintSeconds (){
    // Prints the measured times of the receiver on standard output
    // in seconds
    // Vincent Keller add-on

# ifdef REQUIRE
    Require("inactive timer", ! isActive) ;
    InMethod("Timer::Print()") ;
# endif
    // clock_t	
    float tot_cpu, tot_sys, tot_ela ;
	
    tot_cpu = (float) (tot_cstime.tms_utime + tot_cstime.tms_cutime) / cpu_clktck ;
    tot_sys = (float) (tot_cstime.tms_stime + tot_cstime.tms_cstime) / cpu_clktck ;
    tot_ela = (float) tot_etime / cpu_clktck ;
	
    printf("\n ooo Timer '%s'\n",nameTimer.c_str()) ;
    printf("   CPU: %6.1f | SYSTEM: %6.1f | ELAPSED: %6.1f\n\n",tot_cpu,tot_sys,tot_ela) ;

}


#ifdef PARALLEL

ParallelTimer :: ParallelTimer ()
{
    start_time = 0.0;
    stop_time = 0.0;
    inter_time = 0.0;
    Start();
}


ParallelTimer :: ~ParallelTimer ()
{
    Stop();
}


void ParallelTimer :: Start() 
{
    start_time = MPI_Wtime();
}


void ParallelTimer :: Stop() 
{
    stop_time = MPI_Wtime();
}


void ParallelTimer :: Init()
{
    start_time = 0.0;
}


double ParallelTimer :: GetIntermediateTime()
{
    inter_time = MPI_Wtime();
    return (inter_time - start_time);
}


double ParallelTimer :: GetIntermediateTimeInit()
{
    inter_time = MPI_Wtime();
    double tmp = start_time;
    start_time = 0.0;
    return (inter_time - tmp);
}


void ParallelTimer :: WriteEndLine() 
{
    if (MPI_Comm_rank == 0){
	std::ofstream file( "results.txt", std::ios_base::app );
	file << "\n";
    }
}

void ParallelTimer :: WriteStartLine()
{
    if (MPI_Comm_rank == 0){
	std::ofstream file( "results.txt", std::ios_base::app );
	file << MPI_Comm_size(MPI_COMM_WORLD, &Mpi_nbProcessors) << "\t";
    }
}


void ParallelTimer :: WriteTimeToFile(double time_to_add)
{
    if (MPI_Comm_rank == 0){
	std::ofstream file( "results.txt", std::ios_base::app );
	file << time_to_add << "\t";
    }
}


void ParallelTimer :: WriteTimeToFile()
{
    if (MPI_Comm_rank == 0){
	std::ofstream file( "results.txt", std::ios_base::app );
	file << MPI_Wtime() << "\t";
    }
}

#endif


