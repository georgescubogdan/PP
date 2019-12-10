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

// parallel.cxx

#include "core/parallel.hxx"
#include <new>


#define SIZE_2D 20000    // hopefully the max nb of values in any communication
                         // buffer. Can be increased
#define MAX_PROC 4096       // can be increased

#ifdef PARALLEL

// temporary initialization. True initialization will take place in 
// Mpi_Initialize()

int          Mpi_rank               = -1 ;
int          Mpi_nbProcessors       = -1 ;
Mpi_Buffer** Mpi_sendBuffers        = NULL ;
Mpi_Buffer** Mpi_receiveBuffers     = NULL ;
int*         Mpi_receiveBufferSizes = NULL ;

#else   // scalar implementation

// true initialization. Allows the program to execute correctly even if the
// users forgets to invoke Mpi_Initialize

int          Mpi_rank               = 0 ;
int          Mpi_nbProcessors       = 1 ;
Mpi_Buffer** Mpi_sendBuffers        = NULL ;     // not used
Mpi_Buffer** Mpi_receiveBuffers     = NULL ;     //   id.
int*         Mpi_receiveBufferSizes = NULL ;     //   id.

#endif


#ifdef PARALLEL

void Mpi_Barrier ()
{
    // Synchronizes all processors.

    int ierr, ierr2 ;

    ierr  = MPI_Barrier(MPI_COMM_WORLD) ;
    ierr2 = MPI_SUCCESS ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Barrier()","MPI failure 1") ;
    if (ierr2 != MPI_SUCCESS)
	Error("Mpi_Barrier()","MPI failure 2") ;
# endif
}


void Mpi_Communicate ()
{
    // 1) Sends to every processor p the values in the buffer Mpi_sendBuffers[p].
    // 2) Initializes every buffer Mpi_receiveBuffers[p] from values sent by
    //    processor p.

    int i, n, p ;
    static int  rrr=1, sss=1 ;
    int qqq=1, ppp=1;

    int taille ;
    int myTag;

    MPI_Request request[2 * MAX_PROC] ;
    MPI_Status  status[2 * MAX_PROC] ;
    Mpi_Buffer *sbuffer ;
    Mpi_Buffer *rbuffer ;

    int count;
    int ierr;

# ifdef CHECK
    if (Mpi_nbProcessors > MAX_PROC)
	Error("Mpi_Communicate()","too small value of MAX_PROC") ;
# endif

//Mpi_PrintBuffers() ;

// 1. Send the buffer to the connected processors

    count = 0;

    myTag = Mpi_rank ;

    for (p=0 ; p< Mpi_nbProcessors ; p++) 
    {
	sbuffer = Mpi_sendBuffers[p] ;
	taille = sbuffer->GetSize();
	if (taille != 0)
	{

	    ierr = MPI_Isend(sbuffer->GetValues(), taille, MPI_DOUBLE, 
			     p, myTag, MPI_COMM_WORLD, &(request[count])) ;
	    count++;

# ifdef CHECK
	    if (ierr != MPI_SUCCESS)
		Error("Mpi_ISend","error status returned") ;
# endif

	}
    }
  
// 2. Receive the buffers from the connected processors

    for (p=0 ; p<Mpi_nbProcessors ; p++) {
	taille = Mpi_receiveBufferSizes[p] ;
	if (taille > 0) {
	    rbuffer = Mpi_receiveBuffers[p] ;
	    rbuffer->Resize(taille) ;
      
	    ierr = MPI_Irecv(rbuffer->GetValues(), taille, MPI_DOUBLE, p, p, MPI_COMM_WORLD, &(request[count]) ) ;
	    count++;

//      Mpi_IReceive(taille,rbuffer->GetValues(),1,p,&rrequest[rrequestcount++]) ;
# ifdef CHECK
	    if (ierr != MPI_SUCCESS)
		Error("Mpi_IReceive","error status returned") ;
# endif

	}
    }
  
    MPI_Waitall(count,request, status) ;
}


void Mpi_Finalize ()
{
    // Finalizes MPI. 
    // This function must be the last MPI function to be called in a program, so
    // it is typically called at the end of the main() function.

    int ierr ;

    ierr = MPI_Finalize() ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Finalize()","MPI failure") ;
# endif
}

void Mpi_Initialize (int* argc, char*** argv)
{
    // Initializes MPI.
    // This function must be the first MPI function to be called in a program, so
    // it is typically called at the beginning of the main() function.
    //
    // It also initialize function set_new_handler.

    int p, ierr, np ;

    // 2. Initialize the communication

    ierr = MPI_Init(argc, argv) ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Initialize(argc,argv)","MPI failure") ;
# endif

    // 3. Initialize Mpi_nbProcessors

    ierr = MPI_Comm_size(MPI_COMM_WORLD, &Mpi_nbProcessors) ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Initialize()","Failure when getting nb processors") ;
# endif

    // 4. Initialize Mpi_rank

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &Mpi_rank) ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Initialize()","Failure when getting the processor's rank") ;
# endif

    np = Mpi_nbProcessors ;
    if (Mpi_rank==0) printf("\nOK Mpi_Initialize()    nb Processors=%d\n",np) ;

    // 5. Initialize Mpi_sendBuffers, Mpi_receiveBuffers, etc

    Mpi_sendBuffers        = new Mpi_Buffer* [Mpi_nbProcessors] ;
    Mpi_receiveBuffers     = new Mpi_Buffer* [Mpi_nbProcessors] ;
    Mpi_receiveBufferSizes = new int         [Mpi_nbProcessors] ;
    for (p=0 ; p<Mpi_nbProcessors ; p++) {
	Mpi_sendBuffers[p]        = new Mpi_Buffer() ;
	Mpi_receiveBuffers[p]     = new Mpi_Buffer() ;
	Mpi_receiveBufferSizes[p] = 0 ;
    }

    // 7. Initialize function set_new_handler

    set_new_handler(&ExhaustedMemory) ;
}


real Mpi_Max (real x)
{
    // Returns the largest among the contributions 'x' of every processor.
    // The same value is returned to all processors.

    real answer ;
    int  ierr ;

    ierr = MPI_Allreduce(&x, &answer, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Max(x)","MPI failure") ;
# endif

    return answer ;
}

void Mpi_PrintBuffers (int p)
{
    // Prints the buffer information for processor p (all processors if p = -1).

    int i ;

    if (p == Mpi_rank || p == -1) {
	printf("\nBuffers of processor %d:\n",Mpi_rank) ;

	printf("  Send buffers:\n") ;
	for (i=0 ; i<Mpi_nbProcessors ; i++)
	    if ((Mpi_sendBuffers[i])->GetSize() != 0) {
		printf("    to proc %d: ",i); (Mpi_sendBuffers[i])->Print() ;
	    }

	printf("  Size of receive buffers:\n") ;
	for (i=0 ; i<Mpi_nbProcessors ; i++)
	    if (Mpi_receiveBufferSizes[i] != 0)
		printf("    from proc %d: %d\n",i,Mpi_receiveBufferSizes[i]) ;

	printf("  Receive buffers:\n") ;
	for (i=0 ; i<Mpi_nbProcessors ; i++)
	    if ((Mpi_receiveBuffers[i])->GetSize() != 0) {
		printf("    proc %d: ",i); (Mpi_receiveBuffers[i])->Print() ;
	    }
    }
}


void Mpi_Receive (int size, real* values, int incr, int proc)
{
    // Receives 'size' real values from processor 'proc' and initializes 'values'
    // to them. The initialized coeffs of 'values' are:
    //   values[0], values[incr], values[2*incr], etc. 
    // The stride 'incr' can be negative. 

# ifdef REQUIRE
    Require("involves 2 processors", Mpi_rank != proc) ;
    Require("valid processor number", proc >= 0 && proc < Mpi_nbProcessors) ;
    Require("valid argument 'incr'", incr != 0) ;
    Require("'size' not too large", size <= SIZE_2D) ;
    InMethod("Mpi_Receive(size,values,incr,proc)") ;
# endif

    MPI_Status   status ;
    real         buffer[SIZE_2D] ;
    int          ione, ierr ;

    if (incr==1)
	ierr = MPI_Recv(values, size, MPI_DOUBLE, proc, 99, MPI_COMM_WORLD, 
			&status) ;
    else {
	ierr = MPI_Recv(buffer, size, MPI_DOUBLE, proc, 99, MPI_COMM_WORLD, 
			&status) ;
	ione = IONE ;
	DCOPY(&size, &buffer[0], &ione, &values[0], &incr) ;
    }

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Receive","error status returned") ;
# endif
}

void Mpi_ReceiveInt (int* value, int proc)
{
    // Receives 1 integer value from processor 'proc' and initializes 'value' to 
    // it.

# ifdef REQUIRE
    Require("involves 2 processors", Mpi_rank != proc) ;
    Require("valid processor number", proc >= 0 && proc < Mpi_nbProcessors) ;
    InMethod("Mpi_ReceiveInt(value,proc)") ;
# endif

    MPI_Status   status ;
    int          ione, ierr ;

    ierr = MPI_Recv(value, 1, MPI_INT, proc, 99, MPI_COMM_WORLD, &status) ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Receive","MPI error 1") ;
# endif
}

void Mpi_ResetCommunication ()
{
    // Clears out the inter-processor buffers.

    int p ;

# ifdef CHECK
    // check that last time communication occured, the receive buffer was
    // completely read
    int size, start ;
    for (p=0 ; p<Mpi_nbProcessors ; p++) {
	size  = (Mpi_receiveBuffers[p])->GetSize() ;
	start = (Mpi_receiveBuffers[p])->GetStart() ;
	if (start != size)
	    Error("Mpi_ResetCommunication()",
		  "error last time Mpi_receiveBuffers[p] was read") ;
    }
# endif

    for (p=0 ; p<Mpi_nbProcessors ; p++) {
	(Mpi_sendBuffers[p])   ->Reset() ;
	(Mpi_receiveBuffers[p])->Reset() ;
	Mpi_receiveBufferSizes[p] = 0 ;
    }
}


void Mpi_Send (int size, real* values, int incr, int proc)
{
    // Sends to processor 'proc' 'size' coeffs from 'values'.
    // The involved coeffs are:
    //   values[0], values[incr], values[2*incr], etc.
    // The stride 'incr' can be negative. 

# ifdef REQUIRE
    Require("involves 2 processors", Mpi_rank != proc) ;
    Require("valid processor number", proc >= 0 && proc < Mpi_nbProcessors) ;
    Require("valid argument 'incr'", incr != 0) ;
    Require("'size' not too large", size <= SIZE_2D) ;
    InMethod("Mpi_Send(size,values,incr,proc)") ;
# endif

    real buffer[SIZE_2D] ;
    int  ione, ierr ;

    if (incr == 1) 
	ierr = MPI_Send(values, size, MPI_DOUBLE, proc, 99, MPI_COMM_WORLD) ;
    else {
	ione = IONE ;
	DCOPY(&size, &values[0], &incr, &buffer[0], &ione) ;
	ierr = MPI_Send(buffer, size, MPI_DOUBLE, proc, 99, MPI_COMM_WORLD) ;
    }

# ifdef CHECK
    if (ierr)
	Error("Mpi_Send","error status reported") ;
# endif
}


void Mpi_SendReceive(int sendsize, real* sendvalues, int incr, int destproc, 
                     int recvsize, real* recvvalues, int sourceproc)
{
    // Sends to processor 'proc' 'size' coeffs from 'values'.
    // The involved coeffs are:
    //   values[0], values[incr], values[2*incr], etc.
    // The stride 'incr' can be negative.

# ifdef REQUIRE
//  Require("involves 2 diferent processors", destproc != sourceproc) ;
    Require("valid send    processor number", destproc >= 0 && destproc < Mpi_nbProcessors) ;
    Require("valid receive processor number", sourceproc >= 0 && sourceproc < Mpi_nbProcessors) ;
    Require("valid argument 'incr'", incr != 0) ;
    Require("'sendsize' not too large", sendsize <= SIZE_2D) ;
    Require("'recvsize' not too large", recvsize <= SIZE_2D) ;
    InMethod("Mpi_SendReceive(size,sendvalues,incr,destproc,recvvalues,sourceproc)") ;
# endif

    real sendbuffer[SIZE_2D] ;
    real recvbuffer[SIZE_2D] ;
    int  ione, ierr ;
    MPI_Status   status ;

    if (incr == 1)
	ierr = MPI_Sendrecv(sendvalues, sendsize, MPI_DOUBLE, destproc, 99, 
			    recvvalues, recvsize, MPI_DOUBLE, sourceproc, 99, MPI_COMM_WORLD, &status); 
    else {
	ione = IONE ;
	DCOPY(&sendsize, &sendvalues[0], &incr, &sendbuffer[0], &ione) ;
	ierr = MPI_Sendrecv(sendbuffer, sendsize, MPI_DOUBLE, destproc, 99, 
			    recvbuffer, recvsize, MPI_DOUBLE, sourceproc, 99, MPI_COMM_WORLD, &status);
	DCOPY(&recvsize, &recvbuffer[0], &ione, &recvvalues[0], &incr) ;
    }

# ifdef CHECK
    if (ierr)
	Error("Mpi_SendReceive","error status reported") ;
# endif
}


void Mpi_SendInt (int* value, int proc)
{
    // Sends to processor 'proc' 1 coeff from 'value'.

# ifdef REQUIRE
    Require("involves 2 processors", Mpi_rank != proc) ;
    Require("valid processor number", proc >= 0 && proc < Mpi_nbProcessors) ;
    InMethod("Mpi_SendInt(value,proc)") ;
# endif

    real buffer[SIZE_2D] ;
    int  ione, ierr ;
 
    ierr = MPI_Send(value, 1, MPI_INT, proc, 99, MPI_COMM_WORLD) ;

# ifdef CHECK
    if (ierr)
	Error("Mpi_SendInt","error status reported") ;
# endif
}

void Mpi_SendVectInt (int size, int* value,int incr, int proc)
{
    // Sends to processor 'proc' 'size' coeffs from 'values'.
    // The involved coeffs are:
    //   values[0], values[incr], values[2*incr], etc.
    // The stride 'incr' can be negative. 

# ifdef REQUIRE
    Require("involves 2 processors", Mpi_rank != proc) ;
    Require("valid processor number", proc >= 0 && proc < Mpi_nbProcessors) ;
    InMethod("Mpi_SendInt(value,proc)") ;
# endif

    //  real buffer[SIZE_2D] ;
    int  ione, ierr ;
 
    ierr = MPI_Send(value, size, MPI_INT, proc, 99, MPI_COMM_WORLD) ;

# ifdef CHECK
    if (ierr)
	Error("Mpi_SendInt","error status reported") ;
# endif
}


int Mpi_Sum (int n)
{
    // Returns the sum of the contributions 'n' of every processor.
    // The same value is returned to all processors.

    int answer, ierr ;

    ierr = MPI_Allreduce(&n, &answer, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Sum(n)","MPI failure") ;
# endif

    return answer ;
}

real Mpi_Sum (real x)
{
    // Returns the sum of the contributions 'x' of every processor.
    // The same value is returned to all processors.

    real answer ;
    int  ierr ;

    ierr = MPI_Allreduce(&x, &answer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;

# ifdef CHECK
    if (ierr != MPI_SUCCESS)
	Error("Mpi_Sum(x)","MPI failure") ;
# endif

    return answer ;
}

real Mpi_Wtime ()
{
    // Returns the elapsed wall-clock time (in seconds) since some time in the
    // past.
    // In order to measure an elapsed time, call this function at the beginning
    // and at the end of the code to measure, then make the difference of the two
    // obtained values.

    return MPI_Wtime() ;
}

//-----------------------------------------------------------------------------

#else                 //******* scalar implementation *********//

//-----------------------------------------------------------------------------

void Mpi_Barrier ()
{
    // Default implementation for a scalar computer.
    // Does nothing.
}

void Mpi_Bcast (fpos_t*, int)
{

    // Default implementation for a scalar computer.
    // Does nothing.
}

void Mpi_BcastAll (real*, int)
{

    // Default implementation for a scalar computer.
    // Does nothing.
}

real* Mpi_GatherTo0(real* distmini)
{
    return distmini;
}

real* Mpi_Rassemble0(real* data, int tailledata)
{
    return data;
}

void Mpi_Communicate ()
{
    // Default implementation for a scalar computer.
    // Does nothing.
}



void Mpi_Finalize ()
{
    // Default implementation for a scalar computer.
    // Does nothing.
}


void Mpi_Initialize (int*, char***)
{
    // Default implementation for a scalar computer.
    // This function is preferably called at the beginning of function "main".
    // Failing to call this function is harmless, except that the new_handler is
    // not defined.

    set_new_handler(&ExhaustedMemory) ;
}


real Mpi_Max (real x)
{
    // Default implementation for a scalar computer.
    // Returns 'x'.

    return x ;
}


void Mpi_PrintBuffers (int)
{ 
}

void Mpi_Receive (int, real*, int, int)
{
    // Issues an error message, because a scalar implementation should never
    // invoke this function.

    Error("Mpi_Receive(int,real*,int,int)","Nonsense") ;
}

void Mpi_ReceiveInt (int* value, int proc)
{
// Issues an error message, because a scalar implementation should never
    // invoke this function.

    Error("Mpi_ReceiveInt (int* value, int proc)","Nonsense") ;
}
void Mpi_ReceiveVectInt (int size, int* value, int incr, int proc)
{
    // Issues an error message, because a scalar implementation should never
    // invoke this function.

    Error("Mpi_ReceiveVectInt (int size, int* value, int incr, int proc)","Nonsense") ;
}

void Mpi_SendInt (int* , int )
{
    // Issues an error message, because a scalar implementation should never
    // invoke this function.

    Error("Mpi_ReceiveInt(real*,int)","Nonsense") ;
}


void Mpi_ResetCommunication ()
{
    // Default implementation for a scalar computer.
    // Does nothing.
}


void Mpi_Send (int, real*, int, int)
{
    // Issues an error message, because a scalar implementation should never
    // invoke this function.

    Error("Mpi_Send(int,real*,int,int)","Nonsense") ;
}


int Mpi_Sum (int n)
{
    // Default implementation for a scalar computer.
    // Returns 'n'.

    return n ;
}


real Mpi_Sum (real x)
{
    // Default implementation for a scalar computer.
    // Returns 'x'.

    return x ;
}



real Mpi_Wtime ()
{
    // Default implementation for a scalar computer.
    // Returns 0.

    return 0. ;
}


#endif

