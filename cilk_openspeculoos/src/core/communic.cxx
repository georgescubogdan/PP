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

// communic.cxx

#include "core/communic.hxx"

//_______________________________ Communicator ________________________________


Communicator :: Communicator (Mesh* m1, Mesh* m2, int typ)
{
    // Constructor.

# ifdef REQUIRE
    Require("valid argument 'typ'", typ == INTERNAL_COMMUNICATOR ||
	    typ == SEND_RECEIVE_COMMUNICATOR ||
	    typ == RECEIVE_COMMUNICATOR) ;
    InMethod("Communicator::Communicator(m1,m2,typ)") ;
# endif

    type  = typ ;
    mesh1 = m1 ;
    mesh2 = m2 ;
    pairs = new Vector<Element**>() ;

    CreatePairs() ;
}


Communicator :: ~Communicator ()
{
    // Destructor.

    int i ;

    for (i=1 ; i <= pairs->GetSize() ; i++)
	delete pairs->At(i) ;
    delete pairs ;
}


void Communicator :: CreatePairs ()
{
    // Identifies all elements of mesh1 and mesh2 which require inter-processor
    // communication.
    // Creates the pairs elem1/elem2 and stores them.

    Element *elem1, *elem2 ;
    Element **pair ;
    int     k, kk, p1, p2, dim1, dim2 ;
    boolean createPair ;

    dim1 = mesh1->GetDimension() ;
    dim2 = mesh2->GetDimension() ;

    if (dim1 >= dim2)

	for (k=1 ; k <= mesh1->GetSize() ; k++) {
	    elem1 = mesh1->At(k) ;

	    for (kk=1 ; kk <= mesh2->GetSize() ; kk++) {
		elem2 = mesh2->At(kk) ;

		if ((dim1==dim2 && elem1==elem2) || (dim1!=dim2 && elem1->Has(elem2))){
		    p1 = elem1->GetProcessorNumber() ; 
		    p2 = elem2->GetProcessorNumber() ; 
		    if (p1 == UNINITIALIZED_INTEGER || p2 == UNINITIALIZED_INTEGER)
			Error("Communicator::CreatePairs()",
			      "element with no assigned proc");

		    if (p1 == Mpi_rank || p2 == Mpi_rank) {

			if (type == INTERNAL_COMMUNICATOR)
			    createPair = (p1 == p2) ;

			else if (type == SEND_RECEIVE_COMMUNICATOR)
			    createPair = (p1 != p2) ;

			else if (type == RECEIVE_COMMUNICATOR)
			    createPair = (p1 != p2 && p2 == Mpi_rank) ;

			if (createPair) {
			    pair    = new Element*[2] ;
			    pair[0] = elem1 ;
			    pair[1] = elem2 ;
			    pairs->Put(pair) ;
			}
		    }
		}
	    }
	}

    else 

	for (k=1 ; k <= mesh2->GetSize() ; k++) {
	    elem2 = mesh2->At(k) ;

	    for (kk=1 ; kk <= mesh1->GetSize() ; kk++) {
		elem1 = mesh1->At(kk) ;

		if (elem2->Has(elem1)) {
		    p1 = elem1->GetProcessorNumber() ; 
		    p2 = elem2->GetProcessorNumber() ; 
		    if (p1 == UNINITIALIZED_INTEGER || p2 == UNINITIALIZED_INTEGER)
			Error("Communicator::CreatePairs()",
			      "element with no assigned proc");

		    if (p1 == Mpi_rank || p2 == Mpi_rank) {

			if (type == INTERNAL_COMMUNICATOR)
			    createPair = (p1 == p2) ;

			else if (type == SEND_RECEIVE_COMMUNICATOR)
			    createPair = (p1 != p2) ;

			else if (type == RECEIVE_COMMUNICATOR)
			    createPair = (p1 != p2 && p2 == Mpi_rank) ;

			if (createPair) {
			    pair    = new Element*[2] ;
			    pair[0] = elem1 ;
			    pair[1] = elem2 ;
			    pairs->Put(pair) ;
			}
		    }
		}
	    }
	}
}

    
void Communicator :: Print (int p)
{
    // Prints the receiver on standard output.
    // p = -1 means: for all processors.

    int i ;

    if (p == -1)
	p = Mpi_rank ;

    if (p == Mpi_rank) {

	if (type == INTERNAL_COMMUNICATOR)
	    printf ("\nInternal") ;
	else if (type == SEND_RECEIVE_COMMUNICATOR)
	    printf ("\nSend-receive") ;
	else
	    printf ("\nReceive") ;

	printf(" communicator (proc %d) from %dD mesh %X to %dD mesh %X:\n",
	       Mpi_rank, mesh1->GetDimension(),mesh1,mesh2->GetDimension(),mesh2);

	for (i=1 ; i <= pairs->GetSize() ; i++)
	    PrintPair(i) ;
    }
}


void Communicator :: PrintPair (int i)
{
    // Prints on standard output the 'i'-th pair of the receiver.

# ifdef REQUIRE
    Require("valid argument 'i'", i>=1 && i<=pairs->GetSize()) ;
    InMethod("Communicator::PrintPair(i)") ;
# endif

    Element **pair ;
    Point   *bary ;

    if (Mpi_rank >= 0) {
	pair = pairs->At(i) ;
	printf("  pair %d:\n",i) ;

	// elem1
	bary = pair[0]->GetBarycenter() ;
	printf ("    from elem on proc. %d with barycenter: ",
		pair[0]->GetProcessorNumber()) ;
	bary->Print() ;
	delete bary ;

	// elem2
	bary = pair[1]->GetBarycenter() ;
	printf ("\n    to   elem on proc. %d with barycenter: ",
		pair[1]->GetProcessorNumber()) ;
	bary->Print() ;
	delete bary ;
    }
}


