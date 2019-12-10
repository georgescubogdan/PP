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

// postpro.cxx

#include "core/postpro.hxx"
#include "core/point.hxx"
#include "core/mesh.hxx"
#include "core/element.hxx"
#include "core/field.hxx"
#include "core/parent.hxx"
#include "core/options.hxx"
#include "core/parallel.hxx"
#include <cstring>


//______________________________ PostProcessor ________________________________


PostProcessor :: PostProcessor (FlatField* coord, boolean writeCoord)
{
    // Constructor. Initializes the receiver to a postprocessor with 'coord' as
    // coordinate field.

    string name ;
    int    i ;

    coordinates    = coord ;
    fields         = new Vector<FlatField*>() ;
    components     = new Vector<int>() ;
    componentNames = new Vector<string>() ;

    // store the coordinate components in 'fields'
    if (writeCoord){
	for (i=1 ; i<=coord->GetNbComponents() ; i++) {
	    if (i == 1)
		name="X" ;
	    else if (i == 2)
		name="Y" ;
	    else if (i == 3)
		name="Z" ;
	    else
		Error("PostProcessor::PostProcessor(coord)","nbComponents(coord) > 3") ;
	    Put(coord,i,name) ;
	}
    }

}


PostProcessor :: ~PostProcessor ()
{
    // Destructor.

    int i ;

    for (i=1 ; i<=fields->GetSize() ; i++) {
	unref(fields->At(i)) ;
    }

    delete fields ;
    delete components ;
}


void PostProcessor :: Put (FlatField* x, int icx, string cName)
{
    // Adds the 'icx'-th component of 'x' to the list of field components to be
    // written. Assigns 'cName' as name to that component.

# ifdef REQUIRE
    Require("'x' has valid mesh", x->GetMesh() == coordinates->GetMesh()) ;
    Require("'x' has valid interpolation", 
	    x->HasSameInterpolationAs(coordinates)) ;
    Require("valid argument 'icx'", icx >= 1 && icx <= x->GetNbComponents()) ;
    Require("Tecplot admits no blanks in 'cName'", strchr(cName,' ') == NULL) ;
    InMethod("PostProcessor::Put(x,icx,cName)") ;
# endif

    ref(x) ;
    fields->Put(x) ;
    components->Put(icx) ;
    componentNames->Put(cName) ;
}


//__________________________________ Tecplot __________________________________


Tecplot :: Tecplot (string fname, FlatField* coord)
    : PostProcessor (coord)
{
    // Constructor. Initializes the receiver to a postprocessor with file name
    // 'fname' and coordinate field 'coor'. The file format default is Ascii.

    fileName = fname ;
    title    = "" ;
}


Tecplot :: ~Tecplot ()
{
    // Destructor.

}


void Tecplot :: GenerateOutput ()
{
    // Output the fields to the file.

    Initialize() ;
    WriteZones() ; 
}


void Tecplot :: Initialize ()
{
    // Prints header of preplot file.
   
    FILE   *file ;
    string name ;
    int  i, n ;

    if (Mpi_rank == 0) {

      file = Fopen(fileName.c_str(),"w") ;

	if (!title.empty())
	  fprintf(file,"TITLE = \"%s\"\n",title.c_str()) ;

	fprintf(file,"VARIABLES = ") ;

	n = componentNames->GetSize() ;
	for (i=1 ; i<=n ; i++) { 
	    name = componentNames->At(i) ;
	    if (strcmp(name.c_str(),""))
	      fprintf(file," %s",name.c_str()) ; 
	    else
		fprintf(file," (no name)") ;
	    if (i != n)
		fprintf(file,",") ;
	    else
		fprintf(file,"\n") ;
	} 

	fclose(file) ;
    }

    Mpi_Barrier() ;
}



void Tecplot :: SetTitle (string aTitle)
{
    // Sets the title of the receiver to 'aTitle'.

  if (!aTitle.empty())
	title = aTitle ;
}


void Tecplot :: WriteZones () 
{
    // Output all zones in the file.
    // The file is opened and closed at every element, in order to allow the
    // processors to have access to it sequentially, otherwise a zone written by
    // a processor might be overwritten by another one.

    Mesh                *mesh ;
    Element             *elem ;
    ElementaryField     *e ;
    Vector<RealVector*> *arrays ;
    RealVector          *val ;
    FILE                *file ;
    int                 dim, i, j, k, n, ic, index, nbFields, nbElem ;

    mesh   = coordinates->GetMesh() ;
    dim    = mesh->GetDimension() ;
    nbElem = mesh->GetSize() ;

    for (k=1 ; k<=nbElem ; k++) {
	elem = mesh->At(k) ;

	// print zone header

	if (elem->GetProcessorNumber() == Mpi_rank) {

	    e = elem->GetElementaryField(coordinates->GetIndex()) ;
      
	    file = Fopen(fileName.c_str(),"a") ;
	    fprintf(file,"ZONE T = \"field\"\n") ;

	    if (dim >= 1)
		fprintf(file,"I = %d",e->GetParentElement(1)
			->GetNbCollocationPoints()) ;
	    if (dim >= 2)
		fprintf(file,", J = %i",e->GetParentElement(2)
			->GetNbCollocationPoints()) ;
	    if (dim >= 3)
		fprintf(file,", K = %i",e->GetParentElement(3)
			->GetNbCollocationPoints()) ;
     
	    fprintf(file,", F = POINT\n") ;
  
	    // get the components of the elementary field (including the coordinates)
  
	    nbFields = fields->GetSize() ;
	    arrays   = new Vector<RealVector*>(nbFields) ;
	    for (i=1 ; i<=nbFields ; i++) {
		index = fields->At(i)->GetIndex() ;
		ic    = components->At(i) ;
		val   = elem->GetElementaryField(index)->GetComponent(ic) ;
		arrays->At(i) = val ;
	    }
  
	    // print the components, dof by dof
  
	    n = arrays->At(1)->GetSize() ;
	    for (j=1 ; j<=n ; j++) {
		for (i=1 ; i<=nbFields ; i++)
		    fprintf(file,"%12.4E",arrays->At(i)->At(j)) ;
		fprintf(file,"\n") ;
	    }
  
	    delete arrays ;

	    fclose(file) ;
	}

	Mpi_Barrier() ;
    }
}

void Tecplot :: ReadOutput ()
{
    // Output the fields to the file.

    FILE *file ;

    file = Fopen(fileName.c_str(),"r") ;

    InitializeRead(file) ;
    ReadZones(file) ; 
  
    fclose(file) ;
}

void Tecplot :: InitializeRead (FILE *file)
{
    // Prints header of preplot file.
   
    string name ;
    int    i, n ;

    // saute les 2 premieres lignes
    char ligne[MAXWORD];
    fgets(ligne,MAXWORD-1,file);  // avoid the first line
    fgets(ligne,MAXWORD-1,file);  // avoid the 2nd line

}

void Tecplot :: ReadZones (FILE *file) 
{
    // Output all zones in the file.
    // The file is opened and closed at every element, in order to allow the
    // processors to have access to it sequentially, otherwise a zone written by
    // a processor might be overwritten by another one.

    Mesh                *mesh ;
    Element             *elem ;
    ElementaryField     *e ;
    Vector<RealVector*> *arrays ;
    RealVector          *val ;
    int                 dim, i, j, k, n, ic, index, nbFields, nbElem ;
    char ligne[MAXWORD];

    mesh   = coordinates->GetMesh() ;
    dim    = mesh->GetDimension() ;
    nbElem = mesh->GetSize() ;

    for (k=1 ; k<=nbElem ; k++) {
	elem = mesh->At(k) ;

	// print zone header

	if (elem->GetProcessorNumber() == Mpi_rank) {

	    e = elem->GetElementaryField(coordinates->GetIndex()) ;
      
	    fgets(ligne,MAXWORD-1,file);  // avoid the line ZONE T= etc

	    int tmp1, tmp2, tmp3;
	    if (dim == 1){
		fscanf(file,"I = %d, F = POINT\n", &tmp1);  printf("tmp1 = %d \n", tmp1);
		if (tmp1 != e->GetParentElement(1)->GetNbCollocationPoints()){
		    Error("Tecplot:: readZone","tmp1 != nbcoolloc pts ");
		}
	    }

	    if (dim == 2){
		fscanf(file,"I = %d, J = %d, F = POINT\n", &tmp1, &tmp2);

		if (tmp1 != e->GetParentElement(1)->GetNbCollocationPoints())
		    Error("Tecplot:: readZone","tmp1 != nbcoolloc pts ");
		if (tmp2 != e->GetParentElement(2)->GetNbCollocationPoints())
		    Error("Tecplot:: readZone","tmp2 != nbcoolloc pts ");
	    }

	    if (dim == 3){
		fscanf(file,"I = %d, J = %d, K = %d, F = POINT\n", &tmp1, &tmp2, &tmp3);

		if (tmp1 != e->GetParentElement(1)->GetNbCollocationPoints())
		    Error("Tecplot:: readZone","tmp1 != nbcoolloc pts ");
		if (tmp2 != e->GetParentElement(2)->GetNbCollocationPoints())
		    Error("Tecplot:: readZone","tmp2 != nbcoolloc pts ");
		if (tmp3 != e->GetParentElement(3)->GetNbCollocationPoints())
		    Error("Tecplot:: readZone","tmp3 != nbcoolloc pts ");
	    }
     
 
	    // get the components of the elementary field (including the coordinates)
  
	    nbFields = fields->GetSize() ;
	    arrays   = new Vector<RealVector*>(nbFields) ;
	    for (i=1 ; i<=nbFields ; i++) {
		index = fields->At(i)->GetIndex() ;
		ic    = components->At(i) ;
		val   = elem->GetElementaryField(index)->GetComponent(ic) ;
		arrays->At(i) = val ;
	    }
  
	    // print the components, dof by dof
	    real aa;

	    n = arrays->At(1)->GetSize() ;
	    for (j=1 ; j<=n ; j++) {
		for (i=1 ; i<=nbFields ; i++){
		    //          fscanf(file,"%12.4E",aa);
		    fscanf(file,"%lf",&aa);
		    arrays->At(i)->At(j) = aa ;
		}
		fscanf(file,"\n") ;
	    }
  
	    delete arrays ;

	}

	Mpi_Barrier() ;
    }
}
