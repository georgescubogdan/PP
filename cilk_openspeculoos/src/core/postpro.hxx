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

#ifndef POSTPRO
#define POSTPRO

/*!
 * \file postpro.hxx
 * \brief Manages post-processing
 * \author {Lots of authors}
 * \version 1.0
 */

#include "core/vector.hxx"
#include "core/vector_templates.hxx"
#include <stdio.h>
#include <string>

class Element ;
class FlatField ;
enum DomainType {MONO, MULTI} ;

/*! \class PostProcessor
   * \brief Manages post-processors output.
   * 
   * A post processor is essentially a list of fields. These field must be
   * defined on the same mesh and must have identical interpolation.
   * 
   * The main purpose of a post-processor object is to write down a file
   * containing the fields in a format readable by a post-processing 
   * program.
   * 
   * Attribute <{coordinates}> stores the coordinate field of the mesh.
   * 
   * To increase flexibility, fields are printed by component. Information
   * on what to print is contained in 3 arrays of equal size:
   * - <{fields}> contains the fields to be printed. A field may appear
   * more than once;
   * - <{components}>: components(i) is the component number associated
   * with the i-th field in <{fields}>;
   * - <{componentNames}>: componentNames(i) is the title associated with
   * the i-th component field.
   */ 


class PostProcessor
{ 

protected :

    FlatField*          coordinates ;
    Vector<FlatField*>* fields ;
    Vector<int>*        components ;
    Vector<string>*      componentNames ;

public :

    PostProcessor (FlatField* coord, boolean writeX = true) ;
    virtual ~PostProcessor () ;

    void          Put (FlatField* f, int ic, string compName) ;
    virtual void  GenerateOutput () = 0 ;
} ;

/*! \class Tecplot
   * \brief Manages Tecplot preprocessor (Preplot) output.
   * 
   * Each spectral element defines a zone. Coordinates and fields are 
   * written in a single Ascii file <{fileName}>. A <{title}> can be
   * supplied.
   * 
   * To create a Tecplot data file:
   * 1) create and objet of type Tecplot, 
   * 2) assign it one or more components of one or more fields, using
   * method 'Put',
   * 3) dump the values, using method 'GenerateOutput',
   * 4) delete the Tecplot object.
   * 
   * If the values of only one field are to be written in the file, you can
   * equivalently (and more simply) use method FlatField::TecplotDump.
   * 
   * This class is suited for both serial and parallel implementations. In
   * the latter case, all processors obviously write their data into the
   * same file.
   */ 
class Tecplot : public PostProcessor
{

protected :

    string fileName ;
    string title ;

    void Initialize () ;
    void WriteZones () ;
    void InitializeRead (FILE *file) ;
    void ReadZones (FILE *file) ;
    
public :

    Tecplot (string fname, FlatField* coord) ;
    ~Tecplot () ;
   
    void GenerateOutput () ;
    void ReadOutput () ;
    void SetTitle (string aTitle) ;
} ;


#endif
