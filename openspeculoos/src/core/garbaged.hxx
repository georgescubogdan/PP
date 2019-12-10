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

#ifndef GARBAGED
#define GARBAGED

/*!
 * \file grabaged.hxx
 * \brief Reference counting objects
 * \author {Lots of authors}
 * \version 1.0
 */


/*! \class GarbagedObject
   * \brief Manages reference counting of objects.
   * 
   * Functions ref and unref are not declared "inline" only for the sake of
   * debugging with the DDE debugger. They are good candidates for in-line
   * expansion.
   */ 
class GarbagedObject
{ 

public :

    int counter ;

public:

    // = CONSTRUCTORS

    inline GarbagedObject () ;
    inline virtual ~GarbagedObject () ;

    // = REFERENCE COUNTING
 
    friend  void  ref (GarbagedObject* obj) ;       // a function, not a method
    friend  void  unref (GarbagedObject* obj) ;     // a function, not a method
    virtual void BreakCyclicReferences () ;

    // = DEBUGGING
 
    void         PrintCounter () ;
} ; 


//___________________________ GarbagedObject (inline) _________________________


GarbagedObject :: GarbagedObject ()
{
    // Constructor.

    counter = 1 ;
}


GarbagedObject :: ~GarbagedObject ()
{
    // Destructor.
}


#endif

