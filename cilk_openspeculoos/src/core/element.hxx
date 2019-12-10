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

#ifndef ELEMENT
#define ELEMENT
/*!
 * \file element.hxx
 * \brief Management of spectral elements
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/garbaged.hxx"
#include "core/field.hxx"
#include "core/param.hxx"

class Vertex ;
class Edge ;
class Face ;
class Volume ;
class Point ;
class Distribution ;
class ElementaryField ;

enum ReceiveMode { SIZE_RM, VALUES_RM, UNDEFINED_RM } ;


/*! \class Element
   * \brief Base class managing spectral elements of any dimension.
   * 
   * Objects of derived types (edges, faces, hexas) can be used either as
   * geometrical supports of spectral elements or as "macroelements" for
   * generation. In the former case, attribute <{distributions>} is not
   * used. In the latter case, it contains a set 1, 2 or 3 distributions.
   * 
   * The element stores a number of elementary fields, which are its 
   * contributions to the fields defined on the domains the element belongs 
   * to.
   * 
   * Q operators deriving from the mortar element technique are also
   * managed internally within spectral elements (implemented in 
   * subclasses).
   * 
   * Attribute <{procNumber}> stores the number of the processor to which
   * the element is attached. All processors store every element. This 
   * means that every processor has complete knowledge of the meshes, 
   * fields and elementary fields. The only restriction is the VALUES of
   * the elementary fields, which, for a given element, are stored only by
   * one processor, attribute <{procNumber}> of that element. 
   * On a scalar computer, <{procNumber}> is equal to 0 for all elements.
   * 
   * Remark. The element has no parent elements; since these may differ 
   * from field to field, they are defined in every elementary field.
   * 
   * \see {Vertex, Edge, Face, Volume}
   */

class Element : public GarbagedObject
{ 

protected :

    Vector<Distribution*>*    distributions ;
    Vector<ElementaryField*>* elementaryFields ; 
    int                       procNumber ;

public :

    // = CONSTRUCTORS

    Element () ; 
    virtual ~Element () ;

    // = NATURE

    virtual boolean  IsVertex () ;
    virtual boolean  IsEdge () ;
    virtual boolean  IsFace () ;
    virtual boolean  IsVolume () ;

    // = POINT DISTRIBUTIONS

    Distribution*    GetDistribution (int) ;
    void             SetDistribution (int, Distribution*) ;
    void             SetDistribution (int, int) ;
    void             InitializeDistributions () ;

    // = MANAGEMENT OF THE SUBELEMENTS

    virtual int      GetNbSubelements () ;
    virtual int      GetNbVertices () = 0 ;
    virtual Element* GetSubelement (int i) ;
    virtual Vertex*  GetVertex (int i) = 0 ;

    // = OTHER TOPOLOGICAL OPERATIONS

    virtual int      GetDimension () ;
    virtual boolean  Has (Element* elem) ;
    virtual boolean  IsIdenticalTo (Element* elem) ;
    inline boolean   IsIncludedIn (Element* elem) ;
    virtual int      SearchLocationOf (Element* elem) = 0 ;
   
    // = MANAGEMENT OF THE ELEMENTARY FIELDS

    inline ElementaryField* GetElementaryField (int index) ;
    void             CreateElementaryField (int index, int ncomp, 
                                            Vector<ParentElement*>* vp);
    void             DeleteElementaryField (int index) ;
    void             DuplicateElementaryField (int index1, int index2) ;
    void             DuplicateEmptyElementaryField (int index1, int index2,
                                                    int ncomp) ;
    inline boolean   HasElementaryField (int index) ;

    // = COMPUTATION OF THE ELEMENTARY BASIC FIELDS

    virtual void     SetToCoordinates (int index) = 0 ;

    // = COPY, ADD ELEMENTARY FIELDS FROM ANOTHER ELEMENT
    //     These functions are for use only from methods of class Field.
    //     See the corresponding section in class Field.

    void             CopyAddValuesFrom   (Element* elemX, 
                                          int indX, int indY, int icx, int icy,
                                          CopyAddMode cam, InterpMode im,
                                          NonConform nc, 
                                          ReceiveMode rm=UNDEFINED_RM) ;
    virtual void     CopyAddValuesFromVertex (Vertex* elemX, 
					      int indX, int indY, int icx, int icy,
					      CopyAddMode cam, InterpMode im, 
					      NonConform nc,
					      ReceiveMode rm=UNDEFINED_RM) ;
    virtual void     CopyAddValuesFromEdge (Edge* elemX,
					    int indX, int indY, int icx, int icy,
					    CopyAddMode cam, InterpMode im,
					    NonConform nc,
					    ReceiveMode rm=UNDEFINED_RM) ;
    virtual void     CopyAddValuesFromFace (Face* elemX,
					    int indX, int indY, int icx, int icy,
					    CopyAddMode cam, InterpMode im,
					    NonConform nc,
					    ReceiveMode rm=UNDEFINED_RM) ;
    virtual void     CopyAddValuesFromVolume (Volume* elemX,
					      int indX, int indY, int icx, int icy,
					      CopyAddMode cam, InterpMode im,
					      NonConform nc, 
					      ReceiveMode rm=UNDEFINED_RM) ;

    // = EVALUATION OF USER-SUPPLIED FUNCTIONS
 
    void             SetTo (int index, int indexCoord, Function* f, int ic,
                            int icf, real t=ZERO) ;

    // = PRINTING AND DEBUGGING

    virtual Point*   GetBarycenter () = 0 ;
    virtual void     Print () = 0 ;
    virtual void     PrintPointers () ;
    void             PrintBarycenter() ;

    // = ACCESS TO THE DEGREES OF FREEDOM
    //    These methods are used only in conjunction with direct linear 
    //    solvers. They are irrelevant in the default case (iterative solvers).

    void             GetSetDof (int index, int iop, int idof, real& val) ;
    virtual void     GetSetInteriorDof (int index, int iop, int idof, 
                                        real& val) ;

    // = PARALLEL PROCESSING

    inline int       GetProcessorNumber () ;
    void             SetProcessorNumber (int p) ;

    // = MISCELLANEOUS

    Vector<ParentElement*>* GetCommonParentElements (Element* elemX, int indX);
    virtual real            GetTypicalLength (int dir) ;
  
} ;


//_____________________________ Element (inline) ______________________________


ElementaryField* Element :: GetElementaryField (int index)
{
    // Returns the receiver's elementary field whose index is 'index'.

# ifdef REQUIRE
    Require("'index'-th elementary field exists", HasElementaryField(index)) ;
    InMethod("Element::GetElementaryField(index)") ;
# endif

    return elementaryFields->At(index) ;
}

/*! \fn int Element :: GetProcessorNumber ()
 *  \brief Returns the rank of the processor assigned to the receiver.
 *  \return the rank (an int)
 */
int Element :: GetProcessorNumber ()
{
    return procNumber ;
}


/*! \fn boolean Element :: HasElementaryField (int index) ()
 *  \brief if the receiver already has an elementary field with index
 *  \return true' if the receiver already has an elementary field with index, false otherwise
 */
boolean Element :: HasElementaryField (int index)
{

# ifdef REQUIRE
    Require("'index' > 0", index > 0) ;
    InMethod("Element::HasElementaryField(index)") ;
# endif

    if (index<=elementaryFields->GetSize() && elementaryFields->At(index)!=NULL)
	return true ;
    else
	return false ;
}


/*! \fn boolean Element :: IsIncludedIn (Element* elem)
 *  \brief if the receiver already has an elementary field with index
 *  \param elem an Element
 *  \return 'true' if 'elem' contains the receiver. For example, if the receiver is a face, returns 'true' if 'elem' is the receiver itself or one of its edges. true' if the receiver already has an elementary field with index, false otherwise
 */
boolean Element :: IsIncludedIn (Element* elem)
{
    return (elem==this || elem->Has(this)) ;
}


#endif
