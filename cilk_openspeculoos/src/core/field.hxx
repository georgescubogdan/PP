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

#ifndef FIELD
#define FIELD

/*!
 * \file field.hxx
 * \brief Management of Fields
 * \author {Lots of authors}
 * \version 1.0
 */

#include "core/garbaged.hxx"
#include "core/vector.hxx"
#include "core/options.hxx"
#include "core/param.hxx"
#include <string>

class Mesh ;
class Element ;
class ParentElement ;
class FlatField ;
class Function ;
class FullMatrix ;

enum CopyAddMode      { COPY_CAM,  ADD_CAM } ;
enum InterpMode       { NORMAL_IM, TRANSPOSED_IM } ;
enum NonConform       { POINTWISE_NC, MORTAR_NC } ;
enum PolynomialChoice { POOR_PC,   RICH_PC } ;


/*! \class Field
   * \brief Abstract superclass of field classes.
   * 
   * This class mainly declares methods for polymorphic use for flat/
   * mortared fields.
   * 
   * A field can be given a name, to make printings more telling.
   */ 
class Field : public GarbagedObject
{

protected :

    string name ;

    static PolynomialChoice polynomialChoice ;         // a class attribute

public :

    // = CONSTRUCTORS

    Field () ;
    virtual ~Field () ; 

    // = DUPLICATION

    Field*             FieldDuplicate () ;
    Field*             FieldDuplicateEmpty (int ncomp=ALL) ;

    // = CHOICE OF THE MORTAR (RICH OR POOR)

    static void             SetPolynomialChoice (PolynomialChoice pc) ;
    static PolynomialChoice GetPolynomialChoice () ;

    // = WORK FIELDS

    Field*             FieldGetWork (int ncomp=ALL) ;
    virtual void       Retrieve (Field* w) = 0 ;

    // = INQUIRY AND ATTRIBUTE SETTINGS
  
    virtual FlatField* GetMain () = 0 ;
    virtual Mesh*      GetMesh () = 0 ;
    virtual int        GetNbComponents () = 0 ;
    virtual boolean    HasSameInterpolationAs (Field* x) = 0 ;
    virtual boolean    HasSameTypeAs (Field* x) = 0 ;
    virtual boolean    IsFlatField () ;
    virtual boolean    IsMortaredField () ;
    void               SetName (string s) ;

    // = MODIFYING THE INTERPOLATION
    //     Since these operations reinitialize the values to 0, they are 
    //     usually invoked just after the field has been created and before
    //     arithmetic operations are performed.

    virtual void    SetParentElements (ParentElement* p) = 0 ;
    virtual void    SetParentElements (ParentElement* p, int dir) = 0 ;
    virtual void    SetParentElement  (ParentElement* p, int dir, 
                                       Element* elem) = 0 ;
    virtual void    SetPolynomialDegree (int n) = 0 ;
    virtual void    SetPolynomialDegree (int n, int dir) = 0 ;
    virtual void    SetPolynomialDegree (int n, int dir, Element* elem) = 0 ;
    virtual void    ShiftPolynomialDegree (int shift) = 0 ;

    // = FAST ARITHMETICAL OPERATIONS

    virtual void    Add (Field* x) = 0 ;                      // summation
    virtual void    Add (Field* x, real alpha) = 0 ;
    virtual void    Add (Field* x, int icx, int icy) = 0 ;
    virtual void    Add (Field* x, real alpha, int icx, int icy) = 0 ;
    virtual void    Add (real val, int icy=ALL) = 0 ;

    inline  void    Subtract (Field* x) ;

    virtual void    Multiply (real alpha) = 0 ;               // multiplication
    virtual void    Multiply (Field* x) = 0 ;
    virtual void    Multiply (Field* x, int icx, int icy) = 0 ;
    virtual void    MultiplyAndAdd (real alpha, Field* x,real beta=ONE) = 0 ;

    virtual void    Inverse () = 0 ;                          // inverse

    virtual real    Dot (Field* x) = 0 ;                      // dot product
    virtual void    CopyFrom (Field* x) = 0 ;                 // copy
    virtual void    CopyFrom (Field* x, int icx, int icy) = 0 ;

    virtual real    GetEuclidianNorm () = 0 ;                 // norm
    virtual real    GetInfiniteNorm () = 0 ;

    inline void     SwitchSigns () ;

    // = VERY SLOW COPY/ADD OPERATIONS

    virtual void    CopyValuesFrom (FlatField* x, int icx, int icy) = 0 ;

    // = ADDITIONAL OPERATIONS

    virtual void    SetTo (Function* f, real t=ZERO) = 0 ;
    virtual void    SetTo (Function* f, int ic, int icf, real t=ZERO) = 0 ;
    virtual void    SetValues (real val) = 0 ;
    virtual void    SetValues (real val, int icy) = 0 ;

    // = C0 (AND MORTAR) CONTINUITY

    virtual void    Assemble () ;
    virtual void    Distribute () ;
    virtual void    EnforceMortarConstraints () ;

    // = ACCESS TO DEGREES OF FREEDOM (USED WITH DIRECT SOLVERS)

    void            CopyDofsToVector (RealVector* x) ;
    void            CopyDofsFromVector (RealVector* x) ;
    void            CopyDofsToMatrixColumn (FullMatrix* A, int j) ;
    inline real     GetDof (int idof) ;
    virtual int     GetNbDofs () = 0 ;
    virtual void    GetSetDof (int iop, int idof, real& val) = 0 ;
    inline void     SetDof (int idof, real val) ;

    // = DEBUGGING

    virtual void    Print () = 0 ;
    virtual void    Print (boolean showValues) = 0 ;
    virtual void    PrintComponent (int icy) = 0 ;

} ;


/*! \class FlatField
   * \brief Manages a field of values on a list of spectral elements.
   * 
   * A flat field is a list of values attached to the collocation points of
   * the spectral elements of a <{mesh}>. The values themselves are not
   * stored by the field but by the elementary fields of the spectral 
   * elements.
   * 
   * The field is characterized by <{nbComponents}>, e.g., 3 for a velocity
   * field in a 3D space, or 1 for a pressure field.
   * 
   * A flat field y has an <{index}>; every element uses this index to 
   * point, among its elementary fields, to the elementary field 
   * corresponding to y.
   * When a flat field is created, it gets its index from the static
   * attribute <{freeIndices>}, which is an array of available indices. 
   * Thanks to this array the index of the field is always chosen as small 
   * as possible, in order to maintain the table of elementary fields of 
   * every element as short as possible (this table is an array, rather 
   * than a list, in order to speed up access to elementary fields).
   * 
   * A few internal attributes (<{nbDofs}>, <{nbInteriorDofs}>) are defined
   * for efficiency purposes in direct-solver operations.
   * 
   * To a large extent, parallelism is handled in this class. This is done
   * at the price of a slight violation of encapsulation, because it relies
   * on some knowledge of the contents of methods of ElementaryField. But
   * it makes execution a bit faster, by avoiding sending a few messages.
   * 
   * Copies of fields are frequently needed. However, allocating/
   * deallocating a field is highly time consuming. Therefore, attributes
   * <{workFields}> and <{availability}> have been introduced.
   * <{workFields}> is a list of duplicates of the field; <{availability}>
   * states if the corresponding entry in <{workFields}> is available (not
   * currently used) or unavailable (currently used).
   * Instead of invoking method Duplicate, the user can therefore invoke
   * method GetWork. Then the field looks up for an available work field,
   * and returns it; dynamic allocation is thus avoided. If all work fields
   * are unavailable, it creates one.
   * To get rid of a work field w, the user invokes method Retrieve(w) -
   * not "delete w".
   * Attribute <{availability}> may seem useless; it is useful in debugging
   * phase, for tracing the work fields that the user forgot to retrieve.
   * 
   * Fields are the most memory consuming objects in Speculoos. Therefore a
   * record of the number of created fields is maintained, as a debugging
   * tool, via the static attributes <{nbFields0D}>, <{nbFields1D}>, etc.
   * These can be printed by invoking the class method <{PrintNbFields}>.
   */ 
class FlatField : public Field
{

protected :
  
    Mesh*               mesh ;
    int                 nbComponents ;
    int                 index ;
    Vector<FlatField*>* workFields ;
    Vector<boolean>*    availability ;
    FlatField*          owner ;
    int                 nbDofs ;                   // 2 ugly attributes for
    int                 nbInteriorDofs ;           // direct solvers

    static Vector<int>*     freeIndices ;          // 10 class attributes
    static int              nbFields0D ;
    static int              nbFields1D ;
    static int              nbFields2D ;
    static int              nbFields3D ;
    static int              nbComponents0D ;
    static int              nbComponents1D ;
    static int              nbComponents2D ;
    static int              nbComponents3D ;

    FlatField () ;                                 // protected methods
    void Count (FlatField* x) ;
    void Decount (FlatField* x) ;
    int  GetFreeIndex () ;
    
public :

    // = CONSTRUCTORS

    FlatField (int ncomp, Mesh* m, ParentElement* p) ;
    FlatField (int ncomp, Mesh* m, Vector<ParentElement*>* vp) ;
    FlatField (int ncomp, FlatField* x, Mesh* m) ;
    ~FlatField () ;

    // = INQUIRY AND ATTRIBUTE SETTINGS

    inline FlatField* GetMain () ;
    inline Mesh*      GetMesh () ;
    inline int        GetIndex () ;
    inline int        GetNbComponents () ;
    boolean           HasSameTypeAs (Field* x) ;
    boolean           HasSameInterpolationAs (Field* x) ;
    boolean           IsEqual (FlatField* x) ;
    boolean           IsFlatField () ;

    // = MODIFYING THE INTERPOLATION
    //     Since these operations reinitialize the values to 0, they are 
    //     usually invoked just after the field has been created and before
    //     arithmetic operations are performed.

    void        SetParentElements (ParentElement* p) ;
    void        SetParentElements (ParentElement* p, int dir) ;
    void        SetParentElement  (ParentElement* p, int dir, Element* elem) ;
    void        SetPolynomialDegree (int n) ;
    void        SetPolynomialDegree (int n, int dir) ;
    void        SetPolynomialDegree (int n, int dir, Element* elem) ;
    void        ShiftPolynomialDegree (int shift) ;

    // = DUPLICATION

    FlatField*  Duplicate () ;
    FlatField*  DuplicateEmpty (int ncomp=ALL) ;

    // = WORK FIELDS

    FlatField*  GetWork (int ncomp=ALL) ;
    FlatField*  GetOwner () ;
    void        Retrieve (Field* w) ;
    void        SetOwner (FlatField* x) ;

    // = FAST ARITHMETICAL OPERATIONS
    //     Most of these methods resort to the BLAS package.
    //     When two fields are involved, they must have the same mesh; also,
    //     the elements must have the same parent elements, because 
    //     interpolation is not provided.
    //     Unless specified, these methods involve all components.
    
    void        Add (Field* x) ;                       // summation
    void        Add (Field* x, real alpha) ;
    void        Add (Field* x, int icx, int icy) ;
    void        Add (Field* x, real alpha, int icx, int icy) ;
    void        Add (real val, int icy=ALL) ;

    void        Multiply (real alpha) ;                // multiplication
    void        Multiply (Field* x) ;
    void        Multiply (Field* x, int icx, int icy) ;
    void        MultiplyAndAdd (real alpha, Field* x,real beta=ONE) ;
    void        Inverse () ;                           // inverse
    real        Dot (Field* x) ;                    // dot product
    real        InteriorDot (Field* x) ;
    real        LocalDot (Field* x);                // dot product
    void        Sqrt () ;  
    void        Square () ;
    void        Power (real x) ; // MAH 04.2006
    real        GetEuclidianNorm () ;                  // norm
    real        GetInfiniteNorm () ;
    real        GetInteriorInfiniteNorm () ;
    real        GetL2Norm () ;
    real        GetH1Norm () ;
    FlatField*  GetNormByElement (string type) ;
    void        SetToNormByElement (FlatField* x, string type) ;
    void        CopyFrom (Field* x) ;                  // copy
    void        CopyFrom (Field* x, int icx, int icy) ;
    void        ScalarProductWithTensors(FlatField* x1,FlatField* x2) ;// MAH 07.12.2005
    void        SetToTrace(FlatField* x) ; // MAH 07.12.2005

    // = SLOW ARITHMETICAL OPERATIONS
    //     The 2 fields must be defined on the same mesh (and thus also have
    //     the same geometrical dimension).
    //     However, the elementary fields need not be defined on the same
    //     parent elements, because interpolation is provided.

    void        CopyInterpolateFrom  (FlatField* x) ;
    void        CopyInterpolateFrom  (FlatField* x, int icx, int icy) ;
    void        CopyInterpolateTFrom (FlatField* x) ;
    void        CopyInterpolateTFrom (FlatField* x, int icx, int icy) ;
    
    // = VERY SLOW COPY/ADD OPERATIONS
    //     Copies and/or adds values between fields defined on different
    //     geometrical support that may differ in geometrical dimension.
    //     Typical operations are copying the values of a field defined on the 
    //     border of the mesh into the unknown field.
    //     These operations allow also interpolation between spectral elements
    //     of different polynomial degrees.

    void        AddValuesTo    (FlatField* x, NonConform nc) ;
    void        CopyValuesFrom (FlatField* x, NonConform nc=MORTAR_NC) ;
    void        CopyValuesFrom (FlatField* x, int icx, int icy) ;
    void        CopyAddValuesFrom (FlatField* x, int icx, int icy, 
				   CopyAddMode cam, InterpMode im, NonConform nc) ;

    // = DIFFERENTIAL OPERATIONS
    void        SetToConvection (FlatField* x1, FlatField* x2);
    void        SetToDivergence (FlatField* x, boolean vector);
    void        SetToDyadicProduct (FlatField* x); 
    void        SetToGradient  (FlatField* x) ;
    void        SetToGradient  (FlatField* x, int icx, int icy, int dir) ;
    void        SetToGradientT (FlatField* x, int icx, int icy, int dir) ;
    void        SetToGradientOnParents  (FlatField* x,int icx,int icy,int dir);
    void        SetToGradientOnParentsT (FlatField* x,int icx,int icy,int dir);
    void        SetToLaplacian (FlatField* x) ;
    void        SetToRotational (FlatField* x) ;
    void        SetToRotRot (FlatField* x);
    void        SetToRotRot_old (FlatField* x);

    // = INTEGRAL OPERATIONS
    //     The differential operations included in these methods are performed
    //     element by element, not on the whole field. This is uglier, but
    //     saves a lot of memory space by avoiding the allocation of temporary
    //     copies of the whole field.

    void        SetToWeakIdentity   (FlatField* x, FlatField* rule) ;
    void        SetToWeakDivergence (FlatField* x, FlatField* rule) ;
    void        SetToWeakGradient   (FlatField* x, FlatField* rule) ;
    void        SetToWeakLaplacian  (FlatField* x, FlatField* rule) ;
    void        SetToWeakConvection (FlatField* x1, FlatField* x2,FlatField* rule) ;
    void        SetToWeakTensorDivergence (FlatField* x, FlatField* rule, boolean b) ; 


    // = SPECIAL INITIALIZATION (USER-DEFINED, COORDINATES)
    //     Field values can be initialized through evaluation of user-supplied
    //     functions. A field can also be initialized to one of the basic
    //     fields (coordinates, jacobian, inverse jacobian) of the mesh.

    void        SetTo (Function* f, real t=ZERO) ;
    void        SetTo (Function* f, int ic, int icf, real t=ZERO) ;
    void        SetToCoordinates () ;
    void        SetToJacobian () ;
    void        SetToInvJacobian () ;
    void        SetToNormal () ;

    // = ADDITIONAL OPERATIONS 

    void        DivideByWeights () ;
    void        MultiplyByWeights () ;
    void        OrientLike (RealVector* x) ;
    void        SetToAbsoluteValue () ;
    void        SetValues (real val) ;
    void        SetValues (real val, int icy) ;
    void        SetToZeroMean (int icy) ;
    void        SetToZeroAlgebraicMean (int icy) ;
    real        GetIntegral (int icy) ;
    real        GetMeanValue (int icy) ;
    real        SumValues () ;

    // = ACCESS TO DEGREES OF FREEDOM, ONE BY ONE
    //     Methods related to direct solvers (default solvers are iterative).
    //     Degrees of freedom are defined as field values that are not 
    //     constrained by inter element continuity conditions (strong or
    //     mortar).
    //     A point located on the boundary where Dirichlet b.c applies is 
    //     considered as a degree of freedom. Functionalities are given to 
    //     retrieve and set those degrees of freedom. This is expensive but 
    //     useful for interfacing tools like Lapack. For a flat field, all
    //     values are degrees of freedom.

    int         GetNbDofs () ;
    int         GetNbInteriorDofs () ;
    void        GetSetDof (int iop, int idof, real& val) ;
    void        GetSetInteriorDof (int iop, int idof, real& val) ;
  
    // = ACCESS TO DEGREES OF FREEDOM, ELEMENT BY ELEMENT
    //     Allows faster operations than dof-by-dof access.

    int         GetMaxNbLocalDofs () ;
    void        SetLocalDofs (int i, real val) ;
    void        SetLocalDofs (int i, FlatField* x) ;

    // = DE-ALLOCATING THE VALUES, FOR SAVING STORAGE

    void        DeallocateValues () ;
    void        ReallocateValues () ;

    // = DEBUGGING TOOLS

    void        Print () ;
    void        Print (boolean showValues) ;
    void        Print (Element* elem) ;
    void        PrintComponent (int icy) ;
    void        PrintWork () ;
    static void PrintNbInstances () ;
    static void PrintNbFields () ;
    void        PrintInProc0 (string fName, int icy); 
    void        ReadFromProc0 (string fName, int icy);
    void        Filter (real alpha, int shift); 
} ;

/*! \class MortaredField
   * \brief Manages fields submitted to continuity conditions.
   * 
   * A mortared field is:
   *  - a main field <{main}> defined on a mesh M, and
   *  - a mortar field <{mortar}> defined on the skeleton of M.
   * 
   * The mortar is used to enforce inter-element continuity constraints. 
   * At any interface between two elements, the field can be either
   * conformal, in which case usual C0 continuity applies, or non-
   * conformal, in which case mortar continuity (weak-sense continuity)
   * applies.
   * The mortar is typically used for performing an assembly (i.e., direct-
   * stiffness) operation.
   * 
   * Since the mortar field is also a MortaredField, the definition is
   * recursive until dimension 0. For example, a 2D mortared field is
   * made of 3 flat fields: 2D, 1D and 0D.
   * 
   * When an operation such as computing the norm of a mortared field or
   * computing the dot product of two mortared fields is performed, only
   * the degrees of freedom are taken into account. These are:
   * - the internal values of the main field, i.e., the values not located
   *   on the interface of the elements, and
   * - the internal values of the mortar, recursively.
   * For example, for a 2D mortared field, this amounts to:
   * - the internal values of the 2D main field,
   * - the internal values of the 1D mortar field, and
   * - the internal values of the 0D mortar field of the mortar field.
   * For a 0D field, i.e., for a field defined on a set of vertices, all 
   * values are considered as internal.
   * 
   * Element-by-element operations such as computing a weak gradient are
   * assigned to the main field <{main}>.
   * 
   * When the mortar field is created, at every interface with incompatible
   * interpolation the master side and the slave side must be chosen. If
   * the static attribute <{polynomialChoice}> of class Field is set to 
   * RICH_PC, then the choice is: the master side is the polynomially
   * richest side. The opposite choice is performed if <{polynomialChoice}>
   * = POOR_PC. By default, the rich choice is used, since it promotes
   * higher accuracy.
   */ 
class MortaredField : public Field
{

protected :

    FlatField*     main ;
    MortaredField* mortar ;

    MortaredField () ;
    void CreateMortar () ;
public :

    // = CONSTRUCTORS

    MortaredField (int ncomp, Mesh* m, ParentElement* p) ;
    MortaredField (int ncomp, Mesh* m, Vector<ParentElement*>* vp) ;
    MortaredField (int ncomp, FlatField* x) ;
    MortaredField (int ncomp, FlatField* x, Mesh* m) ;
    ~MortaredField () ;

    // = INQUIRY AND ATTRIBUTE SETTINGS

    inline FlatField*     GetMain () ;
    inline Mesh*          GetMesh () ;
    inline MortaredField* GetMortar () ;
    inline int            GetNbComponents () ;
    boolean               HasSameInterpolationAs (Field* x) ;
    boolean               HasSameTypeAs (Field* x) ;
    boolean               IsMortaredField () ;

    // = MODIFYING THE INTERPOLATION
    //     Since these operations reinitialize the values to 0, they are 
    //     usually invboked just after the field has been created and before
    //     arithmetic operations are performed.

    void    SetParentElements (ParentElement* p) ;
    void    SetParentElements (ParentElement* p, int dir) ;
    void    SetParentElement  (ParentElement* p, int dir, Element* elem);
    void    SetPolynomialDegree (int n) ;
    void    SetPolynomialDegree (int n, int dir) ;
    void    SetPolynomialDegree (int n, int dir, Element* elem) ;
    void    ShiftPolynomialDegree (int shift) ;

    // = DUPLICATION

    MortaredField*  Duplicate () ;
    MortaredField*  DuplicateEmpty (int ncomp=ALL) ;

    // = WORK FIELDS

    MortaredField*  GetWork (int ncomp=ALL) ;
    void            Retrieve (Field* w) ;

    // = BLAS OPERATIONS
    //     As opposed to class FlatField, operations are performed on the 
    //     cascade of fields through the mortar field.

    void    Add (Field* x) ;
    void    Add (Field* x, real alpha) ;
    void    Add (Field* x, int icx, int icy) ;
    void    Add (Field* x, real alpha, int icx, int icy) ;
    void    Add (real val, int icy=ALL) ;
    void    Multiply (real alpha) ;
    void    Multiply (Field* x) ;
    void    Multiply (Field* x, int icx, int icy) ;
    void    MultiplyAndAdd (real alpha, Field* x, real beta=ONE) ;
    void    Inverse () ;
    real    Dot (Field* x) ;
    real    GetEuclidianNorm () ;
    real    GetInfiniteNorm () ;
    void    CopyFrom (Field* x) ;
    void    CopyFrom (Field* x, int icx, int icy) ;
    void    ReadInterpolateFromFLUENT(real l_star, real u_star, int nb_vertex, int nb_quad, int interpoldegree);

    // = VERY SLOW COPY/ADD

    void    CopyValuesFrom (FlatField* x, int icx, int icy) ;

    // = ADDITIONAL OPERATIONS ON COEFFICIENTS

    void    SetValues (real val) ;
    void    SetValues (real val, int icy) ;
    void    SetTo (Function* f, real t=ZERO) ;
    void    SetTo (Function* f, int ic, int icf, real t=ZERO) ;
    void    SwitchSigns () ;
      
    // = SPECIFIC MORTAR FUNCTIONS
    //     These functions transfer information from top level field to fields
    //     of smaller dimension and vice versa.

    void    Assemble () ;
    void    Distribute () ;
    void    EnforceMortarConstraints () ;

    // = ADDITIONAL OPERATIONS ON VALUES

    void    GetSetDof (int iop, int idof, real& val) ;
    int     GetNbDofs () ;

    // = DEBUGGING TOOLS

    void    Print () ;
    void    Print (boolean showValues) ;
    void    PrintComponent (int icy) ;
    void    PrintWithMortar () ;
} ;


/*! \class ElementaryField
   * \brief Manages the contribution of one element to a flat field.
   * 
   * An elementary field stores the values of the contribution of a
   * spectral element to a flat field. As such, it is intended to be
   * attribute of one element.
   * 
   * Two elementary fields of a same spectral element can be defined on
   * different sets of parent elements; for example the parent elements for
   * a pressure field is typically 2-polynomial-order smaller the parent
   * elements for a velocity field. Therefore every elementary field has
   * its own set of <{parentElements}>, whose size is the dimension of the
   * spectral element (ie, 2 for a 2D element).
   * 
   * The values are stored by components in <{components}>, a vector of
   * size <{nbComponents}>. For example, an elementary velocity field in a
   * 3D problem has 3 components, whereas a pressure field has 1 component.
   * All components have the same number of values: the product of the
   * number of points of the parent elements.
   * Note that the number of components and the number of parent elements
   * are unrelated.
   * 
   * Within a given component, the order in which the values are stored in
   * the array is the following: the higher indices are blocked, then the
   * lower indices run.
   * For example, for an elementary field of a 2D element, the values are
   * stored in the following order (assuming 2 collocation points in the X
   * direction and 3 in the Y direction): val(1,1), val(2,1), val(1,2), val(2,2), val(1,3), val(2,3).
   * Idem in 3D: first all values for k=1, then all values for k=2, etc.
   * In other words the values are stored layer by layer, and then row by 
   * row.
   * 
   * When an elementary field is created, its number of components and a 
   * list of parent elements must be specified; therefore it is completely
   * defined immediately and all operations can be requested safely.
   * 
   * Temporary copies of elementary fields are frequently needed in 
   * efficiency-critical parts of the code. Therefore, in order to avoid
   * creating dynamically these copies, the class attribute <{workFields}>
   * has been introduced. It contains a number of elementary fields
   * with one component. Message eY->GetWork() takes an elementary field
   * from the list (and creates one if the list is empty) and adjusts its
   * interpolation (parent elements) to those of eY. Method <{Retrieve}>
   * releases it for future reuse.
   * An additional attribute <{isWork}> tells if the elementary field is a
   * work field or a normal field. It is used for debugging, so, in order 
   * to save memory space, it is defined only if the compilation option 
   * CHECK is turned on.
   */ 
class ElementaryField
{

protected :

    int                              nbComponents ;
    Vector<RealVector*>*             components ;
    Vector<ParentElement*>*          parentElements ;
#   ifdef CHECK
    boolean                          isWork ;
#   endif

    static Vector<ElementaryField*>* workFields ;

    real InteriorDot1D (ElementaryField* x) ;             // non-public methods
    real InteriorDot2D (ElementaryField* x) ;
    real InteriorDot3D (ElementaryField* x) ;
    real GetInteriorInfiniteNorm1D () ;
    real GetInteriorInfiniteNorm2D () ;
    real GetInteriorInfiniteNorm3D () ;
    static ElementaryField* CreateWork () ;

public :

    // = CONSTRUCTORS

    ElementaryField (int ncomp, boolean alloc, Vector<ParentElement*>* vp) ;
    ElementaryField (int ncomp, boolean alloc, ParentElement* p1) ;
    ElementaryField (int ncomp, boolean alloc, ParentElement* p1,
		     ParentElement* p2) ;
    ~ElementaryField () ;

    // = MANAGEMENT OF THE PARENT ELEMENTS

    boolean                        HasSameInterpolationAs (ElementaryField* x);
    inline Vector<ParentElement*>* GetParentElements () ;
    inline ParentElement*          GetParentElement (int) ;
    void                           SetParentElement (ParentElement* p,int dir);

    // = MANAGEMENT OF THE COMPONENTS

    void               AllocateValues () ;
    void               DeallocateValues () ;
    inline RealVector* GetComponent (int icy) ;
    inline int         GetNbComponents () ;
    inline boolean     HasZeroComponent (int icy) ;
    void               SetZeroComponent (int icy) ;
   

    //

    // = INQUIRY

    int                GetNbDofs () ;
    int                GetNbInteriorDofs () ;

    // = DUPLICATION

    ElementaryField*   Duplicate () ;
    ElementaryField*   DuplicateEmpty (int icy=ALL) ;

    // = FAST ARITHMETIC OPERATIONS
    //
    //     These operations are normally requested only from Field classes.
    //     They resort to the BLAS library whenever possible.
    //     When two fields are involved, they must be defined on the same mesh
    //     and have same interpolation.

    void         Add (real alpha, int icy) ;
    inline void  Add (ElementaryField* x, int icx, int icy) ;
    inline void  Add (ElementaryField* x, real alpha, int icx, int icy) ;

    inline void  Subtract (ElementaryField* x, int icx, int icy) ;

    void         Multiply (real alpha, int icy) ;
    inline void  Multiply (ElementaryField* x) ;
    inline void  Multiply (ElementaryField* x, int icx, int icy) ;
    void         MultiplyAndAdd (real alpha, ElementaryField* x,
                                 int icx, int icy, real beta=ONE) ;
    void         MultiplyAndSwitchSigns (ElementaryField* x, int icx, int icy) ;

    void         Divide (ElementaryField* x, int icx, int icy) ;

    void         Inverse () ;

    void         CopyFrom (ElementaryField* x, int icx, int icy) ;
    void         CopyTranspFrom (ElementaryField* x, int icx, int icy) ;

    real         Dot (ElementaryField* x) ;                    // dot product
    real         InteriorDot (ElementaryField* x) ;

    real         GetInfiniteNorm () ;                          // norm
    real         GetInteriorInfiniteNorm () ;
    real         SumValues () ;
    void         Sqrt () ;
    void         Square () ;
    void         Power (real x) ; // MAH 04.2006
    void         SetValues (real val) ;
    void         SetValues (real val, int icy) ;
    void         SetPolynomialDegree (int n) ;
    void         SetPolynomialDegree (int n, int dir) ;
    void         ShiftPolynomialDegree (int shift) ;
    void         ScalarProductWithTensors(ElementaryField* x1, ElementaryField* x2) ; // MAH 12.2005
    void         SetToTrace(ElementaryField *x) ; // MAH 12.2005

    // = SLOW ARITHMETIC OPERATIONS
    //    When 2 elementary fields are involved, interpolation is provided if
    //    they do not have the same parent elements. 

    void         CopyInterpolateFrom  (ElementaryField* x, int icx, int icy) ;
    void         CopyInterpolateTFrom (ElementaryField* x, int icx, int icy) ;

    // = DIFFERENTIAL OPERATIONS

    void         SetToGradient  (ElementaryField* x, int icx, int icy, int dir,ElementaryField* invJacobian) ;
    void         SetToGradientT (ElementaryField* x, int icx, int icy, int dir,ElementaryField* invJacobian) ;
    void         SetToGradientOnParents  (ElementaryField* x) ;
    void         SetToGradientOnParents  (ElementaryField* x, int icx, int icy,int dir) ;
    void         SetToGradientOnParentsT (ElementaryField* x, int icx, int icy,int dir) ;

    // = JACOBIAN AND INVERSE JACOBIAN

    void                    SetToJacobian    (ElementaryField* coord) ;
    void                    SetToInvJacobian (ElementaryField* coord,
                                              ElementaryField* jacobian) ;
    ElementaryField*        SetToNormal (ElementaryField* coord) ;

    // = WORK ELEMENTARY FIELDS

    ElementaryField*        GetWork () ;
    static ElementaryField* GetWork (ParentElement* p1, ParentElement* p2);
    static void             Retrieve (ElementaryField* work) ;

    // = OTHER OPERATIONS

    void                    DivideByWeights (int icy) ;
    void                    MultiplyByWeights (int icy) ;
    void                    OrientLike (RealVector* x) ;
    void                    SetToAbsoluteValue () ;

    // = PRINTINGS ON STANDARD OUTPUT

    void                    Print () ;
    void                    Print (boolean showValues) ;
    void                    PrintComponent (int icy) ;

} ;


//_______________________________ Field (inline) ______________________________


real Field :: GetDof (int idof) 
{
    // Returns the value of the 'idof'-th degrees of freedom of the receiver.

    real answer ; 

    GetSetDof(0,idof,answer) ; 
    return answer ; 
}


void Field :: SetDof (int idof, real val) 
{
    // Set the 'idof'-th degree of freedom of the receivr to 'val'.

    GetSetDof(1,idof,val) ; 
}


void Field :: Subtract (Field* x)
{
    // y = y - x .
    // Subtracts 'x' from the receiver 'y'.

    Add(x,ONE_MINUS) ;
}


void Field :: SwitchSigns ()
{
    // Replace every value of the receiver by its opposite.

    Multiply(ONE_MINUS) ;
}


//_____________________________ FlatField (inline) ____________________________


int FlatField :: GetIndex ()
{
    // Returns the index of the receiver.

    return index ;
}


FlatField* FlatField :: GetMain ()
{
    // Returns the receiver.

    return this ;
}


Mesh* FlatField :: GetMesh ()
{
    // Returns the mesh of the receiver.

    return mesh ;
}


int FlatField :: GetNbComponents () 
{
    // Returns the number of components of the receiver (typically 1, 2 or 3).

    return nbComponents ;
}


//____________________________ MortaredField (inline) _________________________


FlatField* MortaredField :: GetMain ()
{
    // Returns the main field of the receiver.

    return main ;
}


Mesh* MortaredField :: GetMesh ()
{
    // Returns the mesh of the main field of the receiver.

    return main->GetMesh() ;
}


MortaredField* MortaredField :: GetMortar ()
{
    // Returns the mortar field of the receiver.

    return mortar ;
}


int MortaredField :: GetNbComponents () 
{
    // Returns the number of components of the receiver.

    return main->GetNbComponents() ;
}


//___________________________ ElementaryField (inline) ________________________


void ElementaryField :: Add (ElementaryField* x, real alpha, int icx, int icy)
{
    // y(icy) = y(icy) + alpha*x(icx) .
    // Adds the 'icx'-th component of 'x', multiplied by 'alpha', to the 'icy'-th
    // component of the receiver 'y'.
    // All components are affected.
 
# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    InMethod("ElementaryField::Add(x,alpha,icx,icy)") ;
# endif
 
    components->At(icy)->Add(x->components->At(icx),alpha) ;
}


void ElementaryField :: Add (ElementaryField* x, int icx, int icy)
{
    // y(icy) = y(icy) + x(icx) .
    // Adds the 'icx'-th component of 'x' to the 'icy'-th component of the 
    // receiver 'y'.
    // All components are affected.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    InMethod("ElementaryField::Add(x,alpha,icx,icy)") ;
# endif

    components->At(icy)->Add(x->components->At(icx),ONE) ;
}


RealVector* ElementaryField :: GetComponent (int icy)
{
    // Returns the 'icy'-th component of the receiver.

# ifdef REQUIRE
    Require("correct processor", components != NULL) ;
    InMethod("ElementaryField::GetComponent(icy)") ;
# endif

    return components->At(icy) ;
}


int ElementaryField :: GetNbComponents () 
{
    // Returns the number of components of the receiver.

    return nbComponents ;
}


Vector<ParentElement*>* ElementaryField :: GetParentElements ()
{
    // Returns the vector of parent elements (one parent per space direction) of
    // the receiver.

    return parentElements ;
}


ParentElement* ElementaryField :: GetParentElement (int dir)
{
    // Returns the parent element in the 'dir'-th direction.

    return parentElements->At(dir) ;
}


boolean ElementaryField :: HasZeroComponent (int icy)
{
    // Returns 'true' if the 'icy'-th component of the receiver has all values 
    // equal to 0 (and consequently has not been stored), else returns 'false'.
    // See method ElementaryField::SetZeroComponent(icy).

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::HasZeroComponent(icy)") ;
# endif

    return components->At(icy) == NULL ;
}


void ElementaryField :: Multiply (ElementaryField* x)
{
    // y = y * x .
    // Multiplies every value of every component of the receiver 'y' by its
    // vis-a-vis in 'x'.

# ifdef REQUIRE
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components", nbComponents == x->nbComponents) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::Multiply(x)") ;
# endif

    int i ;

    for (i=1 ; i<=nbComponents ; i++)
	components->At(i)->Multiply(x->components->At(i)) ;
}


void ElementaryField :: Multiply (ElementaryField* x, int icx, int icy)
{
    // y = y * x .
    // Multiplies every value of the 'icy'-th component of the receiver 'y' by
    // its vis-a-vis in the 'icx'-th component of 'x'.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::Multiply(x,icx,icy)") ;
# endif

    components->At(icy)->Multiply(x->components->At(icx)) ;
}


void ElementaryField :: Subtract (ElementaryField* x, int icx, int icy)
{
    // y(icy) = y(icy) - x(icx) .
    // Subtracts the 'icx'-th component of 'x' to the 'icy'-th component of the 
    // receiver 'y'.
    // All components are affected.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    InMethod("ElementaryField::Subtract(x,alpha,icx,icy)") ;
# endif

    components->At(icy)->Add(x->components->At(icx),ONE_MINUS) ;
}


#endif

