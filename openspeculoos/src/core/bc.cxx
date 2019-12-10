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

// bc.cxx

#include "core/bc.hxx"

//___________________________ DirichletCondition ______________________________


DirichletCondition :: DirichletCondition (FlatField* p, int ic)
    : GarbagedObject ()
{

    // Constructor. Initialiazes the receiver to a Dirichlet condition on the
    // 'ic'-th component of the target field, with 'p' as prescribed values.

# ifdef REQUIRE
    Require("valid argument 'ic'", ic > 0) ;
    Require("the prescribed-value field has 1 component",
	    p->GetNbComponents() == 1) ;
    InMethod("DirichletCondition::DirichletCondition(p,ic)") ;
# endif

    ref(p) ;
    prescribedField  = p ;
    icomp            = ic ;
    function         = NULL ;
    homogeneousField = p->DuplicateEmpty(1) ;
}


DirichletCondition :: ~DirichletCondition () 
{
    // Destructor.

    unref(prescribedField) ;
    unref(function) ;
    unref(homogeneousField) ;
}


/**
* Copies into the 'icomp'-th component of 'x'
* - the prescribed value of the receiver (if the receiver has no time function), 
* - the values of the time function at time='t' (if the receiver has a time function).
*
* The copy is deep, i.e., if 'x' is a mortared field, the prescribed values are copied also into the mortar.
*/
void DirichletCondition :: ApplyOn (Field* x, real t)
{
# ifdef REQUIRE
    Require("x has an 'icomp'-th component", x->GetNbComponents() >= icomp) ;
    Require("x is a mortared field", x->IsMortaredField()) ;
    InMethod("DirichletCondition::ApplyOn(x)") ;
# endif

    if (x->GetMesh()->GetDimension()>prescribedField->GetMesh()->GetDimension())
	// since x is assumed to be assembled, it would be a waste of time to copy
	// the prescribed value on it if dim(x) > dim(prescribed)
	ApplyOn(((MortaredField*)x)->GetMortar(),t) ;

    else {
	if (function)
	    prescribedField->SetTo(function,1,icomp,t) ;
	x->CopyValuesFrom(prescribedField,1,icomp) ;       // recursive
    }
}


void DirichletCondition :: ApplyHomogeneousOn (Field* x)
{
    // Copies in a homogeneous way (i.e., with values equal to 0) the prescribed 
    // field of the receiver into the 'icomp'-th component of 'x'.
    //
    // The copy is deep, i.e., if 'x' is a mortared field, the prescribed values
    // are copied also into the mortar.

# ifdef REQUIRE
    Require("x has an 'icomp'-th component", x->GetNbComponents() >= icomp) ;
    InMethod("DirichletCondition::ApplyHomogeneousOn(x)") ;
# endif

    if (x->IsMortaredField()) {            // the usual case

	if (x->GetMesh()->GetDimension() > 
	    homogeneousField->GetMesh()->GetDimension())
	    // since x is assumed to be assembled, it would be a waste of time to 
	    // copy the prescribed value on it if dim(x) > dim(prescribed)
	    ApplyHomogeneousOn(((MortaredField*)x)->GetMortar()) ;

	else
	    x->CopyValuesFrom(homogeneousField,1,icomp) ;   // recursive
    }

    else                                   // a rare case
	x->CopyValuesFrom(homogeneousField,1,icomp) ;
}


void DirichletCondition :: ApplyHomogeneousOn (Field* x, int icx)
{
    // Copies in a homogeneous way (i.e., with values equal to 0) the prescribed 
    // field of the receiver into the 'icx'-th component of 'x'.
    //
    // The copy is deep, i.e., if 'x' is a mortared field, the prescribed values
    // are copied also into the mortar.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1  && icx <= x->GetNbComponents()) ;
    InMethod("DirichletCondition::ApplyHomogeneousOn(x,icx)") ;
# endif

    if (x->IsMortaredField()) {            // the usual case

	if (x->GetMesh()->GetDimension() > 
	    homogeneousField->GetMesh()->GetDimension())
	    // since x is assumed to be assembled, it would be a waste of time to 
	    // copy the prescribed value on it if dim(x) > dim(prescribed)
	    ApplyHomogeneousOn(((MortaredField*)x)->GetMortar(),icx) ;

	else
	    x->CopyValuesFrom(homogeneousField,1,icx) ;
    }

    else                                   // a rare case
	x->CopyValuesFrom(homogeneousField,1,icx) ;
}


DirichletCondition* DirichletCondition :: DuplicateHomogeneous ()
{
    // Returns a copy of the receiver, but with the values of the prescribed 
    // field set to 0.

    DirichletCondition *answer ;
    FlatField          *p ;

    p      = prescribedField->DuplicateEmpty() ;
    answer = new DirichletCondition(p,icomp) ;
    unref(p) ;

    return answer ;
}


void DirichletCondition :: Print ()
{
    // Prints the receiver on the standard output.

    printf("Dirichlet condition on component %d. Prescribed values are:\n",
	   icomp) ;
    prescribedField->Print() ;
}


void DirichletCondition :: SetFunction (Function* f)
{
    // Sets the time-scaling function of the receiver to 'f'.

    ref(f) ;
    unref(function) ;
    function = f ;
}

//_______________________________ NeumannCondition ____________________________


NeumannCondition :: NeumannCondition (FlatField* grad, int ic, real ext1,
                                      real ext2, real ext3)
    : GarbagedObject (),
      integrationRule(0)
{ 
    // Constructor. Initialiazes the receiver to a Neumann condition on the
    // 'ic'-th component of the target field, with 'grad' as prescribed values of
    // the gradient of that component.
    //
    // The boundary mesh (i.e., the mesh on which the Neumann condition applies)
    // is the mesh of 'grad'. This boundary mesh must be on the contour of 
    // 'interior', i.e., the elements of the boundary mesh must also belong to
    // the skin of 'interior'.
    // Supplying 'interior' is necessary: without it, the boundary mesh is unable
    // to provide the corrrct (i.e., external) orientation of the normal vector.

# ifdef REQUIRE
    Require("valid argument 'ic'", ic > 0) ;
    Require("correct number of components",
	    grad->GetNbComponents() == grad->GetMesh()->GetNbCoordinates()) ;
    InMethod("NeumannCondition::NeumannCondition(grad,ic,ext1,ext2,ext3)") ;
# endif

    ref(grad) ;

    g        = grad ;
    icomp    = ic ;
    function = NULL ;

    external = new RealVector(g->GetMesh()->GetNbCoordinates()) ;
    external->At(1) = ext1 ;
    external->At(2) = ext2 ;
    if (external->GetSize() == 3)
	external->At(3) = ext3 ;
}


NeumannCondition :: ~NeumannCondition () 
{
    // Destructor.

    unref(g) ;
    unref(function) ;
    delete external ;
}


void NeumannCondition :: ApplyOn (FlatField* x, real t)
{ 
    // Adds to the 'icomp'-th component of 'x' the following:
    //   integral on border {- g.n dS } 
    // where 'border' is the mesh of 'g', 'n' is the external normal vector at
    // every dof of 'border' and 'g' is evaluated at time 't'.
    //
    // If the receiver has no time-scaling function (attribute <{function}>),
    // then 't' is ignored and 'g' is used as is.

# ifdef REQUIRE
    Require("x has an 'icomp'-th component", x->GetNbComponents() >= icomp) ;
    Require("has an integration rule", integrationRule != NULL) ;
    Require("valid mesh dimensions", x->GetMesh()->GetDimension() ==
	    g->GetMesh()->GetDimension() + 1) ; 
    InMethod("NeumannCondition::ApplyOn(x)") ;
# endif

    FlatField *normal, *work, *gn ;
    int       i ;

    // 1. Get the external normal at every dof of the border

    normal = g->GetMesh()->GetNormal(external) ;

# ifdef CHECK
    if (normal->GetNbComponents() != g->GetNbComponents())
	Error("NeumannCondition::ApplyOn(x,t)","nbComp(n) != nbComp(g)") ;
# endif

    // 2. Get the scalar product -g.normal at every dof

    work = g->GetWork(1) ;
    gn   = g->GetWork(1) ;

    if (function){
	g->SetTo(function,t) ;
	Error(" NeumannCondition :: ApplyOn (FlatField* x, real t)"," Avec KIO pas de fct Neumann");
    }

    gn->SetValues(ZERO) ;
    for (i=1 ; i <= g->GetNbComponents() ; i++) {
	//    work->SetValues(ZERO);
	work->CopyFrom(g,i,1) ;                      // g(i)
	work->Multiply(normal,i,1) ;                 // g(i) n(i)
	gn->Subtract(work) ;                         // gn = sum { - g(i) n(i) }
    }

    // 3. Compute B g, the contour integral of (psi g.n)

    work->SetToWeakIdentity(gn,integrationRule) ;

    // 4. Add B g to x

    x->CopyAddValuesFrom(work,1,icomp,ADD_CAM,NORMAL_IM,MORTAR_NC) ;

    // 5. Clean up

    g->Retrieve(work) ;
    g->Retrieve(gn) ;
}


FlatField* NeumannCondition :: GetIntegrationRule ()
{ 
    // Returns the integration of the receiver.
    // If the receiver has no integration rule, it creates it as a copy of "g".

    if (integrationRule == NULL)
	integrationRule = g->DuplicateEmpty(1) ;

    return integrationRule ;
}


void NeumannCondition :: SetFunction (Function* f)
{ // Sets the time-scaling function of the receiver to 'f'.

    ref(f) ;
    unref(function) ;
    function = f ;
}


void NeumannCondition :: SetIntegrationRule (FlatField* rule)
{ // Sets the integration rule (for numerical quadratures) of the receiver to
  // 'rule'.

# ifdef REQUIRE
    Require("'rule' as same mesh as g", rule->GetMesh() == g->GetMesh()) ;
    Require("'rule' has 1 component", rule->GetNbComponents() == 1) ;
    InMethod("NeumannCondition::SetIntegrationRule(rule)") ;
# endif

    ref(rule) ;
    unref(integrationRule) ;

    integrationRule = rule ;
}


void NeumannCondition :: Print(real temps)
{
    // Prints the receiver on the standard output.

    printf("Neumann condition on component %d. Prescribed values are:\n",
	   icomp) ;
    if (function)
	g->SetTo(function,temps) ;

    g->Print() ;

}

