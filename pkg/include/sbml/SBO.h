/**
 * @file    SBO.h
 * @brief   SBO utility Function
 * @author  Ben Bornstein
 *
 * $Id: SBO.h 12780 2011-02-08 04:12:54Z mhucka $
 * $HeadURL: http://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/SBO.h $
 *
 * <!--------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright (C) 2009-2011 jointly by the following organizations: 
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. EMBL European Bioinformatics Institute (EBML-EBI), Hinxton, UK
 *  
 * Copyright (C) 2006-2008 by the California Institute of Technology,
 *     Pasadena, CA, USA 
 *  
 * Copyright (C) 2002-2005 jointly by the following organizations: 
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. Japan Science and Technology Agency, Japan
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 * ---------------------------------------------------------------------- -->
 *
 * @class SBO
 * @brief Methods for interacting with Systems Biology Ontology terms.
 *
 * @htmlinclude not-sbml-warning.html
 *
 * The values of "id" attributes on SBML components allow the components to
 * be cross-referenced within a model. The values of "name" attributes on
 * SBML components provide the opportunity to assign them meaningful labels
 * suitable for display to humans.  The specific identifiers and labels
 * used in a model necessarily must be unrestricted by SBML, so that
 * software and users are free to pick whatever they need.  However, this
 * freedom makes it more difficult for software tools to determine, without
 * additional human intervention, the semantics of models more precisely
 * than the semantics provided by the SBML object classes defined in other
 * sections of this document.  For example, there is nothing inherent in a
 * parameter with identifier <code>k</code> that would indicate to a
 * software tool it is a first-order rate constant (if that's what
 * <code>k</code> happened to be in some given model).  However, one may
 * need to convert a model between different representations (e.g.,
 * Henri-Michaelis-Menten versus elementary steps), or to use it with
 * different modeling approaches (discrete or continuous).  One may also
 * need to relate the model components with other description formats such
 * as SBGN (<a target="_blank"
 * href="http://www.sbgn.org/">http://www.sbgn.org/</a>) using deeper
 * semantics.  Although an advanced software tool <em>might</em> be able to
 * deduce the semantics of some model components through detailed analysis
 * of the kinetic rate expressions and other parts of the model, this
 * quickly becomes infeasible for any but the simplest of models.
 *
 * An approach to solving this problem is to associate model components
 * with terms from carefully curated controlled vocabularies (CVs).  This
 * is the purpose of the optional "sboTerm" attribute provided on the SBML
 * class SBase.  The "sboTerm" attribute always refers to terms belonging
 * to the Systems Biology Ontology (SBO).
 * 
 * @section use Use of SBO
 * 
 * Labeling model components with terms from shared controlled vocabularies
 * allows a software tool to identify each component using identifiers that
 * are not tool-specific.  An example of where this is useful is the desire
 * by many software developers to provide users with meaningful names for
 * reaction rate equations.  Software tools with editing interfaces
 * frequently provide these names in menus or lists of choices for users.
 * However, without a standardized set of names or identifiers shared
 * between developers, a given software package cannot reliably interpret
 * the names or identifiers of reactions used in models written by other
 * tools.
 * 
 * The first solution that might come to mind is to stipulate that certain
 * common reactions always have the same name (e.g., "Michaelis-Menten"), but
 * this is simply impossible to do: not only do humans often disagree on
 * the names themselves, but it would not allow for correction of errors or
 * updates to the list of predefined names except by issuing new releases
 * of the SBML specification&mdash;to say nothing of many other limitations
 * with this approach.  Moreover, the parameters and variables that appear
 * in rate expressions also need to be identified in a way that software
 * tools can interpret mechanically, implying that the names of these
 * entities would also need to be regulated.
 * 
 * The Systems Biology Ontology (SBO) provides terms for identifying most
 * elements of SBML. The relationship implied by an "sboTerm" on an SBML
 * model component is <em>is-a</em> between the characteristic of the
 * component meant to be described by %SBO on this element and the %SBO
 * term identified by the value of the "sboTerm". By adding %SBO term
 * references on the components of a model, a software tool can provide
 * additional details using independent, shared vocabularies that can
 * enable <em>other</em> software tools to recognize precisely what the
 * component is meant to be.  Those tools can then act on that information.
 * For example, if the %SBO identifier @c "SBO:0000049" is assigned
 * to the concept of "first-order irreversible mass-action kinetics,
 * continuous framework", and a given KineticLaw object in a model has an
 * "sboTerm" attribute with this value, then regardless of the identifier
 * and name given to the reaction itself, a software tool could use this to
 * inform users that the reaction is a first-order irreversible mass-action
 * reaction.  This kind of reverse engineering of the meaning of reactions
 * in a model would be difficult to do otherwise, especially for more
 * complex reaction types.
 * 
 * The presence of %SBO labels on Compartment, Species, and Reaction
 * objects in SBML can help map those entities to equivalent concepts in
 * other standards, such as (but not limited to) BioPAX (<a target="_blank"
 * href="http://www.biopax.org/">http://www.biopax.org/</a>), PSI-MI (<a
 * target="_blank"
 * href="http://www.psidev.info/index.php?q=node/60">http://www.psidev.info</a>),
 * or the Systems Biology Graphical Notation (SBGN, <a target="_blank"
 * href="http://www.sbgn.org/">http://www.sbgn.org/</a>).  Such mappings
 * can be used in conversion procedures, or to build interfaces, with %SBO
 * becoming a kind of "glue" between standards of representation.
 * 
 * The presence of the label on a kinetic expression can also allow
 * software tools to make more intelligent decisions about reaction rate
 * expressions.  For example, an application could recognize certain types
 * of reaction formulas as being ones it knows how to solve with optimized
 * procedures.  The application could then use internal, optimized code
 * implementing the rate formula indexed by identifiers such as
 * @c "SBO:0000049" appearing in SBML models.
 * 
 * Finally, %SBO labels may be very valuable when it comes to model
 * integration, by helping identify interfaces, convert mathematical
 * expressions and parameters etc.
 * 
 * Although the use of %SBO can be beneficial, it is critical to keep in
 * mind that the presence of an "sboTerm" value on an object <em>must not
 * change the fundamental mathematical meaning</em> of the model.  An SBML
 * model must be defined such that it stands on its own and does not depend
 * on additional information added by %SBO terms for a correct mathematical
 * interpretation.  %SBO term definitions will not imply any alternative
 * mathematical semantics for any SBML object labeled with that term.  Two
 * important reasons motivate this principle.  First, it would be too
 * limiting to require all software tools to be able to understand the %SBO
 * vocabularies in addition to understanding SBML.  Supporting %SBO is not
 * only additional work for the software developer; for some kinds of
 * applications, it may not make sense.  If %SBO terms on a model are
 * optional, it follows that the SBML model <em>must</em> remain
 * unambiguous and fully interpretable without them, because an application
 * reading the model may ignore the terms.  Second, we believe allowing the
 * use of "sboTerm" to alter the mathematical meaning of a model would
 * allow too much leeway to shoehorn inconsistent concepts into SBML
 * objects, ultimately reducing the interoperability of the models.
 *
 * @section relationship Relationships between SBO and SBML
 *
 * The goal of %SBO labeling for SBML is to clarify to the fullest extent
 * possible the nature of each element in a model.  The approach taken in
 * %SBO begins with a hierarchically-structured set of controlled
 * vocabularies with six main divisions: (1) entity, (2) participant role,
 * (3) quantitative parameter, (4) modeling framework, (5) mathematical
 * expression, and (6) interaction.  The web site for %SBO (<a
 * target="_blank"
 * href="http://biomodels.net/sbo">http://biomodels.net</a>) should be
 * consulted for the current version of the ontology.
 * 
 * The Systems Biology Ontology (SBO) is not part of SBML; it is being
 * developed separately, to allow the modeling community to evolve the
 * ontology independently of SBML.  However, the terms in the ontology are
 * being designed keeping SBML components in mind, and are classified into
 * subsets that can be directly related with SBML components such as
 * reaction rate expressions, parameters, and others.  The use of "sboTerm"
 * attributes is optional, and the presence of "sboTerm" on an element does
 * not change the way the model is <em>interpreted</em>.  Annotating SBML
 * elements with %SBO terms adds additional semantic information that may
 * be used to <em>convert</em> the model into another model, or another
 * format.  Although %SBO support provides an important source of
 * information to understand the meaning of a model, software does not need
 * to support "sboTerm" to be considered SBML-compliant.
 *
 */

#ifndef SBML_h
#define SBML_h


#include <sbml/common/extern.h>
#include <sbml/SBMLNamespaces.h>


#ifdef __cplusplus


#include <string>

#include <cassert>
#include <algorithm>
#include <deque>
#include <map>


LIBSBML_CPP_NAMESPACE_BEGIN

class XMLAttributes;
class XMLOutputStream;
class SBMLErrorLog;

/* create a map of parent-child sbo terms */
typedef std::multimap<int, int>           ParentMap;
typedef ParentMap::iterator               ParentIter;
typedef std::pair<ParentIter, ParentIter> ParentRange;

static ParentMap mParent;


class LIBSBML_EXTERN SBO
{
public:
  /** @cond doxygen-libsbml-internal */

  /**
   * Reads (and checks) sboTerm from the given XMLAttributes set.
   *
   * @return the sboTerm as an integer or -1 if the sboTerm was not in the
   * correct format or not found.
   */
  static int readTerm (const XMLAttributes& attributes, SBMLErrorLog* log,
    unsigned int level = SBML_DEFAULT_LEVEL, 
    unsigned int version = SBML_DEFAULT_VERSION);

  /**
   * Writes sboTerm as an XMLAttribute to the given XMLOutputStream.
   */
  static void writeTerm (XMLOutputStream& stream, int sboTerm);

  /** @endcond */


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"quantiative parameter"</em>, @c false
   * otherwise.
   */
  static bool isQuantitativeParameter  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"participant role"</em>, @c false otherwise.
   */
  static bool isParticipantRole  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"modeling framework"</em>, @c false otherwise.
   */
  static bool isModellingFramework  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"mathematical expression"</em>, @c false otherwise.
   */
  static bool isMathematicalExpression  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"kinetic constant"</em>, @c false otherwise.
   */
  static bool isKineticConstant  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"reactant"</em>, @c false otherwise.
   */
  static bool isReactant  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"product"</em>, @c false otherwise.
   */
  static bool isProduct  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"modifier"</em>, @c false otherwise.
   */
  static bool isModifier  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"rate law"</em>, @c false otherwise.
   */
  static bool isRateLaw  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"event"</em>, @c false otherwise.
   */
  static bool isEvent  (unsigned int term);


  /**
    * Returns @c true if the given term identifier comes from the stated branch of %SBO.
    *
    * @return @c true if @p term is-a %SBO <em>"physical participant</em>, @c false otherwise.
    */
  static bool isPhysicalParticipant  (unsigned int term);


  /**
    * Returns @c true if the given term identifier comes from the stated branch of %SBO.
    *
    * @return @c true if @p term is-a %SBO <em>"participant"</em>, @c false otherwise.
    */
  static bool isParticipant  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"interaction"</em>, @c false otherwise.
   */
  static bool isInteraction  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"entity"</em>, @c false otherwise.
   */
  static bool isEntity  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"functional entity"</em>, @c false otherwise.
   */
  static bool isFunctionalEntity  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"material entity"</em>, @c false otherwise.
   */
  static bool isMaterialEntity  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"conservation law"</em>, @c false otherwise.
   */
  static bool isConservationLaw  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"steady state expression"</em>, @c false otherwise.
   */
  static bool isSteadyStateExpression  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"functional compartment"</em>, @c false otherwise.
   */
  static bool isFunctionalCompartment  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"continuous framework"</em>, @c false otherwise.
   */
  static bool isContinuousFramework  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"discrete framework"</em>, @c false otherwise.
   */
  static bool isDiscreteFramework  (unsigned int term);


  /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   *
   * @return @c true if @p term is-a %SBO <em>"logical framework"</em>, @c false otherwise.
   */
  static bool isLogicalFramework  (unsigned int term);


  /**
   * Predicate for checking whether the given term is obsolete.
   *
   * @return @c true if @p term is-a %SBO <em>"obsolete"</em> term, @c false otherwise.
   */
  static bool isObselete  (unsigned int term);


  /**
   * Returns the integer as a correctly formatted %SBO identifier string.
   *
   * @return the given integer sboTerm as a zero-padded seven digit string.
   *
   * @note If the sboTerm is not in the correct range
   * (0000000&ndash;9999999), an empty string is returned.
   */
  static std::string intToString (int sboTerm);


  /**
   * Returns the string as a correctly formatted %SBO integer portion.
   *
   * @return the given string sboTerm as an integer.  If the sboTerm is not
   * in the correct format (a zero-padded, seven digit string), <code>-1</code> is
   * returned.
   */
  static int stringToInt (const std::string& sboTerm);


  /**
   * Checks the format of the given %SBO identifier string.
   *
   * @return @c true if sboTerm is in the correct format (a zero-padded, seven
   * digit string), @c false otherwise.
   */
  static bool checkTerm (const std::string& sboTerm);


  /**
   * Checks the format of the given %SBO identifier, given in the form of
   * the integer portion alone.
   *
   * @return @c true if sboTerm is in the range (0000000&ndash;9999999), @c false
   * otherwise.
   */
  static bool checkTerm (int sboTerm);


protected:  
  /** @cond doxygen-libsbml-internal */

 /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   * * @return true if the term is-a parent, false otherwise
   */
  static bool isChildOf(unsigned int term, unsigned int parent);

 /**
   * Returns @c true if the given term identifier comes from the stated branch of %SBO.
   * populates the parent-child map
   */
  static void populateSBOTree();

  /** @endcond */
};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */


#ifndef SWIG


BEGIN_C_DECLS


END_C_DECLS


#endif  /* !SWIG */
#endif  /* SBML_h */
