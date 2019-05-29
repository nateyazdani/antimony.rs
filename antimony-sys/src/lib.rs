//! The API for the Antimony parser
//!
//! libAntimony uses a bison parser, libSBML, and internal C++ objects to read,
//! convert, store, and output abstracted models of biological systems.
//! Information about creating antimony-formatted input files is available from
//! [http://sys-bio.org/](http://sys-bio.org/). The functions described in this
//! document are a plain C API (application programming interface) to be used in
//! programs that want to convert Antimony models to their own internal formats,
//! and/or to convert Antimony models to and from SBML models.
//!
//! Note:  It is not currently possible to convert an internally-formatted model
//! into an Antimony model (the API has several `get` functions, but no `set`
//! functions). This restriction will be relaxed in future versions of the
//! library.
//!
//! # Model Structure
//!
//! Antimony models are modular; that is, one model may be defined and then used
//! as a subset of a larger model. As such, a single Antimony file may contain a
//! number of different models, each referred to as a 'module'. SBML core files
//! are not modular, and only contain a single model per file, but SBML L3
//! models with [Hierarchical Model Composition](http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/comp)
//! package are so translations from Antimony may either be 'flattened' to SBML
//! core, or may retain their modularity by using the hierarchy package.
//!
//! Converting files may be accomplished fairly straightforwardly using
//! `loadFile` to read in a file (of Antimony, SBML, or CellML format),
//! `getNumModules` and `getModuleNames` to obtain a list of the modules in that
//! file, and finally the `writeAntimonyFile`, `writeSBMLFile`,
//! `writeCompSBMLFile`, and `writeCellMLFile` routines to write them out again.
//!
//! Converting the models inside Antimony or SBML files to your own internal
//! format for your software is a bit more involved.  An example program is
//! provided that will read in a list of files, write out all the included
//! modules to individual Antimony and SBML files, and write out a list of the
//! components to the screen in a few different formats which will hopefully
//! prove useful as an example.  The procedure `printAllDataFor` may also be
//! used as an example, as it only uses routines provided in the API, and is
//! thus fully exportable and modifiable in other C or C++ programs.
//!
//! In general, Antimony models may contain:
//! - Species
//! - Reactions
//! - Reaction rates and other formulas
//! - Compartments
//! - Events
//! - Interactions (not modelled explicitly in SBML, only implied)
//! - DNA strands (not modelled explicitly in SBML)
//!
//! Note that there are several concepts modelled in SBML that are not modelled
//! in Antimony, units, annotations, and algebraic rules in particular. As such,
//! libAntimony is not the ideal tool to use to convert SBML models when those
//! elements are vital to your modelling; the excellent [libSBML](http://sbml.org/Software/libSBML)
//! should be used instead. Similarly, the concepts in Antimony not in SBML (the
//! main ones being interactions, and DNA strands) will be lost when converting
//! to an SBML file, so converting an Antimony file to SBML and back again may
//! be lossy. Do note that as of Antimony 2.1-beta, modularity may indeed be
//! preserved in SBML files using the new constructs of Hierarchical Model
//! Composition.
//!
//! # Symbol Kinds
//!
//! Many of the functions listed below ask for an enum value to determine what kind of symbol you are asking about.  This list is declared in enums.h, and is as follows:
//! - 0: `allSymbols`:         Every symbol of every type in Antimony
//! - 1: `allSpecies`:         All species, both const (boundary) and variable.
//! - 2: `allFormulas`:        All formulas (values defined by an equation), both const and variable.
//! - 3: `allDNA`:             All symbols defined to be DNA (operators and genes, but not strands).
//! - 4: `allOperators`:       All symbols defined to be operators (formulas embeddable in a DNA strand).
//! - 5: `allGenes`:           All symbols defined to be genes (reactions embeddable in a DNA strand).
//! - 6: `allReactions`:       All reactions (species being converted or created).
//! - 7: `allInteractions`:    All interactions (species involved in reaction rates).
//! - 8: `allEvents`:          All events.
//! - 9: `allCompartments`:    All compartments, both const and variable.
//! - 10: `allUnknown`:        All symbols whose type has never been defined or used.
//! - 11: `varSpecies`:        Variable species.
//! - 12: `varFormulas`:       Formulas (equations) that can change (including as a result of events).
//! - 13: `varOperators`:      Operators with variable values.
//! - 14: `varCompartments`:   Compartments with variable sizes.
//! - 15: `constSpecies`:      Constant species (aka 'border species').
//! - 16: `constFormulas`:     Formulas with constant values.
//! - 17: `constOperators`:    Operators with constant values.
//! - 18: `constCompartments`: Compartments with constant sizes.
//! - 19: `subModules`:        Submodules used within the current module.
//! - 20: `expandedStrands`:   DNA strands containing nothing but operators and genes--any sub-strands have been expanded to their component DNA objects, and those sub-strands are not included in any lists.
//! - 21: `modularStrands`:    All defined DNA strands, with some being subparts of the others.
//!
//! # Pointers and Memory Management
//!
//! The majority of the functions described below return pointers to arrays
//! and/or strings. These pointers you then own, and are created with `malloc`;
//! you must `free` them yourself to release the allocated memory. Some
//! programming environments will handle this automatically for you, and others
//! will not.  If you want to not bother with it, the function `freeAll` is
//! provided, which will free every pointer created by this library. In order
//! for this to work, however, you must have not freed a single provided pointer
//! yourself, and you must not subsequently try to reference any data provided
//! by the library (your own copies of the data will be fine, of course).
//!
//! If the library runs out of memory when trying to return a pointer, it will
//! return `NULL` instead and attempt to set an error message, retrievable with
//! `getLastError`.

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use std::os::raw::*;

pub const LIBANTIMONY_VERSION_STRING: &'static [u8; 7usize] = b"v2.7.0\0";

/// The different types of reactions and interactions.
///
/// Corresponds to `rd_type` in the C API.
#[repr(u32)]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Interaction {
    /// A reversible reaction: `->` or `<=>`
    ///
    /// Corresponds to `rdBecomes` in the C API.
    Becomes = 0,
    /// An activation interaction: `-o`
    ///
    /// Corresponds to `rdActivates` in the C API.
    Activates = 1,
    /// An inhibition interaction: `-|`
    ///
    /// Corresponds to `rdInhibits` in the C API.
    Inhibits = 2,
    /// A generic interaction: `-(`
    ///
    /// Corresponds to `rdInfluences` in the C API.
    Influences = 3,
    /// An irreversible reaction: `=>`
    ///
    /// Corresponds to `rdBecomesIrreversibly` in the C API.
    Transforms = 4,
}

/// Each kind classifies a different group of symbols, which may overlap (e.g., a single symbol
/// can be included in `allGenes` and also in `allReactions`).
///
/// Corresponds to `return_type` in the C API.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum SymbolKind {
    /// Every symbol of every type in Antimony
    ///
    /// Corresponds to `allSymbols` in the C API.
    Any = 0,
    /// All compartments, both const and variable.
    ///
    /// Corresponds to `allCompartments` in the C API.
    Compartment = 9,
    /// All submodel elements that have been deleted from the containing model.
    ///
    /// Corresponds to `allDeleted` in the C API.
    Deleted = 23,
    /// All symbols defined to be DNA (operators and genes, but not strands).
    ///
    /// Corresponds to `allDNA` in the C API.
    DNA = 3,
    /// All events.
    ///
    /// Corresponds to `allEvents` in the C API.
    Event = 8,
    /// All formulas (values defined by an equation), both const and variable.
    ///
    /// Corresponds to `allFormulas` in the C API.
    Formula = 2,
    /// All symbols defined to be genes (reactions embeddable in a DNA strand).
    ///
    /// Corresponds to `allGenes` in the C API.
    Gene = 5,
    /// All interactions (species involved in reaction rates).
    ///
    /// Corresponds to `allInteractions` in the C API.
    Interaction = 7,
    /// All symbols defined to be operators (formulas embeddable in a DNA strand).
    ///
    /// Corresponds to `allOperators` in the C API.
    Operator = 4,
    /// All reactions (species being converted or created).
    ///
    /// Corresponds to `allReactions` in the C API.
    Reaction = 6,
    /// All species, both const (border) and variable.
    ///
    /// Corresponds to `allSpecies` in the C API.
    Species = 1,
    /// All unit definitions.
    ///
    /// Corresponds to `allUnits` in the C API.
    Unit = 22,
    /// All symbols whose type has never been defined or used.
    ///
    /// Corresponds to `allUnknown` in the C API.
    Unknown = 10,
    /// Compartments with constant sizes.
    ///
    /// Corresponds to `constCompartments` in the C API.
    CompartmentConstant = 18,
    /// Formulas with constant values.
    ///
    /// Corresponds to `constFormulas` in the C API.
    FormulaConstant = 16,
    /// Operators with constant values.
    ///
    /// Corresponds to `constOperators` in the C API.
    OperatorConstant = 17,
    /// Constant species (aka 'border species').
    ///
    /// Corresponds to `constSpecies` in the C API.
    SpeciesConstant = 15,
    /// DNA strands containing nothing but operators and genes--any sub-strands have been expanded
    /// to their component DNA objects, and those sub-strands are not included in any lists.
    ///
    /// Corresponds to `expandedStrands` in the C API.
    StrandExpanded = 20,
    /// All defined DNA strands, with some being subparts of the others.
    ///
    /// Corresponds to `modularStrands` in the C API.
    StrandModular = 21,
    /// Submodules used within the current module.
    ///
    /// Corresponds to `subModules` in the C API.
    Module = 19,
    /// Compartments with variable sizes.
    ///
    /// Corresponds to `varCompartments` in the C API.
    CompartmentVariable = 14,
    /// Formulas (equations) that can change (including as a result of events).
    ///
    /// Corresponds to `varFormulas` in the C API.
    FormulaVariable = 12,
    /// Operators with variable values.
    ///
    /// Corresponds to `varOperators` in the C API.
    OperatorVariable = 13,
    /// Variable species.
    ///
    /// Corresponds to `varSpecies` in the C API.
    SpeciesVariable = 11,
}

/// Every symbol starts off with a default type based on its type (most things are
/// `formulaINITIAL`; reactions are `formulaKINETIC`, and events are `formulaTRIGGER`), and those
/// symbols that aren't reactions, events, or modules may be `formulaASSIGNMENT` (for those symbols
/// that have assignment rules) or `formulaRATE` (for those symbols that have rate rules).
/// `formulaASSIGNMENT` symbols have only one formula, but `formulaRATE` have two: one for their
/// initial condition, and one for how it changes with time.
///
/// Corresponds to `formula_type` in the C API.
#[repr(u32)]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum FormulaKind {
    /// Corresponds to `formulaINITIAL` in the C API.
    Initial = 0,
    /// Corresponds to `formulaASSIGNMENT` in the C API.
    Assignment = 1,
    /// Corresponds to `formulaRATE` in the C API.
    Rate = 2,
    /// Corresponds to `formulaKINETIC` in the C API.
    Kinetic = 3,
    /// Corresponds to `formulaTRIGGER` in the C API.
    Trigger = 4,
}

#[link(name = "antimony")]
extern "C" {
    /// Load a file of any format libAntimony knows about (potentially Antimony, SBML, or CellML).  If all attempts fail, the errors from the attempt to read the file in the Antimony format are saved, so if the file is actually SBML or CellML, the error is likely to be "but contains errors, the reported errors will be from the attempt to read it as Antimony, and a '-1' is returned.
    ///
    /// NOTE:  This function will not attempt to parse the file as SBML if libAntimony is compiled with the '-NSBML' flag, and will not attempt to parse the file as CellML if compiled with the '-NCELLML' flag.
    ///
    /// Returns a long integer indicating the index of the file read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadFile(filename: *const c_char) -> c_long;

    /// Load a string of any format libAntimony knows about (potentially Antimony, SBML, or CellML).  The first attempts to read the string as SBML, and if this results in an error, then reads it as Antimony.  If this, too, results in an error, the second error is saved, and a '-1' is returned.
    ///
    /// NOTE:  This function will not attempt to parse the string as SBML if libAntimony is compiled with the '-NSBML' flag, and will not attempt to parse the string as CellML if compiled with the '-NCELLML' flag.
    ///
    /// Returns a long integer indicating the index of the string read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadString(model: *const c_char) -> c_long;

    /// Loads a file and parses it as an Antimony file.  On an error, the error is saved, -1 is returned, and no information is stored.
    ///
    /// Returns a long integer indicating the index of the file read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadAntimonyFile(filename: *const c_char) -> c_long;

    /// Loads a string and parses it as an Antimony set of modules.  On an error, the error is saved, -1 is returned, and no information is stored.
    ///
    /// Returns a long integer indicating the index of the string read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadAntimonyString(model: *const c_char) -> c_long;

    /// Load a file known to be SBML.
    ///
    /// Loads a file and parses it (using libSBML) as an SBML file. On an error, the
    /// error is saved, -1 is returned, and no information is stored.
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NSBML` flag.
    ///
    /// Returns a long integer indicating the index of the file read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadSBMLFile(filename: *const c_char) -> c_long;

    /// Load a string known to be SBML.
    ///
    /// Loads a string and parses it (using libSBML) as an SBML file. On an error,
    /// the error is saved, -1 is returned, and no information is stored.
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NSBML` flag.
    ///
    /// Returns a long integer indicating the index of the string read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadSBMLString(model: *const c_char) -> c_long;

    /// Load a string known to be SBML with its file location.
    ///
    /// Loads a string and parses it (using libSBML) as an SBML file. On an error,
    /// the error is saved, -1 is returned, and no information is stored. This
    /// function additionally allows you to set the location of the string, in case
    /// there are relative file references in the file (as there can be in some
    /// hierarchical models).
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NSBML` flag.
    ///
    /// Returns a long integer indicating the index of the string read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadSBMLStringWithLocation(
        model: *const c_char,
        location: *const c_char,
    ) -> c_long;

    /// Load a file known to be CellML.
    ///
    /// Loads a file and parses it (using libCellML) as a CellML file. On an error,
    /// the error is saved, -1 is returned, and no information is stored.
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NCELLML` flag.
    ///
    /// Returns a long integer indicating the index of the file read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadCellMLFile(filename: *const c_char) -> c_long;

    /// Load a string known to be CellML.
    ///
    /// Loads a string and parses it (using libCellML) as a CellML file. On an error,
    /// the error is saved, -1 is returned, and no information is stored.
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NCELLML` flag.
    ///
    /// Returns a long integer indicating the index of the string read and stored.  On an error, returns -1 and no information is stored.
    ///
    /// See also `getLastError`.
    pub fn loadCellMLString(model: *const c_char) -> c_long;

    /// Returns the number of files loaded into memory so far.
    ///
    /// Every time `load{File,String}` is called successfully, the module(s) in
    /// those files are saved.  This function will tell you how many sets of modules
    /// from successful reads are resident in memory.
    ///
    /// Returns the number of files currently stored in memory.
    pub fn getNumFiles() -> c_ulong;

    /// Change the 'active' set of modules to the ones from the given index (as
    /// received from `load{File,String}`). Attempting to revert to a negative or
    /// nonexistent index returns `false` and the previous active set of modules is
    /// retained. A successful change return `true`.
    pub fn revertTo(index: c_long) -> bool;

    /// Clears memory of all files loaded. The next successful call to
    /// `load{File,String}` will return 0 as the first valid index.
    pub fn clearPreviousLoads();

    /// Add a directory in which imported files may be found, and in which to look
    /// for a `*.antimony` file (which contains rules about where to look locally
    /// for imported antimony and `*.sbml` files).
    pub fn addDirectory(directory: *const c_char);

    /// Clears the list of directories added with the `addDirectory` function.
    pub fn clearDirectories();

    /// Writes out an antimony-formatted file containing the given module. If no
    /// module name is given, all modules in the current set are returned. If the
    /// module depends on any sub-modules, those modules are written out as well,
    /// also in the antimony format.  Returns 0 on failure (and sets an error), 1 on
    /// success.
    pub fn writeAntimonyFile(
        filename: *const c_char,
        moduleName: *const c_char,
    ) -> c_int;

    /// Returns the same output as `writeAntimonyFile`, but to a string instead of a
    /// file.  Returns `NULL` on failure, and sets an error.
    pub fn getAntimonyString(
        moduleName: *const c_char,
    ) -> *mut c_char;

    /// Writes out a SBML-formatted XML file to the file indicated. The output is
    /// 'flattened', that is, all components of sub-modules are re-named and placed
    /// in a single model. Returns the output of libSBML's `writeSBML`, which
    /// "Returns non-zero on success and zero if the filename could not be opened
    /// for writing."  An error indicating this is set on returning zero.
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NSBML` flag.
    ///
    /// See also `getSBMLString`.
    pub fn writeSBMLFile(
        filename: *const c_char,
        moduleName: *const c_char,
    ) -> c_int;

    /// Returns the same output as `writeSBMLFile`, but to a string instead of a
    /// file. The output is 'flattened', that is, all components of sub-modules are
    /// re-named and placed in a single model. Returns the output of libSBML's
    /// `writeSBMLToString`, which "Returns the string on success and NULL if one of
    /// the underlying parser components fail (rare)."
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NSBML` flag.
    ///
    /// See also `writeSBMLFile`.
    pub fn getSBMLString(moduleName: *const c_char) -> *mut c_char;

    /// Writes out a CellML-formatted XML file to the file indicated, retaining the same Antimony hierarchy using the CellML 'component' hieararchy.  Returns one on success and zero on failure.
    /// NOTE:  This function is unavailable when libAntimony is compiled with the '-NCELLML' flag.
    ///
    /// See also `getCellMLString`.
    pub fn writeCellMLFile(
        filename: *const c_char,
        moduleName: *const c_char,
    ) -> c_int;

    /// Returns the same output as writeCellMLFile, but to a char* array instead of to a file.  Returns the string on success (as translated to 'char' from CellML's native 'wchar') and NULL on failure."
    /// NOTE:  This function is unavailable when libAntimony is compiled with the '-NCELLML' flag.
    ///
    /// See also `writeCellMLToString`.
    pub fn getCellMLString(
        moduleName: *const c_char,
    ) -> *mut c_char;

    /// An example function that will print to stdout all the information in the
    /// given module. This function probably isn't as useful to call as it is to
    /// examine and copy for your own purposes; it only calls procedures in the API.
    pub fn printAllDataFor(moduleName: *const c_char);

    /// Returns 'true' if the submitted module name exists in the current active set,
    /// 'false' if not.
    pub fn checkModule(moduleName: *const c_char) -> bool;

    /// When any function returns an error condition, a longer description of the
    /// problem is stored in memory, and is obtainable with this function. In most
    /// cases, this means that a call that returns a pointer returned `NULL` (or 0).
    pub fn getLastError() -> *mut c_char;

    /// When translating some other format to Antimony, elements that are unable to
    /// be translated are saved as warnings, retrievable with this function (returns
    /// `NULL` if no warnings present). Warnings may also be generated by problems
    /// discovered in `.antimony` files.
    pub fn getWarnings() -> *mut c_char;

    /// Returns the 'info' messages from libSBML. libAntimony always translates its
    /// modules into SBML to check for errors. If SBML finds errors, libAntimony
    /// gives up, passes on the error message, and does not save the model. However,
    /// libSBML may discover other things about your model it wants to tell you
    /// about, in 'info' and 'warning' messages. Info messages are just things it
    /// found it thinks you might want to know; warning messages are things it found
    /// which it feels violates 'best practices' in biological modelling, but not to
    /// the extent that it feels you did something actually wrong. Since Antimony is
    /// unitless, for example, you will always find warnings about how you didn't
    /// set any units. This function returns the 'info' messages from libSBML. If
    /// there are no info messages, returns an empty string.
    ///
    /// NOTE:  This function is unavailable when libAntimony is compiled with the
    /// `-NSBML` flag.
    ///
    /// See also `getSBMLWarnings`.
    pub fn getSBMLInfoMessages(
        moduleName: *const c_char,
    ) -> *mut c_char;

    /// Returns the 'warning' messages from libSBML.  If there are no warning
    /// messages (an unlikely occurrence), returns an empty string.
    ///
    /// NOTE: This function is unavailable when libAntimony is compiled with the
    /// `-NSBML` flag.
    ///
    /// See also `getSBMLInfoMessages`.
    pub fn getSBMLWarnings(
        moduleName: *const c_char,
    ) -> *mut c_char;

    /// Returns the number of modules in the current active set (the last file
    /// successfully loaded, or whichever file was returned to with `revertTo`).
    pub fn getNumModules() -> c_ulong;

    /// Returns an array of all the current module names.
    pub fn getModuleNames() -> *mut *mut c_char;

    /// Returns the nth module name. Returns `NULL` and sets an error if there is no
    /// such module `n`.
    pub fn getNthModuleName(n: c_ulong) -> *mut c_char;

    /// Returns the 'main' module name. In Antimony, this is either the module
    /// marked by an asterisk (`model *mainModel()`)  or the last module defined in
    /// the file. In translated SBML models, this is the model child of the `<sbml>`
    /// object. In translated CellML models, this is the 'containing' model that the
    /// translator creates to hold all the CellML components. Returns `NULL` only if
    /// there are no modules at all.
    pub fn getMainModuleName() -> *mut c_char;

    /// Returns the number of symbols defined to be in the interface of the given
    /// module. In other words, if a module is defined `module M(x, y, z)`, this
    /// returns 3. (Modules with no interface symbols return 0.)
    pub fn getNumSymbolsInInterfaceOf(
        moduleName: *const c_char,
    ) -> c_ulong;

    /// Returns the names of the symbols defined to be in the interface of the given
    /// module. In other words, if a module is defined `module M(x, y, z)`, this
    /// returns the list 'x, y, z'.  A module with no symbols defined in its
    /// interface would return a pointer to an empty string.
    pub fn getSymbolNamesInInterfaceOf(
        moduleName: *const c_char,
    ) -> *mut *mut c_char;

    /// Returns the Nth symbol name defined to be in the interface of the given
    /// module. If a module is defined `module M(x, y, z)`, calling this with n=0
    /// returns "x". If no such symbol is found, `NULL` is returned and an error is
    /// set.
    pub fn getNthSymbolNameInInterfaceOf(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the Nth replacement symbol name of a symbol that has replaced a different symbol in the given module, through the use of an 'is' construct, or through the use of a module's interface.
    /// See also `getNthFormerSymbolName`.
    /// See also `getNthReplacementSymbolName`.
    pub fn getNumReplacedSymbolNames(
        moduleName: *const c_char,
    ) -> c_ulong;

    /// Returns a list of pairs of symbol names that have been synchronized with each other--the first the symbol that was replaced, and the second the symbol used as the replacement.  These replacements are created when 'is' is used, and when a module's 'interface' (the symbols listed in parentheses) is used.
    /// See also `getNthFormerSymbolName`.
    /// See also `getNthReplacementSymbolName`.
    /// See also `getNthReplacementSymbolPair`.
    pub fn getAllReplacementSymbolPairs(
        moduleName: *const c_char,
    ) -> *mut *mut *mut c_char;

    /// Returns the Nth pair of symbol names that have been synchronized with each other--the first the symbol that was replaced, and the second the symbol used as the replacement.  These replacements are created when 'is' is used, and when a module's 'interface' (the symbols listed in parentheses) is used.
    /// See also `getNthFormerSymbolName`.
    /// See also `getNthReplacementSymbolName`.
    pub fn getNthReplacementSymbolPair(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns the Nth symbol name that has been replaced by a new symbol name in the given module, through the use of an 'is' construct, or through the use of a module's interface.
    /// See also `getNthReplacementSymbolName`.
    /// See also `GetNumReplacedSymbolNames`.
    pub fn getNthFormerSymbolName(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the Nth replacement symbol name of a symbol that has replaced a different symbol in the given module, through the use of an 'is' construct, or through the use of a module's interface.
    /// See also `getNthFormerSymbolName`.
    /// See also `GetNumReplacedSymbolNames`.
    pub fn getNthReplacementSymbolName(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the Nth replacement symbol name of a symbol that has replaced a different symbol in the given module, through the use of an 'is' construct, or through the use of a module's interface, between the given submodules, with the variable in the first submodule being the former variable name, and the variable in the second being the replacement variable name.  If an empty string is used as one of the submodule names, those synchronized variables that are not part of any submodule are searched for.
    /// See also `getNthFormerSymbolName`.
    /// See also `getNthReplacementSymbolName`.
    pub fn getNumReplacedSymbolNamesBetween(
        moduleName: *const c_char,
        formerSubmodName: *const c_char,
        replacementSubmodName: *const c_char,
    ) -> c_ulong;

    /// Returns a list of pairs of symbol names that have been synchronized with each other--the first the symbol that was replaced, and the second the symbol used as the replacement, between the given submodules, with the variable in the first submodule being the former variable name, and the variable in the second being the replacement variable name.  These replacements are created when 'is' is used, and when a module's 'interface' (the symbols listed in parentheses) is used.
    /// See also `getNthFormerSymbolName`.
    /// See also `getNthReplacementSymbolName`.
    /// See also `getNthReplacementSymbolPair`.
    pub fn getAllReplacementSymbolPairsBetween(
        moduleName: *const c_char,
        formerSubmodName: *const c_char,
        replacementSubmodName: *const c_char,
        n: c_ulong,
    ) -> *mut *mut *mut c_char;

    /// Returns the Nth pair of symbol names that have been synchronized with each other--the first the symbol that was replaced, and the second the symbol used as the replacement, between the given submodules, with the variable in the first submodule being the former variable name, and the variable in the second being the replacement variable name.  These replacements are created when 'is' is used, and when a module's 'interface' (the symbols listed in parentheses) is used.
    /// See also `getNthFormerSymbolName`.
    /// See also `getNthReplacementSymbolName`.
    pub fn getNthReplacementSymbolPairBetween(
        moduleName: *const c_char,
        formerSubmodName: *const c_char,
        replacementSubmodName: *const c_char,
        n: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns the Nth symbol name that has been replaced by a new symbol name in the given module, through the use of an 'is' construct, or through the use of a module's interface, between the given submodules, with the variable in the first submodule being the former variable name, and the variable in the second being the replacement variable name.
    /// See also `getNthReplacementSymbolName`.
    /// See also `GetNumReplacedSymbolNames`.
    pub fn getNthFormerSymbolNameBetween(
        moduleName: *const c_char,
        formerSubmodName: *const c_char,
        replacementSubmodName: *const c_char,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the Nth replacement symbol name of a symbol that has replaced a different symbol in the given module, through the use of an 'is' construct, or through the use of a module's interface, between the given submodules, with the variable in the first submodule being the former variable name, and the variable in the second being the replacement variable name.
    /// See also `getNthFormerSymbolName`.
    /// See also `GetNumReplacedSymbolNames`.
    pub fn getNthReplacementSymbolNameBetween(
        moduleName: *const c_char,
        formerSubmodName: *const c_char,
        replacementSubmodName: *const c_char,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the number of symbols of the given return type.  Useful when looping over the arrays in the subsequent functions.
    /// See also `get`.
    pub fn getNumSymbolsOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> c_ulong;

    /// Returns the names of the symbols of the given return type.  (In SBML, these are the 'id's.)
    pub fn getSymbolNamesOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> *mut *mut c_char;

    /// Returns the 'display names' of the symbols of the given return type.  (In SBML, these are the 'name's.)
    pub fn getSymbolDisplayNamesOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> *mut *mut c_char;

    /// Returns the equations associated with the symbols of the given return type.
    /// - Species:                 The initial assignment or assignment rule for the species in question
    /// - Formulas and operators:  The initial assignment or assignment rule for the formula in question
    /// - Compartments:            The initial assignment or assignment rule for the compartment in question
    /// - DNA elements:            The assignment rule or reaction rate of the element in question (no DNA element is defined by an initial assignment or by a rate rule with an initial assignment)
    /// - DNA Strands:             The assignment rule or reaction rate for the last element of the strand
    /// - Reactions and genes:     The reaction rate
    /// - Events:                  The trigger condition
    /// - Interactions:            Nothing
    /// - Modules:                 Nothing
    ///
    /// For elements that could have either initial assignments or assignment rules, use getTypeOfEquationForSymbol, or just use getSymbolInitialAssignmentsOfType and getSymbolAssignmentRulesOfType explicitly.
    pub fn getSymbolEquationsOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> *mut *mut c_char;

    /// Returns the equations associated with the initial assignment for symbols of the given return type.
    /// - Species:                 The initial assignment for the species in question
    /// - Formulas and operators:  The initial assignment of the formula in question
    /// - Compartments:            The initial assignment for the compartment
    ///
    /// - DNA Strands:             Nothing
    /// - Reactions and genes:     Nothing
    /// - Events:                  Nothing
    /// - Interactions:            Nothing
    /// - Modules:                 Nothing
    pub fn getSymbolInitialAssignmentsOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> *mut *mut c_char;

    /// Returns the equations associated with the assignment rule for symbols of the given return type.
    /// - Species:                 The assignment rule for the species in question
    /// - Formulas and operators:  The assignment rule of the formula in question
    /// - Compartments:            The assignment rule for the compartment
    /// - DNA Strands:             The assignment rule or reaction rate at the end of the strand.
    /// - Reactions and genes:     The reaction rate (for consistency with DNA strands)
    ///
    /// - Events:                  Nothing
    /// - Interactions:            Nothing
    /// - Modules:                 Nothing
    pub fn getSymbolAssignmentRulesOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> *mut *mut c_char;

    /// Returns the equations associated with the rate rule for symbols of the given return type.
    /// - Species:                 The rate rule for the species in question
    /// - Formulas and operators:  The rate rule of the formula in question
    /// - Compartments:            The rate rule for the compartment
    /// - DNA Strands:             The rate rule or reaction rate at the end of the strand.
    /// - Reactions and genes:     Nothing
    /// - Events:                  Nothing
    /// - Interactions:            Nothing
    /// - Modules:                 Nothing
    pub fn getSymbolRateRulesOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> *mut *mut c_char;

    /// Returns the compartments associated with the symbols of the given return type.  Note that unlike in SBML, any symbol of any type may have an associated compartment, including compartments themselves.  Rules about compartments in Antimony can be found in the <A class="el" target="_top" HREF="Tutorial.pdf">Tutorial.pdf</a> document included with this documentation.
    pub fn getSymbolCompartmentsOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
    ) -> *mut *mut c_char;

    /// Returns the name of the Nth symbol of the given type.  If no such symbol exists, NULL is returned and an error is set.  (In SBML, this is the 'id' of the element.)
    pub fn getNthSymbolNameOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the 'display name' of the Nth symbol of the given type.  If no such symbol exists, NULL is returned and an error is set.  (In SBML, this is the 'name' of the element.)
    pub fn getNthSymbolDisplayNameOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the equation associated with the Nth symbol of the given type.  If no equation is set for the symbol in question, an empty string is returned.  If no symbol can be found, NULL is returned and an error is set.
    pub fn getNthSymbolEquationOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the initial assignment associated with the Nth symbol of the given type.  If no initial assignment is set for the symbol in question, an empty string is returned.  If no symbol can be found, NULL is returned and an error is set.
    pub fn getNthSymbolInitialAssignmentOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the assignment rule associated with the Nth symbol of the given type.  If no assignment rule is set for the symbol in question, an empty string is returned.  If no symbol can be found, NULL is returned and an error is set.
    pub fn getNthSymbolAssignmentRuleOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the rate rule associated with the Nth symbol of the given type.  If no rate rule is set for the symbol in question, an empty string is returned.  If no symbol can be found, NULL is returned and an error is set.
    pub fn getNthSymbolRateRuleOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the name of the compartment associated with the nth symbol of the given type.  If no compartment is explicitly set in the file, the string "default_compartment" is returned.  If no symbol can be found, NULL is returned and an error is set.
    pub fn getNthSymbolCompartmentOfType(
        moduleName: *const c_char,
        rtype: SymbolKind,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the most specific return type available for the given symbolName.  A symbol defined to be a gene, for example, will return 'allGenes' and not 'allReactions', though the symbol does indeed qualify as a reaction.
    pub fn getTypeOfSymbol(
        moduleName: *const c_char,
        symbolName: *const c_char,
    ) -> SymbolKind;

    /// Returns the type of the 'main' equation associated with the given symbolName. All reactions
    /// will return `formulaKINETIC`, and all events will return `formulaTRIGGER`. All DNA elements
    /// that are not genes will return `formulaASSIGNMENT`, as DNA elements are defined by
    /// assignment rules and kinetic laws. All other symbols will return `formulaINITIAL` by
    /// default (i.e. in the case where no equation at all is associated with the symbol in
    /// question), and otherwise will return `formulaINITIAL` for symbols defined by initial
    /// assignments only, `formulaASSIGNMENT` for symbols defined by assignment rules, and
    /// `formulaRATE` for symbols defined by both initial assignments and rate rules (or just rate
    /// rules; it is valid though not simulatable to have a symbol with a rate rule but no initial
    /// assignment). In the case of rate rules, the initial assignment is found in the 'Equation'
    /// associated with the symbol, and the rate rule is found in the 'RateRule' associated with
    /// the symbol.
    pub fn getTypeOfEquationForSymbol(
        moduleName: *const c_char,
        symbolName: *const c_char,
    ) -> FormulaKind;

    /// Returns the name of the compartment the given symbol is a member of.  In antimony, all symbols may have compartments, not just species.  If a symbol has no set compartment, and is not a member of a symbol with a set compartment, this will return "default_compartment"
    pub fn getCompartmentForSymbol(
        moduleName: *const c_char,
        symbolName: *const c_char,
    ) -> *mut c_char;

    /// Returns the number of reactions (including genes) in the named module.  Useful when looping over all reactions in the arrays returned by subsequent functions.
    pub fn getNumReactions(moduleName: *const c_char) -> c_ulong;

    /// Returns the number of reactants (species on the left side of the reaction) for the given reaction.  If no such reaction is present, '0' is returned and an error is set.  Sadly, if there are no reactants, '0' is also returned, though no error is set.  So you'll have to keep track of this one on your own, most likely.
    pub fn getNumReactants(
        moduleName: *const c_char,
        rxn: c_ulong,
    ) -> c_ulong;

    /// Returns the number of products (species on the right side of the reaction) for the given reaction.  If no such reaction is present, '0' is returned and an error is set.  Sadly, if there are no products, '0' is also returned, though no error is set.  So you'll have to keep track of this one on your own, too.
    pub fn getNumProducts(
        moduleName: *const c_char,
        rxn: c_ulong,
    ) -> c_ulong;

    /// Returns all the reactant names for all reactions in the given module.  The dimensions of the included arrays can be found with 'getNumReactions' and 'getNumReactants' (the array is not 'square'--each sub array may have a different length).
    pub fn getReactantNames(
        moduleName: *const c_char,
    ) -> *mut *mut *mut c_char;

    /// Returns an array of all the reactant names for the given reaction.  The length of the array can be obtained with 'getNumReactants'.  If no such reaction is present, NULL is returned and an error is set.
    pub fn getNthReactionReactantNames(
        modulename: *const c_char,
        rxn: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns the mth reactant name of the mth reaction.  If no such reaction is present, NULL is returned and an error is set.
    pub fn getNthReactionMthReactantName(
        modulename: *const c_char,
        rxn: c_ulong,
        reactant: c_ulong,
    ) -> *mut c_char;

    /// Returns all the product names for all reactions in the given module.  The dimensions of the included arrays can be found with 'getNumReactions' and 'getNumProducts' (the array is not 'square'--each sub array may have a different length).
    pub fn getProductNames(
        moduleName: *const c_char,
    ) -> *mut *mut *mut c_char;

    /// Returns an array of all the product names for the given reaction.  The length of the array can be obtained with 'getNumProducts'.  If no such reaction is present, NULL is returned and an error is set.
    pub fn getNthReactionProductNames(
        modulename: *const c_char,
        rxn: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns the mth product name of the given reaction.  If no such reaction or product is present, NULL is returned and an error is set.
    pub fn getNthReactionMthProductName(
        modulename: *const c_char,
        rxn: c_ulong,
        product: c_ulong,
    ) -> *mut c_char;

    /// Returns a two-dimensional array of the stoichiometries for all reactants in all reactions in the given module.
    pub fn getReactantStoichiometries(moduleName: *const c_char) -> *mut *mut f64;

    /// Returns a two-dimensional array of the stoichiometries for all products in all reactions in the given module.
    pub fn getProductStoichiometries(moduleName: *const c_char) -> *mut *mut f64;

    /// Returns an array of the stoichiometries for the reactants of the Nth reaction in the module.  If no such reaction exists, an error is set and NULL is returned.
    pub fn getNthReactionReactantStoichiometries(
        moduleName: *const c_char,
        rxn: c_ulong,
    ) -> *mut f64;

    /// Returns an array of the stoichiometries for the products of the Nth reaction in the module.  If no such reaction exists, an error is set and NULL is returned.
    pub fn getNthReactionProductStoichiometries(
        moduleName: *const c_char,
        rxn: c_ulong,
    ) -> *mut f64;

    /// Returns the stoichiometry for the Mth reactant of the Nth reaction in the module.  If no such reactant or reaction exists, an error is set and 0 is returned.
    pub fn getNthReactionMthReactantStoichiometries(
        moduleName: *const c_char,
        rxn: c_ulong,
        reactant: c_ulong,
    ) -> f64;

    /// Returns the stoichiometries for the Mth product of the Nth reaction in the module.  If no such product or reaction exists, an error is set and 0 is returned.
    pub fn getNthReactionMthProductStoichiometries(
        moduleName: *const c_char,
        rxn: c_ulong,
        product: c_ulong,
    ) -> f64;

    /// Returns the number of interactions in the named module.  Useful when looping over all interactions in the arrays returned by subsequent functions.
    pub fn getNumInteractions(moduleName: *const c_char)
        -> c_ulong;

    /// Returns the number of interactors (species on the left side of the interaction) for the given interaction.  If no such interaction is present, '0' is returned and an error is set.  Sadly, if there are no interactors, '0' is also returned, though no error is set.  So you'll have to keep track of this one on your own, most likely.
    pub fn getNumInteractors(
        moduleName: *const c_char,
        rxn: c_ulong,
    ) -> c_ulong;

    /// Returns the number of interactees (reactions on the right side of the interaction) for the given interaction.  If no such interaction is present, '0' is returned and an error is set.  Sadly, if there are no interactees, '0' is also returned, though no error is set.  So you'll have to keep track of this one on your own, too.
    pub fn getNumInteractees(
        moduleName: *const c_char,
        rxn: c_ulong,
    ) -> c_ulong;

    /// Returns all the interactor names for all interactions in the given module.  The dimensions of the included arrays can be found with 'getNumInteractions' and 'getNumInteractors' (the array is not 'square'--each sub array may have a different length).
    pub fn getInteractorNames(
        moduleName: *const c_char,
    ) -> *mut *mut *mut c_char;

    /// Returns an array of all the interactor names for the given interaction.  The length of the array can be obtained with 'getNumInteractors'.  If no such interaction is present, NULL is returned and an error is set.
    pub fn getNthInteractionInteractorNames(
        modulename: *const c_char,
        rxn: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns the Mth interactor names for the given interaction.  If no such interactor or interaction is present, NULL is returned and an error is set.
    pub fn getNthInteractionMthInteractorName(
        modulename: *const c_char,
        interaction: c_ulong,
        interactor: c_ulong,
    ) -> *mut c_char;

    /// Returns all the interactee names for all interactions in the given module.  The dimensions of the included arrays can be found with 'getNumInteractions' and 'getNumInteractees' (the array is not 'square'--each sub array may have a different length).
    pub fn getInteracteeNames(
        moduleName: *const c_char,
    ) -> *mut *mut *mut c_char;

    /// Returns an array of all the interactee names for the given interaction.  The length of the array can be obtained with 'getNumInteractees'.  If no such interaction is present, NULL is returned and an error is set.
    pub fn getNthInteractionInteracteeNames(
        modulename: *const c_char,
        rxn: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns the Mth interactee name for the given interaction.  If no such interactee or interaction is present, NULL is returned and an error is set.
    pub fn getNthInteractionMthInteracteeName(
        modulename: *const c_char,
        interaction: c_ulong,
        interactee: c_ulong,
    ) -> *mut c_char;

    /// Returns an array of all the interaction dividers in the given module.  The length of the array can be obtained with 'getNumInteractions'.
    pub fn getInteractionDividers(moduleName: *const c_char) -> *mut Interaction;

    /// Returns the Nth interaction divider in the module.  If no such interaction is present, 0 is returned, which is 'rdBecomes, which is an invalid Interaction divider (since it's used for reactions instead).
    pub fn getNthInteractionDivider(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> Interaction;

    /// Returns an N x M stoichiometry matrix where N is the number of reactions in the model, and M is the number of variable species (or 'floating species').
    pub fn getStoichiometryMatrix(moduleName: *const c_char) -> *mut *mut f64;

    /// The row labels for the stoichiometry matrix.  Is exactly the same as calling 'getSymbolNamesOfType(moduleName, varSpecies)', but provided here so you don't have to think about it.
    pub fn getStoichiometryMatrixRowLabels(
        moduleName: *const c_char,
    ) -> *mut *mut c_char;

    /// The column labels for the stoichiometry matrix.  Is exactly the same as calling 'getSymbolNamesOfType(moduleName, allReactions)' but provided here so you don't have to think about it.
    pub fn getStoichiometryMatrixColumnLabels(
        moduleName: *const c_char,
    ) -> *mut *mut c_char;

    /// The number of rows in the stoichiometry matrix (or, the number of 'varSpecies').
    pub fn getStoichiometryMatrixNumRows(
        moduleName: *const c_char,
    ) -> c_ulong;

    /// The number of columns in the stoichiometry matrix (or, the number of 'allReactions').
    pub fn getStoichiometryMatrixNumColumns(
        moduleName: *const c_char,
    ) -> c_ulong;

    /// Returns the number of reactions (and hence reaction rates) in the module.  Useful for looping over all reaction rates in the following function.
    pub fn getNumReactionRates(
        moduleName: *const c_char,
    ) -> c_ulong;

    /// Returns an array of the reaction rates for the given module.  Is the same as 'getSymbolEquationsOfType(moduleName, allReactions)', but is provided for convenience.
    pub fn getReactionRates(
        moduleName: *const c_char,
    ) -> *mut *mut c_char;

    /// Returns the reaction rate for the Nth reaction in the module.  If the reaction exists, but its reaction rate has not been set, returns an empty string.  If the reaction does not exist, an error is set, and NULL is returned.
    pub fn getNthReactionRate(
        moduleName: *const c_char,
        rxn: c_ulong,
    ) -> *mut c_char;

    /// Returns the number of events in the given module.  Useful for subsequent functions that return arrays of information for all events.
    pub fn getNumEvents(moduleName: *const c_char) -> c_ulong;

    /// Returns the names of the events in the module.  Is the same as 'getSymbolNamesOfType(moduleName, allEvents)', but is provided for convenience.
    pub fn getEventNames(
        moduleName: *const c_char,
    ) -> *mut *mut c_char;

    /// Returns the name of the nth event in the module.
    pub fn getNthEventName(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> *mut c_char;

    /// Returns the number of assignments stored in the given event.  Useful when looping through those assignements in functions below.
    pub fn getNumAssignmentsForEvent(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> c_ulong;

    /// Returns the trigger for the given event, as an equation that can be interpreted in a boolean context.
    pub fn getTriggerForEvent(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> *mut c_char;

    /// Returns the delay for the given event, as an equation (if present; if the event has no
    /// delay, "" is returned.  If no such module or event is present, `NULL` is returned and an
    /// error is set.).
    pub fn getDelayForEvent(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> *mut c_char;

    /// Returns `true` if the given event has a delay; `false` otherwise.
    pub fn getEventHasDelay(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> bool;

    /// Returns the priority for the given event, as an equation (if present; if the event has no
    /// priority, "" is returned.  If no such module or event is present, NULL is returned and an
    /// error is set.).
    pub fn getPriorityForEvent(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> *mut c_char;

    /// Returns `true` if the given event has a priority; `false` otherwise.
    pub fn getEventHasPriority(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> bool;

    /// Returns the value of the persistence flag for the given event (default is `false`).
    /// Unable to return an error if there is no such event or module, so will simply return
    /// `false` in those situations, as well.
    pub fn getPersistenceForEvent(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> bool;

    /// Returns the value at time 0 for the given event trigger (default is `true`). Unable to
    /// return an error if there is no such event or module, so will simply return `true` in those
    /// situations, as well.
    pub fn getT0ForEvent(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> bool;

    /// Returns the value of the `fromTrigger` flag for the given event trigger (default is
    /// `true`). Unable to return an error if there is no such event or module, so will simply
    /// return `true` in those situations, as well.
    pub fn getFromTriggerForEvent(
        moduleName: *const c_char,
        event: c_ulong,
    ) -> bool;

    /// Each assignment for an event assigns a formula to a variable. This function returns the
    /// variable in question for the given event and assignment.
    pub fn getNthAssignmentVariableForEvent(
        moduleName: *const c_char,
        event: c_ulong,
        n: c_ulong,
    ) -> *mut c_char;

    /// Each assignment for an event assigns a formula to a variable. This function returns the in
    /// question in question for the given event and assignment.
    pub fn getNthAssignmentEquationForEvent(
        moduleName: *const c_char,
        event: c_ulong,
        n: c_ulong,
    ) -> *mut c_char;

    /// Returns the number of unique DNA strands in the module, as defined in the Antimony
    /// documentation (or, the number of physical cassettes of DNA present in the module). Useful
    /// in looping over the arrays returned by functions below.
    pub fn getNumDNAStrands(moduleName: *const c_char) -> c_ulong;

    /// Returns an array of DNA strand sizes for all strands in the module. Useful for looping over
    /// the arrays returned by `getDNAStrands`
    pub fn getDNAStrandSizes(
        moduleName: *const c_char,
    ) -> *mut c_ulong;

    /// Returns just the size (in number of components) of the nth DNA strand in the given module.
    /// If no such strand exists, sets an error and returns 0. This is actually useful here, since
    /// all DNA strands otherwise have a size of at least 1.
    pub fn getSizeOfNthDNAStrand(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> c_ulong;

    /// Returns an array of all DNA strands in the given module as lists of their components. All
    /// components are either Operator objects or Gene objects, depending on whether they have an
    /// associated reaction.
    pub fn getDNAStrands(
        moduleName: *const c_char,
    ) -> *mut *mut *mut c_char;

    /// Returns an array of names of the components in the nth DNA strand in the given module. If
    /// no such strand exists, sets an error and returns `NULL`.
    pub fn getNthDNAStrand(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns whether the given DNA strand was defined to be 'open' (that is, have an attachable
    /// end) at the upstream end (if 'upstream' is true) or at the downstream end (if 'upstream' is
    /// false).  This allows reproduction of a strand defined by "--X--Y--" vs. "X--Y", etc.
    pub fn getIsNthDNAStrandOpen(
        moduleName: *const c_char,
        n: c_ulong,
        upstream: bool,
    ) -> bool;

    /// Returns the sizes (in number of components) of all modular (separately-defined) DNA
    /// strands. Modular strands may contain genes, operators, and other DNA strands. Useful for
    /// looping over the strands in the array returned by `getModularDNAStrands`.
    pub fn getNumModularDNAStrands(
        moduleName: *const c_char,
    ) -> c_ulong;

    /// Returns an array of Modular DNA strand sizes for the given module. Useful for looping over
    /// the components in the sub-arrays returned by `getModularDNAStrands`.
    pub fn getModularDNAStrandSizes(
        moduleName: *const c_char,
    ) -> *mut c_ulong;

    /// Returns an array of strands, each of which has an array of the names of the components of
    /// that strand.  The components may be operators, genes, and other modular DNA strands.
    pub fn getModularDNAStrands(
        moduleName: *const c_char,
    ) -> *mut *mut *mut c_char;

    /// Returns an array of names of the components in the nth modular DNA strand in the given
    /// module. If no such strand exists, an error is set, and `NULL` is returned.
    pub fn getNthModularDNAStrand(
        moduleName: *const c_char,
        n: c_ulong,
    ) -> *mut *mut c_char;

    /// Returns whether the given modular DNA strand was defined to be 'open' (that is, have an
    /// attachable end) at the upstream end (if 'upstream' is true) or at the downstream end (if
    /// 'upstream' is false). This allows reproduction of a strand defined by "--X--Y--" vs.
    /// "X--Y", etc.
    pub fn getIsNthModularDNAStrandOpen(
        moduleName: *const c_char,
        n: c_ulong,
        upstream: bool,
    ) -> bool;

    /// Adds default initial values to the named module.
    /// By default, you must provide initial values to all the values in your model. If you call
    /// this function, all parameters and compartments will be given a default value of `1`, and
    /// all your species and reaction rates will be given a default value of `0`.
    ///
    /// Returns `true` if no such moduleName exists, `false` otherwise.
    pub fn addDefaultInitialValues(moduleName: *const c_char) -> bool;

    /// Sets whether bare numbers are dimensionless or undefined.
    /// By default, all numbers in mathematical equations do not have units unless they are
    /// explicitly declared (`"1 second"` vs. `"1"`). If this function is called with a value of
    /// `true`, all numbers without declared units will be assumed to have the units
    /// 'dimensionless'. If called with a value of `false`, the numbers will not have declared
    /// units (the default).
    pub fn setBareNumbersAreDimensionless(dimensionless: bool);
}
