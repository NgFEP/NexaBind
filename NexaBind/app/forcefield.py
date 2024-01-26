"""
forcefield.py: Constructs NexaBind System objects based on a Topology and an XML force field description


"""
from __future__ import absolute_import, print_function


import os
import itertools
import xml.etree.ElementTree as etree
import math
import warnings
from math import sqrt, cos
from copy import deepcopy
from collections import defaultdict
import NexaBind.unit as unit
from . import element as elem
from NexaBind.app.internal.singleton import Singleton
#from NexaBind.app.internal import compiled
from NexaBind.app.internal.argtracker import ArgTracker

# Directories from which to load built in force fields.

_dataDirectories = None

def _getDataDirectories():
    global _dataDirectories
    if _dataDirectories is None:
        _dataDirectories = [os.path.join(os.path.dirname(__file__), 'data')]
        try:
            from importlib_metadata import entry_points
            for entry in entry_points().select(group='NexaBind.forcefielddir'):
                _dataDirectories.append(entry.load()())
        except:
            pass # importlib_metadata is not installed
    return _dataDirectories

def _convertParameterToNumber(param):
    if unit.is_quantity(param):
        if param.unit.is_compatible(unit.bar):
            return param / unit.bar
        return param.value_in_unit_system(unit.md_unit_system)
    return float(param)

def _parseFunctions(element):
    """Parse the attributes on an XML tag to find any tabulated functions it defines."""
    functions = []
    for function in element.findall('Function'):
        values = [float(x) for x in function.text.split()]
        if 'type' in function.attrib:
            functionType = function.attrib['type']
        else:
            functionType = 'Continuous1D'
        params = {}
        for key in function.attrib:
            if key.endswith('size'):
                params[key] = int(function.attrib[key])
            elif key.endswith('min') or key.endswith('max'):
                params[key] = float(function.attrib[key])
        if functionType.startswith('Continuous'):
            periodicStr = function.attrib.get('periodic', 'false').lower()
            if periodicStr in ['true', 'false', 'yes', 'no', '1', '0']:
                params['periodic'] = periodicStr in ['true', 'yes', '1']
            else:
                raise ValueError('ForceField: non-boolean value for periodic attribute in tabulated function definition')
        functions.append((function.attrib['name'], functionType, values, params))
    return functions

def _createFunctions(force, functions):
    """Add TabulatedFunctions to a Force based on the information that was recorded by _parseFunctions()."""
    print("The Forcefield._createFunctions() function is under construction")


# Enumerated values for nonbonded method

class NoCutoff(Singleton):
    def __repr__(self):
        return 'NoCutoff'
NoCutoff = NoCutoff()

class CutoffNonPeriodic(Singleton):
    def __repr__(self):
        return 'CutoffNonPeriodic'
CutoffNonPeriodic = CutoffNonPeriodic()

class CutoffPeriodic(Singleton):
    def __repr__(self):
        return 'CutoffPeriodic'
CutoffPeriodic = CutoffPeriodic()

class Ewald(Singleton):
    def __repr__(self):
        return 'Ewald'
Ewald = Ewald()

class PME(Singleton):
    def __repr__(self):
        return 'PME'
PME = PME()

class LJPME(Singleton):
    def __repr__(self):
        return 'LJPME'
LJPME = LJPME()

# Enumerated values for constraint type

class HBonds(Singleton):
    def __repr__(self):
        return 'HBonds'
HBonds = HBonds()

class AllBonds(Singleton):
    def __repr__(self):
        return 'AllBonds'
AllBonds = AllBonds()

class HAngles(Singleton):
    def __repr__(self):
        return 'HAngles'
HAngles = HAngles()

# A map of functions to parse elements of the XML file.

parsers = {}

class ForceField(object):
    """A ForceField constructs NexaBind System objects based on a Topology."""

    def __init__(self, *files):
        """Load one or more XML files and create a ForceField object based on them.

        Parameters
        ----------
        files : list
            A list of XML files defining the force field.  Each entry may
            be an absolute file path, a path relative to the current working
            directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a
            read() method from which the forcefield XML data can be loaded.
        """
        self._atomTypes = {}
        self._templates = {}
        self._patches = {}
        self._templatePatches = {}
        self._templateSignatures = {None:[]}
        self._atomClasses = {'':set()}
        self._forces = []
        self._scripts = []
        self._templateMatchers = []
        self._templateGenerators = []
        self.loadFile(files)

    def loadFile(self, files, resname_prefix=''):
        """Load an XML file and add the definitions from it to this ForceField.

        Parameters
        ----------
        files : string or file or tuple
            An XML file or tuple of XML files containing force field definitions.
            Each entry may be either an absolute file path, a path relative to the current working
            directory, a path relative to this module's data subdirectory (for
            built in force fields), or an open file-like object with a read()
            method from which the forcefield XML data can be loaded.
        prefix : string
            An optional string to be prepended to each residue name found in the
            loaded files.
        """

        if isinstance(files, tuple):
            files = list(files)
        else:
            files = [files]

        trees = []

        i = 0
        while i < len(files):
            file = files[i]
            tree = None
            try:
                # this handles either filenames or open file-like objects
                tree = etree.parse(file)
            except IOError:
                for dataDir in _getDataDirectories():
                    f = os.path.join(dataDir, file)
                    if os.path.isfile(f):
                        tree = etree.parse(f)
                        break
            except Exception as e:
                # Fail with an error message about which file could not be read.
                # TODO: Also handle case where fallback to 'data' directory encounters problems,
                # but this is much less worrisome because we control those files.
                msg  = str(e) + '\n'
                if hasattr(file, 'name'):
                    filename = file.name
                else:
                    filename = str(file)
                msg += "ForceField.loadFile() encountered an error reading file '%s'\n" % filename
                raise Exception(msg)
            if tree is None:
                raise ValueError('Could not locate file "%s"' % file)

            trees.append(tree)
            i += 1

            # Process includes in this file.

            if isinstance(file, str):
                parentDir = os.path.dirname(file)
            else:
                parentDir = ''
            for included in tree.getroot().findall('Include'):
                includeFile = included.attrib['file']
                joined = os.path.join(parentDir, includeFile)
                if os.path.isfile(joined):
                    includeFile = joined
                if includeFile not in files:
                    files.append(includeFile)

        # Load the atom types.

        for tree in trees:
            if tree.getroot().find('AtomTypes') is not None:
                for type in tree.getroot().find('AtomTypes').findall('Type'):
                    self.registerAtomType(type.attrib)

        # Load the residue templates.

        for tree in trees:
            if tree.getroot().find('Residues') is not None:
                for residue in tree.getroot().find('Residues').findall('Residue'):
                    resName = resname_prefix+residue.attrib['name']
                    template = ForceField._TemplateData(resName)
                    if 'override' in residue.attrib:
                        template.overrideLevel = int(residue.attrib['override'])
                    if 'rigidWater' in residue.attrib:
                        template.rigidWater = (residue.attrib['rigidWater'].lower() == 'true')
                    for key in residue.attrib:
                        template.attributes[key] = residue.attrib[key]
                    atomIndices = template.atomIndices
                    for ia, atom in enumerate(residue.findall('Atom')):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertParameterToNumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in atomIndices:
                            raise ValueError('Residue '+resName+' contains multiple atoms named '+atomName)
                        typeName = atom.attrib['type']
                        atomIndices[atomName] = ia
                        template.atoms.append(ForceField._TemplateAtomData(atomName, typeName, self._atomTypes[typeName].element, params))
                    for site in residue.findall('VirtualSite'):
                        template.virtualSites.append(ForceField._VirtualSiteData(site, atomIndices))
                    for bond in residue.findall('Bond'):
                        if 'atomName1' in bond.attrib:
                            template.addBondByName(bond.attrib['atomName1'], bond.attrib['atomName2'])
                        else:
                            template.addBond(int(bond.attrib['from']), int(bond.attrib['to']))
                    for bond in residue.findall('ExternalBond'):
                        if 'atomName' in bond.attrib:
                            template.addExternalBondByName(bond.attrib['atomName'])
                        else:
                            template.addExternalBond(int(bond.attrib['from']))
                    for patch in residue.findall('AllowPatch'):
                        patchName = patch.attrib['name']
                        if ':' in patchName:
                            colonIndex = patchName.find(':')
                            self.registerTemplatePatch(resName, patchName[:colonIndex], int(patchName[colonIndex+1:])-1)
                        else:
                            self.registerTemplatePatch(resName, patchName, 0)
                    self.registerResidueTemplate(template)

        # Load the patch definitions.

        for tree in trees:
            if tree.getroot().find('Patches') is not None:
                for patch in tree.getroot().find('Patches').findall('Patch'):
                    patchName = patch.attrib['name']
                    if 'residues' in patch.attrib:
                        numResidues = int(patch.attrib['residues'])
                    else:
                        numResidues = 1
                    patchData = ForceField._PatchData(patchName, numResidues)
                    for key in patch.attrib:
                        patchData.attributes[key] = patch.attrib[key]
                    for atom in patch.findall('AddAtom'):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertParameterToNumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = ForceField._PatchAtomData(atomName)
                        typeName = atom.attrib['type']
                        patchData.addedAtoms[atomDescription.residue].append(ForceField._TemplateAtomData(atomDescription.name, typeName, self._atomTypes[typeName].element, params))
                    for atom in patch.findall('ChangeAtom'):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertParameterToNumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = ForceField._PatchAtomData(atomName)
                        typeName = atom.attrib['type']
                        patchData.changedAtoms[atomDescription.residue].append(ForceField._TemplateAtomData(atomDescription.name, typeName, self._atomTypes[typeName].element, params))
                    for atom in patch.findall('RemoveAtom'):
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = ForceField._PatchAtomData(atomName)
                        patchData.deletedAtoms.append(atomDescription)
                    for bond in patch.findall('AddBond'):
                        atom1 = ForceField._PatchAtomData(bond.attrib['atomName1'])
                        atom2 = ForceField._PatchAtomData(bond.attrib['atomName2'])
                        patchData.addedBonds.append((atom1, atom2))
                    for bond in patch.findall('RemoveBond'):
                        atom1 = ForceField._PatchAtomData(bond.attrib['atomName1'])
                        atom2 = ForceField._PatchAtomData(bond.attrib['atomName2'])
                        patchData.deletedBonds.append((atom1, atom2))
                    for bond in patch.findall('AddExternalBond'):
                        atom = ForceField._PatchAtomData(bond.attrib['atomName'])
                        patchData.addedExternalBonds.append(atom)
                    for bond in patch.findall('RemoveExternalBond'):
                        atom = ForceField._PatchAtomData(bond.attrib['atomName'])
                        patchData.deletedExternalBonds.append(atom)
                    atomIndices = dict((atom.name, i) for i, atom in enumerate(patchData.addedAtoms[atomDescription.residue]+patchData.changedAtoms[atomDescription.residue]))
                    for site in patch.findall('VirtualSite'):
                        patchData.virtualSites[atomDescription.residue].append(ForceField._VirtualSiteData(site, atomIndices))
                    for residue in patch.findall('ApplyToResidue'):
                        name = residue.attrib['name']
                        if ':' in name:
                            colonIndex = name.find(':')
                            self.registerTemplatePatch(name[colonIndex+1:], patchName, int(name[:colonIndex])-1)
                        else:
                            self.registerTemplatePatch(name, patchName, 0)
                    self.registerPatch(patchData)

        # Load force definitions

        for tree in trees:
            for child in tree.getroot():
                if child.tag in parsers:
                    parsers[child.tag](child, self)

        # Load scripts

        for tree in trees:
            for node in tree.getroot().findall('Script'):
                self.registerScript(node.text)

        # Execute initialization scripts.

        for tree in trees:
            for node in tree.getroot().findall('InitializationScript'):
                exec(node.text, locals())

    def getGenerators(self):
        """Get the list of all registered generators."""
        return self._forces

    def registerGenerator(self, generator):
        """Register a new generator."""
        self._forces.append(generator)

    def registerAtomType(self, parameters):
        """Register a new atom type."""
        name = parameters['name']
        if name in self._atomTypes:
            raise ValueError('Found multiple definitions for atom type: '+name)
        atomClass = parameters['class']
        mass = _convertParameterToNumber(parameters['mass'])
        element = None
        if 'element' in parameters:
            element = parameters['element']
            if not isinstance(element, elem.Element):
                element = elem.get_by_symbol(element)
        self._atomTypes[name] = ForceField._AtomType(name, atomClass, mass, element)
        if atomClass in self._atomClasses:
            typeSet = self._atomClasses[atomClass]
        else:
            typeSet = set()
            self._atomClasses[atomClass] = typeSet
        typeSet.add(name)
        self._atomClasses[''].add(name)

    def registerResidueTemplate(self, template):
        """Register a new residue template."""
        if template.name in self._templates:
            # There is already a template with this name, so check the override levels.

            existingTemplate = self._templates[template.name]
            if template.overrideLevel < existingTemplate.overrideLevel:
                # The existing one takes precedence, so just return.
                return
            if template.overrideLevel > existingTemplate.overrideLevel:
                # We need to delete the existing template.
                del self._templates[template.name]
                existingSignature = _createResidueSignature([atom.element for atom in existingTemplate.atoms])
                self._templateSignatures[existingSignature].remove(existingTemplate)
            else:
                raise ValueError('Residue template %s with the same override level %d already exists.' % (template.name, template.overrideLevel))

        # Register the template.

        self._templates[template.name] = template
        signature = _createResidueSignature([atom.element for atom in template.atoms])
        if signature in self._templateSignatures:
            self._templateSignatures[signature].append(template)
        else:
            self._templateSignatures[signature] = [template]

    def registerPatch(self, patch):
        """Register a new patch that can be applied to templates."""
        self._patches[patch.name] = patch

    def registerTemplatePatch(self, residue, patch, patchResidueIndex):
        """Register that a particular patch can be used with a particular residue."""
        if residue not in self._templatePatches:
            self._templatePatches[residue] = set()
        self._templatePatches[residue].add((patch, patchResidueIndex))

    def registerScript(self, script):
        """Register a new script to be executed after building the System."""
        self._scripts.append(script)

    def registerTemplateMatcher(self, matcher):
        """Register an object that can override the default logic for matching templates to residues.

        A template matcher is a callable object that can be invoked as::

            template = f(forcefield, residue, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)

        where ``forcefield`` is the ForceField invoking it, ``residue`` is an NexaBind.app.Residue object,
        ``bondedToAtom[i]`` is the set of atoms bonded to atom index i, and ``ignoreExternalBonds`` and
        ``ignoreExtraParticles`` indicate whether external bonds and extra particules should be considered
        in matching templates.

        It should return a _TemplateData object that matches the residue.  Alternatively it may return
        None, in which case the standard logic will be used to find a template for the residue.

        .. CAUTION:: This method is experimental, and its API is subject to change.
        """
        self._templateMatchers.append(matcher)

    def registerTemplateGenerator(self, generator):
        """Register a residue template generator that can be used to parameterize residues that do not match existing forcefield templates.

        This functionality can be used to add handlers to parameterize small molecules or unnatural/modified residues.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        generator : function
            A function that will be called when a residue is encountered that does not match an existing forcefield template.

        When a residue without a template is encountered, the ``generator`` function is called with:

        ::
           success = generator(forcefield, residue)

        where ``forcefield`` is the calling ``ForceField`` object and ``residue`` is a NexaBind.app.topology.Residue object.

        ``generator`` must conform to the following API:

        ::
           generator API

           Parameters
           ----------
           forcefield : NexaBind.app.ForceField
               The ForceField object to which residue templates and/or parameters are to be added.
           residue : NexaBind.app.Topology.Residue
               The residue topology for which a template is to be generated.

           Returns
           -------
           success : bool
               If the generator is able to successfully parameterize the residue, `True` is returned.
               If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.

           The generator should either register a residue template directly with `forcefield.registerResidueTemplate(template)`
           or it should call `forcefield.loadFile(file)` to load residue definitions from an ffxml file.

           It can also use the `ForceField` programmatic API to add additional atom types (via `forcefield.registerAtomType(parameters)`)
           or additional parameters.

        """
        self._templateGenerators.append(generator)

    def _findAtomTypes(self, attrib, num):
        """Parse the attributes on an XML tag to find the set of atom types for each atom it involves.

        Parameters
        ----------
        attrib : dict of attributes
            The dictionary of attributes for an XML parameter tag.
        num : int
            The number of atom specifiers (e.g. 'class1' through 'class4') to extract.

        Returns
        -------
        types : list
            A list of atom types that match.

        """
        types = []
        for i in range(num):
            if num == 1:
                suffix = ''
            else:
                suffix = str(i+1)
            classAttrib = 'class'+suffix
            typeAttrib = 'type'+suffix
            if classAttrib in attrib:
                if typeAttrib in attrib:
                    raise ValueError('Specified both a type and a class for the same atom: '+str(attrib))
                if attrib[classAttrib] not in self._atomClasses:
                    types.append(None) # Unknown atom class
                else:
                    types.append(self._atomClasses[attrib[classAttrib]])
            elif typeAttrib in attrib:
                if attrib[typeAttrib] == '':
                    types.append(self._atomClasses[''])
                elif attrib[typeAttrib] not in self._atomTypes:
                    types.append(None) # Unknown atom type
                else:
                    types.append([attrib[typeAttrib]])
            else:
                types.append(None) # Unknown atom type
        return types

    def _parseTorsion(self, attrib):
        """Parse the node defining a torsion."""
        types = self._findAtomTypes(attrib, 4)
        if None in types:
            return None
        torsion = PeriodicTorsion(types)
        index = 1
        while 'phase%d'%index in attrib:
            torsion.periodicity.append(int(attrib['periodicity%d'%index]))
            torsion.phase.append(_convertParameterToNumber(attrib['phase%d'%index]))
            torsion.k.append(_convertParameterToNumber(attrib['k%d'%index]))
            index += 1
        return torsion

    class _SystemData(object):
        """Inner class used to encapsulate data about the system being created."""
        def __init__(self, topology):
            self.atomType = {}
            self.atomParameters = {}
            self.atomTemplateIndexes = {}
            self.atoms = list(topology.atoms())
            self.excludeAtomWith = [[] for _ in self.atoms]
            self.virtualSites = {}
            self.bonds = [ForceField._BondData(bond[0].index, bond[1].index) for bond in topology.bonds()]
            self.angles = []
            self.propers = []
            self.impropers = []
            self.atomBonds = [[] for _ in self.atoms]
            self.isAngleConstrained = []
            self.constraints = {}
            self.bondedToAtom = [set() for _ in self.atoms]

            # Record which atoms are bonded to each other atom

            for i in range(len(self.bonds)):
                bond = self.bonds[i]
                self.bondedToAtom[bond.atom1].add(bond.atom2)
                self.bondedToAtom[bond.atom2].add(bond.atom1)
                self.atomBonds[bond.atom1].append(i)
                self.atomBonds[bond.atom2].append(i)
            self.bondedToAtom = [sorted(b) for b in self.bondedToAtom]

        def addConstraint(self, system, atom1, atom2, distance):
            """Add a constraint to the system, avoiding duplicate constraints."""
            key = (min(atom1, atom2), max(atom1, atom2))
            if key in self.constraints:
                if self.constraints[key] != distance:
                    raise ValueError('Two constraints were specified between atoms %d and %d with different distances' % (atom1, atom2))
            else:
                self.constraints[key] = distance
                system.addConstraint(atom1, atom2, distance)

        def recordMatchedAtomParameters(self, residue, template, matches):
            """Record parameters for atoms based on having matched a residue to a template."""
            matchAtoms = dict(zip(matches, residue.atoms()))
            for atom, match in zip(residue.atoms(), matches):
                self.atomType[atom] = template.atoms[match].type
                self.atomParameters[atom] = template.atoms[match].parameters
                self.atomTemplateIndexes[atom] = match
                for site in template.virtualSites:
                    if match == site.index:
                        self.virtualSites[atom] = (site, [matchAtoms[i].index for i in site.atoms], matchAtoms[site.excludeWith].index)

    class _TemplateData(object):
        """Inner class used to encapsulate data about a residue template definition."""
        def __init__(self, name):
            self.name = name
            self.atoms = []
            self.atomIndices = {}
            self.virtualSites = []
            self.bonds = []
            self.externalBonds = []
            self.overrideLevel = 0
            self.rigidWater = True
            self.attributes = {}

        def getAtomIndexByName(self, atom_name):
            """Look up an atom index by atom name, providing a helpful error message if not found."""
            index = self.atomIndices.get(atom_name, None)
            if index is not None:
                return index

            # Provide a helpful error message if atom name not found.
            msg =  "Atom name '%s' not found in residue template '%s'.\n" % (atom_name, self.name)
            msg += "Possible atom names are: %s" % str(list(map(lambda x: x.name, self.atoms)))
            raise ValueError(msg)

        def addAtom(self, atom):
            self.atoms.append(atom)
            self.atomIndices[atom.name] = len(self.atoms)-1

        def addBond(self, atom1, atom2):
            """Add a bond between two atoms in a template given their indices in the template."""
            self.bonds.append((atom1, atom2))
            self.atoms[atom1].bondedTo.append(atom2)
            self.atoms[atom2].bondedTo.append(atom1)

        def addBondByName(self, atom1_name, atom2_name):
            """Add a bond between two atoms in a template given their atom names."""
            atom1 = self.getAtomIndexByName(atom1_name)
            atom2 = self.getAtomIndexByName(atom2_name)
            self.addBond(atom1, atom2)

        def addExternalBond(self, atom_index):
            """Designate that an atom in a residue template has an external bond, using atom index within template."""
            self.externalBonds.append(atom_index)
            self.atoms[atom_index].externalBonds += 1

        def addExternalBondByName(self, atom_name):
            """Designate that an atom in a residue template has an external bond, using atom name within template."""
            atom = self.getAtomIndexByName(atom_name)
            self.addExternalBond(atom)

        def areParametersIdentical(self, template2, matchingAtoms, matchingAtoms2):
            """Get whether this template and another one both assign identical atom types and parameters to all atoms.

            Parameters
            ----------
            template2: _TemplateData
                the template to compare this one to
            matchingAtoms: list
                the indices of atoms in this template that atoms of the residue are matched to
            matchingAtoms2: list
                the indices of atoms in template2 that atoms of the residue are matched to
            """
            atoms1 = [self.atoms[m] for m in matchingAtoms]
            atoms2 = [template2.atoms[m] for m in matchingAtoms2]
            if any(a1.type != a2.type or a1.parameters != a2.parameters for a1,a2 in zip(atoms1, atoms2)):
                return False
            # Properly comparing virtual sites really needs a much more complicated analysis.  This simple check
            # could easily fail for templates containing vsites, even if they're actually identical.  Since we
            # currently have no force fields that include both patches and vsites, I'm not going to worry about it now.
            if self.virtualSites != template2.virtualSites:
                return False
            return True

    class _TemplateAtomData(object):
        """Inner class used to encapsulate data about an atom in a residue template definition."""
        def __init__(self, name, type, element, parameters={}):
            self.name = name
            self.type = type
            self.element = element
            self.parameters = parameters
            self.bondedTo = []
            self.externalBonds = 0

    class _BondData(object):
        """Inner class used to encapsulate data about a bond."""
        def __init__(self, atom1, atom2):
            self.atom1 = atom1
            self.atom2 = atom2
            self.isConstrained = False
            self.length = 0.0

    class _VirtualSiteData(object):
        """Inner class used to encapsulate data about a virtual site."""
        def __init__(self, node, atomIndices):
            attrib = node.attrib
            self.type = attrib['type']
            if self.type == 'average2':
                numAtoms = 2
                self.weights = [float(attrib['weight1']), float(attrib['weight2'])]
            elif self.type == 'average3':
                numAtoms = 3
                self.weights = [float(attrib['weight1']), float(attrib['weight2']), float(attrib['weight3'])]
            elif self.type == 'outOfPlane':
                numAtoms = 3
                self.weights = [float(attrib['weight12']), float(attrib['weight13']), float(attrib['weightCross'])]
            elif self.type == 'localCoords':
                numAtoms = 0
                self.originWeights = []
                self.xWeights = []
                self.yWeights = []
                while ('wo%d' % (numAtoms+1)) in attrib:
                    numAtoms += 1
                    self.originWeights.append(float(attrib['wo%d' % numAtoms]))
                    self.xWeights.append(float(attrib['wx%d' % numAtoms]))
                    self.yWeights.append(float(attrib['wy%d' % numAtoms]))
                self.localPos = [float(attrib['p1']), float(attrib['p2']), float(attrib['p3'])]
            else:
                raise ValueError('Unknown virtual site type: %s' % self.type)
            if 'siteName' in attrib:
                self.index = atomIndices[attrib['siteName']]
                self.atoms = [atomIndices[attrib['atomName%d'%(i+1)]] for i in range(numAtoms)]
            else:
                self.index = int(attrib['index'])
                self.atoms = [int(attrib['atom%d'%(i+1)]) for i in range(numAtoms)]
            if 'excludeWith' in attrib:
                self.excludeWith = int(attrib['excludeWith'])
            else:
                self.excludeWith = self.atoms[0]

    class _PatchData(object):
        """Inner class used to encapsulate data about a patch definition."""
        def __init__(self, name, numResidues):
            self.name = name
            self.numResidues = numResidues
            self.addedAtoms = [[] for i in range(numResidues)]
            self.changedAtoms = [[] for i in range(numResidues)]
            self.deletedAtoms = []
            self.addedBonds = []
            self.deletedBonds = []
            self.addedExternalBonds = []
            self.deletedExternalBonds = []
            self.allAtomNames = set()
            self.virtualSites = [[] for i in range(numResidues)]
            self.attributes = {}

        def createPatchedTemplates(self, templates):
            """Apply this patch to a set of templates, creating new modified ones."""
            if len(templates) != self.numResidues:
                raise ValueError("Patch '%s' expected %d templates, received %d", (self.name, self.numResidues, len(templates)))

            # Construct a new version of each template.

            newTemplates = []
            for index, template in enumerate(templates):
                newTemplate = ForceField._TemplateData("%s-%s" % (template.name, self.name))
                newTemplates.append(newTemplate)

                # Build the list of atoms in it.

                for atom in template.atoms:
                    if not any(deleted.name == atom.name and deleted.residue == index for deleted in self.deletedAtoms):
                        newTemplate.addAtom(ForceField._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters))
                for atom in self.addedAtoms[index]:
                    if any(a.name == atom.name for a in newTemplate.atoms):
                        raise ValueError("Patch '%s' adds an atom with the same name as an existing atom: %s" % (self.name, atom.name))
                    newTemplate.addAtom(ForceField._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters))
                oldAtomIndex = dict([(atom.name, i) for i, atom in enumerate(template.atoms)])
                newAtomIndex = dict([(atom.name, i) for i, atom in enumerate(newTemplate.atoms)])
                for atom in self.changedAtoms[index]:
                    if atom.name not in newAtomIndex:
                        raise ValueError("Patch '%s' modifies nonexistent atom '%s' in template '%s'" % (self.name, atom.name, template.name))
                    newTemplate.atoms[newAtomIndex[atom.name]] = ForceField._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters)

                # Copy over the virtual sites, translating the atom indices.

                indexMap = dict([(oldAtomIndex[name], newAtomIndex[name]) for name in newAtomIndex if name in oldAtomIndex])
                for site in template.virtualSites:
                    if site.index in indexMap and all(i in indexMap for i in site.atoms):
                        newSite = deepcopy(site)
                        newSite.index = indexMap[site.index]
                        newSite.atoms = [indexMap[i] for i in site.atoms]
                        newTemplate.virtualSites.append(newSite)

                # Build the lists of bonds and external bonds.

                atomMap = dict([(template.atoms[i], indexMap[i]) for i in indexMap])
                deletedBonds = [(atom1.name, atom2.name) for atom1, atom2 in self.deletedBonds if atom1.residue == index and atom2.residue == index]
                for atom1, atom2 in template.bonds:
                    a1 = template.atoms[atom1]
                    a2 = template.atoms[atom2]
                    if a1 in atomMap and a2 in atomMap and (a1.name, a2.name) not in deletedBonds and (a2.name, a1.name) not in deletedBonds:
                        newTemplate.addBond(atomMap[a1], atomMap[a2])
                deletedExternalBonds = [atom.name for atom in self.deletedExternalBonds if atom.residue == index]
                for atom in template.externalBonds:
                    if template.atoms[atom].name not in deletedExternalBonds:
                        newTemplate.addExternalBond(indexMap[atom])
                for atom1, atom2 in self.addedBonds:
                    if atom1.residue == index and atom2.residue == index:
                        newTemplate.addBondByName(atom1.name, atom2.name)
                    elif atom1.residue == index:
                        newTemplate.addExternalBondByName(atom1.name)
                    elif atom2.residue == index:
                        newTemplate.addExternalBondByName(atom2.name)
                for atom in self.addedExternalBonds:
                    newTemplate.addExternalBondByName(atom.name)

                # Add new virtual sites.

                indexMap = dict((i, newAtomIndex[atom.name]) for i, atom in enumerate(self.addedAtoms[index]+self.changedAtoms[index]))
                for site in self.virtualSites[index]:
                    newSite = deepcopy(site)
                    newSite.index = indexMap[site.index]
                    newSite.atoms = [indexMap[i] for i in site.atoms]
                    newSite.excludeWith = indexMap[site.excludeWith]
                    newTemplate.virtualSites = [site for site in newTemplate.virtualSites if site.index != newSite.index]
                    newTemplate.virtualSites.append(newSite)
            return newTemplates

    class _PatchAtomData(object):
        """Inner class used to encapsulate data about an atom in a patch definition."""
        def __init__(self, description):
            if ':' in description:
                colonIndex = description.find(':')
                self.residue = int(description[:colonIndex])-1
                self.name = description[colonIndex+1:]
            else:
                self.residue = 0
                self.name = description

    class _AtomType(object):
        """Inner class used to record atom types and associated properties."""
        def __init__(self, name, atomClass, mass, element):
            self.name = name
            self.atomClass = atomClass
            self.mass = mass
            self.element = element

    class _AtomTypeParameters(object):
        """Inner class used to record parameter values for atom types."""
        def __init__(self, forcefield, forceName, atomTag, paramNames):
            self.ff = forcefield
            self.forceName = forceName
            self.atomTag = atomTag
            self.paramNames = paramNames
            self.paramsForType = {}
            self.extraParamsForType = {}

        def registerAtom(self, parameters, expectedParams=None):
            if expectedParams is None:
                expectedParams = self.paramNames
            types = self.ff._findAtomTypes(parameters, 1)
            if None not in types:
                values = {}
                extraValues = {}
                for key in parameters:
                    if key in expectedParams:
                        values[key] = _convertParameterToNumber(parameters[key])
                    else:
                        extraValues[key] = parameters[key]
                if len(values) < len(expectedParams):
                    for key in expectedParams:
                        if key not in values:
                            raise ValueError('%s: No value specified for "%s"' % (self.forceName, key))
                for t in types[0]:
                    self.paramsForType[t] = values
                    self.extraParamsForType[t] = extraValues

        def parseDefinitions(self, element):
            """"Load the definitions from an XML element."""
            expectedParams = list(self.paramNames)
            excludedParams = [node.attrib['name'] for node in element.findall('UseAttributeFromResidue')]
            for param in excludedParams:
                if param not in expectedParams:
                    raise ValueError('%s: <UseAttributeFromResidue> specified an invalid attribute: %s' % (self.forceName, param))
                expectedParams.remove(param)
            for atom in element.findall(self.atomTag):
                for param in excludedParams:
                    if param in atom.attrib:
                        raise ValueError('%s: The attribute "%s" appeared in both <%s> and <UseAttributeFromResidue> tags' % (self.forceName, param, self.atomTag))
                self.registerAtom(atom.attrib, expectedParams)

        def getAtomParameters(self, atom, data):
            """Get the parameter values for a particular atom."""
            t = data.atomType[atom]
            p = data.atomParameters[atom]
            if t in self.paramsForType:
                values = self.paramsForType[t]
                result = [None]*len(self.paramNames)
                for i, name in enumerate(self.paramNames):
                    if name in values:
                        result[i] = values[name]
                    elif name in p:
                        result[i] = p[name]
                    else:
                        raise ValueError('%s: No value specified for "%s"' % (self.forceName, name))
                return result
            else:
                raise ValueError('%s: No parameters defined for atom type %s' % (self.forceName, t))

        def getExtraParameters(self, atom, data):
            """Get extra parameter values for an atom that appeared in the <Atom> tag but were not included in paramNames."""
            t = data.atomType[atom]
            if t in self.paramsForType:
                return self.extraParamsForType[t]
            else:
                raise ValueError('%s: No parameters defined for atom type %s' % (self.forceName, t))


    def _getResidueTemplateMatches(self, res, bondedToAtom, templateSignatures=None, ignoreExternalBonds=False, ignoreExtraParticles=False):
        """Return the templates that match a residue, or None if none are found.

        Parameters
        ----------
        res : Topology.Residue
            The residue for which template matches are to be retrieved.
        bondedToAtom : list of set of int
            bondedToAtom[i] is the set of atoms bonded to atom index i

        Returns
        -------
        template : _TemplateData
            The matching forcefield residue template, or None if no matches are found.
        matches : list
            a list specifying which atom of the template each atom of the residue
            corresponds to, or None if it does not match the template

        """
        template = None
        matches = None
        for matcher in self._templateMatchers:
            template = matcher(self, res, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
            if template is not None:
                match = compiled.matchResidueToTemplate(res, template, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                if match is None:
                    raise ValueError('A custom template matcher returned a template for residue %d (%s), but it does not match the residue.' % (res.index, res.name))
                return [template, match]
        if templateSignatures is None:
            templateSignatures = self._templateSignatures
        signature = _createResidueSignature([atom.element for atom in res.atoms()])
        if signature in templateSignatures:
            allMatches = []
            for t in templateSignatures[signature]:
                #Note24: trimmed version: no match checking from compiled for now since we don't have compiled yet
                #match = compiled.matchResidueToTemplate(res, t, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                match=None
                if match is not None:
                    allMatches.append((t, match))
            if len(allMatches) == 1:
                template = allMatches[0][0]
                matches = allMatches[0][1]
            elif len(allMatches) > 1:
                # We found multiple matches.  This is OK if and only if they assign identical types and parameters to all atoms.
                t1, m1 = allMatches[0]
                for t2, m2 in allMatches[1:]:
                    if not t1.areParametersIdentical(t2, m1, m2):
                        raise Exception('Multiple non-identical matching templates found for residue %d (%s): %s.' % (res.index+1, res.name, ', '.join(match[0].name for match in allMatches)))
                template = allMatches[0][0]
                matches = allMatches[0][1]
        return [template, matches]

    def _buildBondedToAtomList(self, topology):
        """Build a list of which atom indices are bonded to each atom.

        Parameters
        ----------
        topology : Topology
            The Topology whose bonds are to be indexed.

        Returns
        -------
        bondedToAtom : list of list of int
            bondedToAtom[index] is the list of atom indices bonded to atom `index`

        """
        bondedToAtom = [set() for _ in topology.atoms()]
        for (atom1, atom2) in topology.bonds():
            bondedToAtom[atom1.index].add(atom2.index)
            bondedToAtom[atom2.index].add(atom1.index)
        bondedToAtom = [sorted(b) for b in bondedToAtom]
        return bondedToAtom

    def getUnmatchedResidues(self, topology):
        """Return a list of Residue objects from specified topology for which no forcefield templates are available.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.

        Returns
        -------
        unmatched_residues : list of Residue
            List of Residue objects from `topology` for which no forcefield residue templates are available.
            Note that multiple instances of the same residue appearing at different points in the topology may be returned.

        This method may be of use in generating missing residue templates or diagnosing parameterization failures.
        """
        # Find the template matching each residue, compiling a list of residues for which no templates are available.
        bondedToAtom = self._buildBondedToAtomList(topology)
        unmatched_residues = list() # list of unmatched residues
        for res in topology.residues():
            # Attempt to match one of the existing templates.
            [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom)
            if matches is None:
                # No existing templates match.
                unmatched_residues.append(res)

        return unmatched_residues

    def getMatchingTemplates(self, topology, ignoreExternalBonds=False):
        """Return a list of forcefield residue templates matching residues in the specified topology.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.
        ignoreExternalBonds : bool=False
            If true, ignore external bonds when matching residues to templates.
        Returns
        -------
        templates : list of _TemplateData
            List of forcefield residue templates corresponding to residues in the topology.
            templates[index] is template corresponding to residue `index` in topology.residues()

        This method may be of use in debugging issues related to parameter assignment.
        """
        # Find the template matching each residue, compiling a list of residues for which no templates are available.
        bondedToAtom = self._buildBondedToAtomList(topology)
        templates = list() # list of templates matching the corresponding residues
        for residue in topology.residues():
            # Attempt to match one of the existing templates.
            [template, matches] = self._getResidueTemplateMatches(residue, bondedToAtom, ignoreExternalBonds=ignoreExternalBonds)
            # Raise an exception if we have found no templates that match.
            if matches is None:
                raise ValueError('No template found for chainid <%s> resid <%s> resname <%s> (residue index within topology %d).\n%s' % (residue.chain.id, residue.id, residue.name, residue.index, _findMatchErrors(self, residue)))
            else:
                templates.append(template)

        return templates

    def generateTemplatesForUnmatchedResidues(self, topology):
        """Generate forcefield residue templates for residues in specified topology for which no forcefield templates are available.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.

        Returns
        -------
        templates : list of _TemplateData
            List of forcefield residue templates corresponding to residues in `topology` for which no forcefield templates are currently available.
            Atom types will be set to `None`, but template name, atom names, elements, and connectivity will be taken from corresponding Residue objects.
        residues : list of Residue
            List of Residue objects that were used to generate the templates.
            `residues[index]` is the Residue that was used to generate the template `templates[index]`

        """
        # Get a non-unique list of unmatched residues.
        unmatched_residues = self.getUnmatchedResidues(topology)
        # Generate a unique list of unmatched residues by comparing fingerprints.
        bondedToAtom = self._buildBondedToAtomList(topology)
        unique_unmatched_residues = list() # list of unique unmatched Residue objects from topology
        templates = list() # corresponding _TemplateData templates
        signatures = set()
        for residue in unmatched_residues:
            signature = _createResidueSignature([ atom.element for atom in residue.atoms() ])
            template = _createResidueTemplate(residue)
            is_unique = True
            if signature in signatures:
                # Signature is the same as an existing residue; check connectivity.
                for check_residue in unique_unmatched_residues:
                    matches = compiled.matchResidueToTemplate(check_residue, template, bondedToAtom, False)
                    if matches is not None:
                        is_unique = False
            if is_unique:
                # Residue is unique.
                unique_unmatched_residues.append(residue)
                signatures.add(signature)
                templates.append(template)

        return [templates, unique_unmatched_residues]

    def createSystem(self, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=None, removeCMMotion=True, hydrogenMass=None, residueTemplates=dict(),
                     ignoreExternalBonds=False, switchDistance=None, flexibleConstraints=False, drudeMass=0.4*unit.amu, **args):
        

        """Construct an NexaBind System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=None
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument.  If None (the default), it uses the
            default behavior for this force field's water model.
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.  If rigidWater is used to make water molecules
            rigid, then water hydrogens are not altered.
        residueTemplates : dict=dict()
            Key: Topology Residue object
            Value: string, name of _TemplateData residue template object to use for (Key) residue.
            This allows user to specify which template to apply to particular Residues
            in the event that multiple matching templates are available (e.g Fe2+ and Fe3+
            templates in the ForceField for a monoatomic iron ion in the topology).
        ignoreExternalBonds : boolean=False
            If true, ignore external bonds when matching residues to templates.  This is
            useful when the Topology represents one piece of a larger molecule, so chains are
            not terminated properly.  This option can create ambiguities where multiple
            templates match the same residue.  If that happens, use the residueTemplates
            argument to specify which one to use.
        switchDistance : float=None
            The distance at which the potential energy switching function is turned on for
            Lennard-Jones interactions. If this is None, no switching function will be used.
        flexibleConstraints : boolean=False
            If True, parameters for constrained degrees of freedom will be added to the System
        drudeMass : mass=0.4*amu
            The mass to use for Drude particles.  Any mass added to a Drude particle is
            subtracted from its parent atom to keep their total mass the same.
        args
            Arbitrary additional keyword arguments may also be specified.
            This allows extra parameters to be specified that are specific to
            particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        print("The Forcefield XML file parser class is under construction!")
        print("The Protein PDB file reader class is under construction!")
        print("The Energy Minimizer class is under construction!")
        print("The Report Generator class is under construction!")
        print("The Simulation class is under construction!")
        print("The Back End XmlSerilizer class is under construction!")
        print("The Back End Computaional Kernels are under construction!")
        print("The Back End VerletIntegrator class is under construction!")
        print("No Report can be generated at the moment!")
        print("New updates are coming soon for the NexaBind MD simulation package! Thank you!")
        


    def _matchAllResiduesToTemplates(self, data, topology, residueTemplates, ignoreExternalBonds, ignoreExtraParticles=False, recordParameters=True):
        """Return a list of which template matches each residue in the topology, and assign atom types."""
        templateForResidue = [None]*topology.getNumResidues()
        unmatchedResidues = []
        for chain in topology.chains():
            for res in chain.residues():
                if res in residueTemplates:
                    tname = residueTemplates[res]
                    template = self._templates[tname]
                    matches = compiled.matchResidueToTemplate(res, template, data.bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                    if matches is None:
                        raise Exception('User-supplied template %s does not match the residue %d (%s)' % (tname, res.index+1, res.name))
                else:
                    # Attempt to match one of the existing templates.
                    [template, matches] = self._getResidueTemplateMatches(res, data.bondedToAtom, ignoreExternalBonds=ignoreExternalBonds, ignoreExtraParticles=ignoreExtraParticles)
                if matches is None:
                    unmatchedResidues.append(res)
                else:
                    if recordParameters:
                        data.recordMatchedAtomParameters(res, template, matches)
                    templateForResidue[res.index] = template

        # Try to apply patches to find matches for any unmatched residues.

        if len(unmatchedResidues) > 0:
            unmatchedResidues = _applyPatchesToMatchResidues(self, data, unmatchedResidues, templateForResidue, data.bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)

        # If we still haven't found a match for a residue, attempt to use residue template generators to create
        # new templates (and potentially atom types/parameters).

        for res in unmatchedResidues:
            # A template might have been generated on an earlier iteration of this loop.
            [template, matches] = self._getResidueTemplateMatches(res, data.bondedToAtom, ignoreExternalBonds=ignoreExternalBonds, ignoreExtraParticles=ignoreExtraParticles)
            if matches is None:
                # Try all generators.
                for generator in self._templateGenerators:
                    if generator(self, res):
                        # This generator has registered a new residue template that should match.
                        [template, matches] = self._getResidueTemplateMatches(res, data.bondedToAtom, ignoreExternalBonds=ignoreExternalBonds, ignoreExtraParticles=ignoreExtraParticles)
                        if matches is None:
                            # Something went wrong because the generated template does not match the residue signature.
                            raise Exception('The residue handler %s indicated it had correctly parameterized residue %s, but the generated template did not match the residue signature.' % (generator.__class__.__name__, str(res)))
                        else:
                            # We successfully generated a residue template.  Break out of the for loop.
                            break
            if matches is None:
                raise ValueError('No template found for residue %d (%s).  %s  For more information, see https://github.com/NexaBind/NexaBind/wiki/Frequently-Asked-Questions#template' % (res.index+1, res.name, _findMatchErrors(self, res)))
            else:
                if recordParameters:
                    data.recordMatchedAtomParameters(res, template, matches)
                templateForResidue[res.index] = template
        return templateForResidue

def _findMatchErrors(forcefield, res):
    """Try to guess why a residue failed to match any template and return an error message."""
    residueCounts = _countResidueAtoms([atom.element for atom in res.atoms()])
    numResidueAtoms = sum(residueCounts.values())
    numResidueHeavyAtoms = sum(residueCounts[element] for element in residueCounts if element not in (None, elem.hydrogen))

    # Loop over templates and see how closely each one might match.

    bestMatchName = None
    numBestMatchAtoms = 3*numResidueAtoms
    numBestMatchHeavyAtoms = 2*numResidueHeavyAtoms
    for templateName in forcefield._templates:
        template = forcefield._templates[templateName]
        templateCounts = _countResidueAtoms([atom.element for atom in template.atoms])

        # Does the residue have any atoms that clearly aren't in the template?

        if any(element not in templateCounts or templateCounts[element] < residueCounts[element] for element in residueCounts):
            continue

        # If there are too many missing atoms, discard this template.

        numTemplateAtoms = sum(templateCounts.values())
        numTemplateHeavyAtoms = sum(templateCounts[element] for element in templateCounts if element not in (None, elem.hydrogen))
        if numTemplateAtoms > numBestMatchAtoms:
            continue
        if numTemplateHeavyAtoms > numBestMatchHeavyAtoms:
            continue

        # If this template has the same number of missing atoms as our previous best one, look at the name
        # to decide which one to use.

        if numTemplateAtoms == numBestMatchAtoms:
            if bestMatchName == res.name or res.name not in templateName:
                continue

        # Accept this as our new best match.

        bestMatchName = templateName
        numBestMatchAtoms = numTemplateAtoms
        numBestMatchHeavyAtoms = numTemplateHeavyAtoms
        numBestMatchExtraParticles = len([atom for atom in template.atoms if atom.element is None])

    # Return an appropriate error message.

    if numBestMatchAtoms == numResidueAtoms:
        chainResidues = list(res.chain.residues())
        if len(chainResidues) > 1 and (res == chainResidues[0] or res == chainResidues[-1]):
            return 'The set of atoms matches %s, but the bonds are different.  Perhaps the chain is missing a terminal group?' % bestMatchName
        return 'The set of atoms matches %s, but the bonds are different.' % bestMatchName
    if bestMatchName is not None:
        if numBestMatchHeavyAtoms == numResidueHeavyAtoms:
            numResidueExtraParticles = len([atom for atom in res.atoms() if atom.element is None])
            if numResidueExtraParticles == 0 and numBestMatchExtraParticles == 0:
                return 'The set of atoms is similar to %s, but it is missing %d hydrogen atoms.' % (bestMatchName, numBestMatchAtoms-numResidueAtoms)
            if numBestMatchExtraParticles-numResidueExtraParticles == numBestMatchAtoms-numResidueAtoms:
                return 'The set of atoms is similar to %s, but it is missing %d extra particles.  You can add them with Modeller.addExtraParticles().' % (bestMatchName, numBestMatchAtoms-numResidueAtoms)
        return 'The set of atoms is similar to %s, but it is missing %d atoms.' % (bestMatchName, numBestMatchAtoms-numResidueAtoms)
    return 'This might mean your input topology is missing some atoms or bonds, or possibly that you are using the wrong force field.'

def _applyPatchesToMatchResidues(forcefield, data, residues, templateForResidue, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles):
    """Try to apply patches to find matches for residues."""
    # Start by creating all templates than can be created by applying a combination of one-residue patches
    # to a single template.  The number of these is usually not too large, and they often cover a large fraction
    # of residues.

    patchedTemplateSignatures = {}
    patchedTemplates = {}
    for name, template in forcefield._templates.items():
        if name in forcefield._templatePatches:
            patches = [forcefield._patches[patchName] for patchName, patchResidueIndex in forcefield._templatePatches[name] if forcefield._patches[patchName].numResidues == 1]
            if len(patches) > 0:
                newTemplates = []
                patchedTemplates[name] = newTemplates
                _generatePatchedSingleResidueTemplates(template, patches, 0, newTemplates, set())
                for patchedTemplate in newTemplates:
                    signature = _createResidueSignature([atom.element for atom in patchedTemplate.atoms])
                    if signature in patchedTemplateSignatures:
                        patchedTemplateSignatures[signature].append(patchedTemplate)
                    else:
                        patchedTemplateSignatures[signature] = [patchedTemplate]

    # Now see if any of those templates matches any of the residues.

    unmatchedResidues = []
    for res in residues:
        [template, matches] = forcefield._getResidueTemplateMatches(res, bondedToAtom, patchedTemplateSignatures, ignoreExternalBonds, ignoreExtraParticles)
        if matches is None:
            unmatchedResidues.append(res)
        else:
            data.recordMatchedAtomParameters(res, template, matches)
            templateForResidue[res.index] = template
    if len(unmatchedResidues) == 0:
        return []

    # We need to consider multi-residue patches.  This can easily lead to a combinatorial explosion, so we make a simplifying
    # assumption: that no residue is affected by more than one multi-residue patch (in addition to any number of single-residue
    # patches).  Record all multi-residue patches, and the templates they can be applied to.

    patches = {}
    maxPatchSize = 0
    for patch in forcefield._patches.values():
        if patch.numResidues > 1:
            patches[patch.name] = [[] for i in range(patch.numResidues)]
            maxPatchSize = max(maxPatchSize, patch.numResidues)
    if maxPatchSize == 0:
        return unmatchedResidues # There aren't any multi-residue patches
    for templateName in forcefield._templatePatches:
        for patchName, patchResidueIndex in forcefield._templatePatches[templateName]:
            if patchName in patches:
                # The patch should accept this template, *and* all patched versions of it generated above.
                patches[patchName][patchResidueIndex].append(forcefield._templates[templateName])
                if templateName in patchedTemplates:
                    patches[patchName][patchResidueIndex] += patchedTemplates[templateName]

    # Record which unmatched residues are bonded to each other.

    bonds = set()
    topology = residues[0].chain.topology
    for atom1, atom2 in topology.bonds():
        if atom1.residue != atom2.residue:
            res1 = atom1.residue
            res2 = atom2.residue
            if res1 in unmatchedResidues and res2 in unmatchedResidues:
                bond = tuple(sorted((res1, res2), key=lambda x: x.index))
                if bond not in bonds:
                    bonds.add(bond)

    # Identify clusters of unmatched residues that are all bonded to each other.  These are the ones we'll
    # try to apply multi-residue patches to.

    clusterSize = 2
    clusters = bonds
    while clusterSize <= maxPatchSize:
        # Try to apply patches to clusters of this size.

        for patchName in patches:
            patch = forcefield._patches[patchName]
            if patch.numResidues == clusterSize:
                matchedClusters = _matchToMultiResiduePatchedTemplates(data, clusters, patch, patches[patchName], bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                for cluster in matchedClusters:
                    for residue in cluster:
                        unmatchedResidues.remove(residue)
                bonds = set(bond for bond in bonds if bond[0] in unmatchedResidues and bond[1] in unmatchedResidues)

        # Now extend the clusters to find ones of the next size up.

        largerClusters = set()
        for cluster in clusters:
            for bond in bonds:
                if bond[0] in cluster and bond[1] not in cluster:
                    newCluster = tuple(sorted(cluster+(bond[1],), key=lambda x: x.index))
                    largerClusters.add(newCluster)
                elif bond[1] in cluster and bond[0] not in cluster:
                    newCluster = tuple(sorted(cluster+(bond[0],), key=lambda x: x.index))
                    largerClusters.add(newCluster)
        if len(largerClusters) == 0:
            # There are no clusters of this size or larger
            break
        clusters = largerClusters
        clusterSize += 1

    return unmatchedResidues

def _countResidueAtoms(elements):
    """Count the number of atoms of each element in a residue."""
    counts = {}
    for element in elements:
        if element in counts:
            counts[element] += 1
        else:
            counts[element] = 1
    return counts

def _createResidueSignature(elements):
    """Create a signature for a residue based on the elements of the atoms it contains."""
    counts = _countResidueAtoms(elements)
    sig = []
    for c in counts:
        if c is not None:
            sig.append((c, counts[c]))
    sig.sort(key=lambda x: -x[0].mass)

    # Convert it to a string.

    s = ''
    for element, count in sig:
        s += element.symbol+str(count)
    return s

# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define two methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System.  The static method should be added to the parsers map.


## @private
class NonbondedGenerator(object):
    """A NonbondedGenerator constructs a NonbondedForce."""

    SCALETOL = 1e-5

    def __init__(self, forcefield, coulomb14scale, lj14scale, useDispersionCorrection):
        self.ff = forcefield
        self.coulomb14scale = coulomb14scale
        self.lj14scale = lj14scale
        self.useDispersionCorrection = useDispersionCorrection
        self.params = ForceField._AtomTypeParameters(forcefield, 'NonbondedForce', 'Atom', ('charge', 'sigma', 'epsilon'))

    def registerAtom(self, parameters):
        self.params.registerAtom(parameters)

    @staticmethod
    #Note24: NonbondedGenerator.parseElement is a static method within the NonbondedGenerator class.
    #Note24: Inside parseElement, it performs various operations to parse and process the content of these XML elements,
    # registering nonbonded parameters, defining force-related information, and setting up force generation 
    # based on the provided data.
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, NonbondedGenerator)]
        if 'useDispersionCorrection' in element.attrib:
            useDispersionCorrection = bool(eval(element.attrib.get('useDispersionCorrection')))
        else:
            useDispersionCorrection = None
        if len(existing) == 0:
            generator = NonbondedGenerator(
                ff,
                float(element.attrib['coulomb14scale']),
                float(element.attrib['lj14scale']),
                useDispersionCorrection
            )
            ff.registerGenerator(generator)
        else:
            # Multiple <NonbondedForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
            if (abs(generator.coulomb14scale - float(element.attrib['coulomb14scale'])) > NonbondedGenerator.SCALETOL
                or abs(generator.lj14scale - float(element.attrib['lj14scale'])) > NonbondedGenerator.SCALETOL
            ):
                raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales')
            if (
                    generator.useDispersionCorrection is not None
                    and useDispersionCorrection is not None
                    and generator.useDispersionCorrection != useDispersionCorrection
            ):
                raise ValueError('Found multiple NonbondedForce tags with different useDispersionCorrection settings.')
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        print("The Forcefield.NonbondedGenerator.createForce() function is under construction")


    #def postprocessSystem(self, sys, data, args):
        # Create the exceptions.
        print("The Forcefield.NonbondedGenerator.postprocessSystem() function is under construction")

#Note24: This line associates the NonbondedGenerator.parseElement method with the "NonbondedForce" tag within a dictionary or a similar data structure named parsers.
parsers["NonbondedForce"] = NonbondedGenerator.parseElement


## @private
class LennardJonesGenerator(object):
    """A NBFix generator to construct the L-J force with NBFIX implemented as a lookup table"""

    def __init__(self, forcefield, lj14scale, useDispersionCorrection):
        self.ff = forcefield
        self.nbfixTypes = {}
        self.lj14scale = lj14scale
        self.useDispersionCorrection = useDispersionCorrection
        self.ljTypes = ForceField._AtomTypeParameters(forcefield, 'LennardJonesForce', 'Atom', ('sigma', 'epsilon'))

    def registerNBFIX(self, parameters):
        types = self.ff._findAtomTypes(parameters, 2)
        if None not in types:
            for type1 in types[0]:
                for type2 in types[1]:
                    epsilon = _convertParameterToNumber(parameters['epsilon'])
                    sigma = _convertParameterToNumber(parameters['sigma'])
                    self.nbfixTypes[(type1, type2)] = [sigma, epsilon]
                    self.nbfixTypes[(type2, type1)] = [sigma, epsilon]

    def registerLennardJones(self, parameters):
        self.ljTypes.registerAtom(parameters)

    @staticmethod
    def parseElement(element, ff):
        #Note24: It looks for existing instances of LennardJonesGenerator within the force field's forces.
        existing = [f for f in ff._forces if isinstance(f, LennardJonesGenerator)]
        if 'useDispersionCorrection' in element.attrib:
            useDispersionCorrection = bool(eval(element.attrib.get('useDispersionCorrection')))
        else:
            useDispersionCorrection = None
        if len(existing) == 0: #Note24: If no existing instance is found, it creates a new LennardJonesGenerator object.

            generator = LennardJonesGenerator(
                ff,
                float(element.attrib['lj14scale']),
                useDispersionCorrection=useDispersionCorrection
            )
            ff.registerGenerator(generator)
        else:
            #Note24: If multiple instances are found, it checks if they have different scales for 1-4 interactions. If they do, it raises an error.
            # Multiple <LennardJonesForce> tags were found, probably in different files
            generator = existing[0]
            if abs(generator.lj14scale - float(element.attrib['lj14scale'])) > NonbondedGenerator.SCALETOL:
                raise ValueError('Found multiple LennardJonesForce tags with different 1-4 scales')
            if (
                    generator.useDispersionCorrection is not None
                    and useDispersionCorrection is not None
                    and generator.useDispersionCorrection != useDispersionCorrection
            ):
                raise ValueError('Found multiple LennardJonesForce tags with different useDispersionCorrection settings.')
        #Note24: It then iterates through XML elements such as <Atom> and <NBFixPair> to register parameters using registerLennardJones and registerNBFIX methods respectively.

        for LJ in element.findall('Atom'):
            generator.registerLennardJones(LJ.attrib)
        for Nbfix in element.findall('NBFixPair'):
            generator.registerNBFIX(Nbfix.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        # First derive the lookup tables.  We need to include entries for every type
        # that a) appears in the system and b) has unique parameters.

    #def postprocessSystem(self, sys, data, args):
        # Create the exceptions.
        print("The Forcefield.LennardJonesGenerator.postprocessSystem() function is under construction")

parsers["LennardJonesForce"] = LennardJonesGenerator.parseElement

