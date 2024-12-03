from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )
import numpy as np
from nomad.config import config
from nomad.datamodel.data import Schema
from nomad.metainfo import Quantity, SchemaPackage, Section, MSection, SubSection
from runschema.run import Run
from runschema.calculation import Calculation

configuration = config.get_plugin_entry_point(
    'apbs_parser.schema_packages:schema_package_entry_point'
)

m_package = SchemaPackage()

class APBSMol(MSection):
    m_def = Section(validate=False)

    name = Quantity(
        type=str,
        description='Pqr file of Protein')

    charge = Quantity(
        type=np.float64, shape=['*'],
        description='Atom charges')

    radius = Quantity(
        type=np.float64, shape=['*'],
        description='Atom radii')

    index = Quantity(
        type=int, shape=['*'],
        description='Atom indexes')

    atom_names = Quantity(
        type=str, shape=['*'],
        description='Atom names')

    resname = Quantity(
        type=str, shape=['*'],
        description='Residue names')

    chain = Quantity(
        type=str, shape=['*'],
        description='Chain IDs')

    resid = Quantity(
        type=int, shape=['*'],
        description='Residue indexes')

    pos = Quantity(
        type=np.float64, shape=['*', 3],
        description='Atom positions')

##################################################################

class PB_file_results(MSection):
    m_def = Section(validate=False)

    salt_concentration = Quantity(
        type=np.float64,
        description='Salt concentration mol/l')

    COMs = Quantity(
        type=np.float64, shape=['*'],
        description='COMs')

    Binding_E = Quantity(
        type=np.float64, shape=['*'],
        description='PB results')


##################################################################

class APBSIon(MSection):
    m_def = Section(validate=False)

    charge = Quantity(
        type=float,
        description='Ion charge')

    conc = Quantity(
        type=np.float64, shape=['*'],
        description='Ion concentration (mol/l)')

    radius = Quantity(
        type=float,
        description='Ion radius')


class APBSElecCalc(MSection):
    m_def = Section(validate=False)

    total_energy = SubSection(
        sub_section=simulation.calculation.EnergyEntry.m_def,
        description='''
        Total electrostatic energy
        ''')

    calc_id = Quantity(
        type=int,
        description='ID')

    grid = Quantity(
        type=float,
        shape=[3],
        description='Grid spacing')

    glen = Quantity(
        type=float,
        shape=[3],
        description='Grid side length')

class APBSElec(MSection):
    m_def = Section(validate=False)

    name = Quantity(
        type=str,
        description='''
        Name
        ''')

    calc_type = Quantity(
        type=str,
        description='''
        Calculation Type
        ''')

    pbe_type = Quantity(
        type=str,
        description='''
        Poisson-Boltzmann type
        ''')

    dime = Quantity(
        type=int,
        shape=[3],
        description='Number of grid points')

    cglen = Quantity(
        type=float,
        shape=[3],
        description='Coarse grid dimensions')

    fglen = Quantity(
        type=float,
        shape=[3],
        description='Fine grid dimensions')


    pdie = Quantity(
        type=float,
        description='Solute dielectric constant')

    sdie = Quantity(
        type=float,
        description='Solvent dielectric constant')

    srfm = Quantity(
        type=str,
        description='Model for dielectric and ion-accessibility coefficients')

    srad = Quantity(
        type=float,
        description='Radius of the solvent molecules')

    bcfl = Quantity(
        type=str,
        description='Boundary conditions for the Poisson-Boltzmann equation')

    temp = Quantity(
        type=float,
        description='Temperature')

    calcs = SubSection(sub_section=APBSElecCalc.m_def, repeats=True)

    ions = SubSection(sub_section=APBSIon.m_def, repeats=True)

##############################################################################################################

class APBSCalculation(Calculation):
    m_def = Section(validate=False, extends_base_section=False)

    elec = SubSection(sub_section=APBSElec.m_def, repeats=True)

    mols = SubSection(sub_section=APBSMol.m_def, repeats=True)

    Binding_energy = SubSection(sub_section=PB_file_results.m_def, repeats=True)


m_package.__init_metainfo__()
