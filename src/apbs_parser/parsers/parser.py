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

import datetime
import numpy as np
from pathlib import Path
import glob
import os
import re

from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg as units
from runschema.run import Run, Program
from runschema.calculation import Calculation, EnergyEntry
from apbs_parser.schema_packages.schema_package import APBSCalculation, APBSElec, APBSElecCalc, APBSIon, APBSMol, PB_file_results

configuration = config.get_plugin_entry_point(
    'apbs_parser.parsers:parser_entry_point'
)
MOL = 6.022140857e+23

C1=np.linspace(0.01 , 0.2, 8)
C2=np.linspace(0.22 , 1, 12)
C=np.concatenate((C1, C2), axis=0)
C = C.astype(float)


def str_to_sites(string):
    sym, pos = string.split('(')
    pos = np.array(pos.split(')')[0].split(',')[:3], dtype=float)
    return sym, pos


calculation_parser = UnstructuredTextFileParser(quantities=[
    Quantity('sites', r'([A-Z]\([\d\.\, \-]+\))', str_operation=str_to_sites, repeats=True),
    Quantity('energy', r'energy: (\d\.\d+)'),
    Quantity('magic_source', r'done with magic source\s*\*{3}\s*\*{3}\s*[^\d]*(\d+)', repeats=False)])

mainfile_parser = UnstructuredTextFileParser(quantities=[
    Quantity('date', r'(\d\d\d\d\/\d\d\/\d\d)', repeats=False),
    Quantity('program_version', r'super\_code\s*v(\d+)\s*', repeats=False),
    Quantity(
        'calculation', r'\s*system \d+([\s\S]+?energy: [\d\.]+)([\s\S]+\*\*\*)*',
        sub_parser=calculation_parser,
        repeats=True)
])

def parse_apbs(parent, filepath):
    pqr_files = []
    with open(filepath) as file:
        mode = 'beginning'
        for i, line in enumerate(file):
            parts = line.split()
            if len(parts) == 0:
                continue

            parts[0] = parts[0].lower()

            if mode == 'beginning':
                if parts[0] == 'read':
                    mode = 'read'
                elif parts[0] == 'elec':

                    elec = APBSElec()
                    parent.elec.append(elec)
                    if len(parts) == 3 and parts[1] == 'name':
                        elec.name = parts[2]

                    mode = 'elec'
            elif mode == 'read':
                if parts[0] == 'mol':
                    if len(parts) < 3 or parts[1] != 'pqr':
                        continue
                    pqr_files.append(parts[2])
                elif parts[0] == 'end':
                    mode = 'beginning'
            elif mode == 'elec':
                if parts[0] == 'mol':
                    elec.mol = parts[1]
                elif parts[0] in ('mg-auto', 'mg-para', 'mg-manual'):
                    elec.calc_type = parts[0]
                elif parts[0] in ('lpbe', 'npbe', 'lrpbe', 'nrpbe'):
                    elec.pbe_type = parts[0]
                elif parts[0] == 'dime':
                    elec.dime = [int(s) for s in parts[1:4]]
                elif parts[0] == 'cglen':
                    elec.cglen = [float(s) for s in parts[1:4]]
                elif parts[0] == 'fglen':
                    elec.fglen = [float(s) for s in parts[1:4]]
                elif parts[0] == 'pdie':
                    elec.pdie = float(parts[1])
                elif parts[0] == 'sdie':
                    elec.sdie = float(parts[1])
                elif parts[0] == 'srfm':
                    elec.srfm = parts[1]
                elif parts[0] == 'srad':
                    elec.srad = float(parts[1])
                elif parts[0] == 'bcfl':
                    elec.bcfl = parts[1]
                elif parts[0] == 'temp':
                    elec.temp = float(parts[1])
                elif parts[0] == 'ion':
                    ion = APBSIon()
                    elec.ions.append(ion)

                    partmode = None
                    for s in parts[1:]:
                        s = s.lower()
                        if partmode == 'charge':
                            ion.charge = float(s)
                            partmode = None
                        elif partmode == 'conc':
                            potentials_directory = os.path.join(filepath.parent, "Potentials")
                            if os.path.exists(potentials_directory):
                                ion.conc=C
                            else:
                                ion.conc = str(s)
                            partmode = None
                        elif partmode == 'radius':
                            ion.radius = float(s)
                            partmode = None
                        elif s in ('charge', 'conc', 'radius'):
                            partmode = s
                elif parts[0] == 'calc':
                    mode = 'calc'
                    elec_calc = APBSElecCalc()
                    elec.calcs.append(elec_calc)
                elif parts[0] == 'end':
                    mode = 'beginning'
            elif mode == 'calc':
                if parts[0] == 'id':
                    elec_calc.calc_id = int(parts[1])
                if parts[0] == 'grid':
                    elec_calc.grid = [float(s) for s in parts[1:4]]
                if parts[0] == 'glen':
                    elec_calc.glen = [float(s) for s in parts[1:4]]
                elif parts[0] == 'totenergy':
                    elec_calc.total_energy = EnergyEntry(value=float(parts[1]) *  units.kJ / MOL)
                elif parts[0] == 'end':
                    mode = 'elec'

    return pqr_files

########################################################################################################

def extract_salt_concentration(filename):
    # Use regular expression to find the salt concentration value in the filename
    match = re.search(r"\d+\.\d+", filename)
    if match:
        return float(match.group())
    else:
        return None

#########################################################################################################

class APBSLogParser:
    def __init__(self):
        pass

    def parse(self, run, filepath):
        calculation = APBSCalculation()
        run.calculation.append(calculation)
        parse_apbs(calculation, filepath)

###########################################################################################################################


class PBFileParser:
    def parse(self, calculation, directory):
        potentials_directory = os.path.join(directory, "Potentials")
        file_list = glob.glob(os.path.join(potentials_directory, "AB_*.dat"))

        for file_path in file_list:
            self.parse_file(calculation, file_path)

    def parse_file(self, calc, file_path):
        free_energy_diff = PB_file_results()
        calculation.Binding_energy.append(free_energy_diff) # Replace `AB_file_results` with the actual metainfo section for AB_.dat data

        salt_conc = extract_salt_concentration(file_path)
        # Parse and populate the AB_.dat file data into the AB_file_results section
        with open(file_path) as file:
            data = np.loadtxt(file)
            free_energy_diff.salt_concentration=salt_conc
            free_energy_diff.COMs = data[:, 0]*10**9
            free_energy_diff.Binding_E = data[:, 1]



###############################################################################################################################


class APBSInputParser():
    def __init__(self):
        self.pqr_parser = APBSPQRParser()

    def parse(self, calculation, filepath):
        pqr_files = parse_apbs(calculation, filepath)

        mainpath = Path(filepath).parent

        for filename in pqr_files:
            p = mainpath / filename
            if p.is_file():
                self.pqr_parser.parse(calculation, mainpath / filename)

class APBSPQRParser():
    def __init__(self):
        pass

    def parse(self, system, filename):
        mols = APBSMol()
        system.mols.append(mols)

        charges = []
        radii=[]
        indexes=[]
        names=[]
        resnames=[]
        chains=[]
        resids=[]
        positions=[]
        with open(filename) as file:
            for line in file:
                parts = line.split()
                if parts[0] == 'ATOM' or parts[0] == 'HETATOM':
                    #atom = mol.m_create(APBSAtom)
                    charges.append(float(parts[-2]))
                    radii.append(float(parts[-1]))
                    indexes.append(int(parts[1]))
                    names.append(parts[2])
                    resnames.append(parts[3])
                    if len(parts) >= 11:
                        chains.append(parts[4])
                        resids.append(int(parts[5]))
                    else:
                        resids.append(int(parts[4]))
                    positions.append([float(s) for s in parts[-5:-2]])

            mols.name=os.path.basename(filename)
            mols.charge = charges
            mols.radius = radii
            mols.index = indexes
            mols.atom_names = names
            mols.resname = resnames
            mols.chain = chains
            mols.resid = resids
            mols.pos = positions


class NewParser(MatchingParser):
    def __init__(self):
        self.input_parser = APBSInputParser()
        self.pqr_parser = APBSPQRParser()
        self.file_parser = PBFileParser()
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None:
        logger.info('NewParser.parse', parameter=configuration.parameter)

        archive.workflow2 = Workflow(name='test')
        mainfile = Path(mainfile)
        run = Run()
        archive.run.append(run)
        run.program = Program(name="NOMAD parser for APBS")

        calculation = APBSCalculation()
        run.calculation.append(calculation)

        self.input_parser.parse(calculation, mainfile)
        self.file_parser.parse(calculation, mainfile.parent)