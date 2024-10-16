# ---------------------------------------------------------------------------------------------------------------
#
# This is a parser for Nanomatch GmbH's kMC software "Lightforge"
# Important: For this parser to work properly, make sure that the very last row in the Lightforge
# settings file starts with any letter, in other words it must not start with space or with a dash.
#


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

from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.parser import MatchingParser

#from nomad.datamodel.metainfo.simulation.run import Run                        
#from nomad.datamodel.metainfo.simulation.calculation import Calculation

from runschema.run import Run, Program                                         
from runschema.calculation import Calculation
#from nomad_simulations.schema_packages.general import Program, Simulation


import yaml
import os
import re
import datetime
import numpy as np
import filecmp
from pathlib import Path


from nomad.parsing.parser import MatchingParser
from runschema.run import Run, Program
from runschema.calculation import Calculation

#from nomad.datamodel import EntryArchive
#from nomad.units import ureg as units
#from nomad.datamodel.results import Results, Properties, Structure
#from nomad.parsing.file_parser import UnstructuredTextFileParser, Quantity
#from nomad.datamodel.optimade import Species
#from . import metainfo  # pylint: disable=unused-import
from nomad.datamodel.metainfo.workflow import Workflow
from lightforge_v2.schema_packages.schema_package import (
                            LightforgeCalculation, IV, IQE2, Current_density, Current_characteristics, LF_experiments, LF_material,
                            LF_input, Settings, Settings_pair_input, Settings_materials, Settings_layers, 
                            Layer_molecule_species, Settings_electrodes, Settings_hole_transfer_integrals, 
                            Settings_electron_transfer_integrals, Settings_dexter_transfer_integrals, 
                            Settings_qp_output_files, Run_lf_slr, Files_for_kmc, Js_homo_mol_pairs, Js_lumo_mol_pairs, Js_dexter_mol_pairs,
                            Sigma_mol_pairs, LF_molecule_pdb_file, LF_vacuum_lambda,
                            LF_add_info, Material_data, Lightforge_data, Runtime_data, LF_experiment_inventory, 
                            LF_particle_positions, Mobility, Particle_densities, 
                            Charge_density_average, Exciton_decay_density_average, 
                            Photon_creation_density_average,
                            Quenching_density_average, Exciton_molpairs, Emitter_emitter_transport_count,
                            Host_emitter_transport_count, Host_host_transport_count, Runtime_analysis,
                            Event_counts_by_type, Device_data, Electrodes, Energy_levels,
                            Exciton_separation, Foerster, Site_energies, Mol_types, Coordinates,
                            Dexter_and_foerster, Foerster_expansion_errors) 

configuration = config.get_plugin_entry_point(
    'lightforge_v2.parsers:parser_entry_point'
)


def DetailedParser(filepath, archive):
    run = Run()
    archive.run.append(run)
    
    calculation = LightforgeCalculation()
    run.calculation.append(calculation)
    
    lf_experiments = LF_experiments()
    calculation.lf_experiments = lf_experiments
    
    lf_material = LF_material()
    calculation.lf_material = lf_material
    
    current_characteristics = Current_characteristics()
    lf_experiments.current_characteristics = current_characteristics
    
    particle_densities = Particle_densities()
    lf_experiments.particle_densities = particle_densities
    
    iqe2 = IQE2()
    current_characteristics.IQE2 = iqe2
    
    iv = IV()
    current_characteristics.IV = iv
    
    lightforge_data = Lightforge_data()
    calculation.lightforge_data = lightforge_data
    
    lf_input = LF_input()
    calculation.lf_input = lf_input

    exciton_molpairs_hasrun = False
    runtime_analysis_hasrun = False
#    device_data_hasrun= False
    foerster_hasrun = False
    material_data_hasrun = False
    files_for_kmc_hasrun = False
    runtime_data_hasrun = False
    _coordinates = []
    coordinates_counter = 0
    coordinates_rows = [0]
    
#    current_density = Current_density()


    for root, dirs, files in sorted(os.walk(filepath.parent)):        
        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]
        
        files = sorted(files, key = natsort)
        i = 0
        while i < len(files):
            if '.png' in files[i] or '.npz' in files[i] or '.zip' in files[i] or '.out' in files[i]:
                files.remove(files[i])
            else:
                i += 1    
        for file in files:    
            with open(root +'/'+ file, 'r') as f:
                
                if 'current_density' in file and 'all_data_points' not in root:
                    current_density = Current_density()
                    current_characteristics.current_density.append(current_density)
                    value = []
                    for i, line in enumerate(f):
                        line = float(line)
                        value.append(line)

                    current_density.value = np.array(value)
                if 'IQE2_all_currents' in file and 'all_data_points' not in root:
                    for i, line in enumerate(f):
                        rows = i+1
                    
                    a = np.zeros((rows,2))
                    iqe2.iqe2_all_currents = a

                    f.seek(0)
                    for i, line in enumerate(f):
                        
                        parts = line.split()
                        
                        iqe2.iqe2_all_currents[i][0] = parts[0]
                        iqe2.iqe2_all_currents[i][1] = parts[1]
                if 'IQE2_all_fields' in file and 'all_data_points' not in root:
                    for i, line in enumerate(f):
                        rows = i + 1
                    b = np.zeros((rows,2))
                    iqe2.iqe2_all_fields = b
                    f.seek(0)
                    for i, line in enumerate(f):
                        parts = line.split()
                        iqe2.iqe2_all_fields[i][0] = parts[0]
                        iqe2.iqe2_all_fields[i][1] = parts[1] 
                if re.search(r'^IV_all_fields.dat$', file) and 'all_data_points' not in root:
                    for i, line in enumerate(f):
                        rows = i + 1
                    c = np.zeros((rows,3))
                    iv.iv_all_fields = c
                    f.seek(0)
                    for i, line in enumerate(f):
                        parts = line.split()
                        
                        iv.iv_all_fields[i][0] = parts[0]
                        iv.iv_all_fields[i][1] = parts[1]                 
                        iv.iv_all_fields[i][2] = parts[2]                          
                if re.search(r'mobilities_\d+', file) and 'all_data_points' not in root:
                    mobility = Mobility()
                    current_characteristics.mobility.append(mobility)
                    value = []
                    for i, line in enumerate(f):
                        line=float(line)
                        value.append(line)
                    mobility.value = np.array(value)
                
                if re.search(r'mobilities_all_fields', file) and 'all_data_points' not in root:
                    mobility = Mobility()
                    current_characteristics.mobility.append(mobility)
                    for i, line in enumerate(f):
                        rows  = i +1 
                    d = np.zeros((rows, 3))
                    mobility.mobilities_all_fields = d
                    f.seek(0)
                    for i, line in enumerate(f):
                        parts = line.split()
                        mobility.mobilities_all_fields[i][0] = parts[0]
                        mobility.mobilities_all_fields[i][1] = parts[1]
                        mobility.mobilities_all_fields[i][2] = parts[2]
                
                
                if re.search(r'charge_density_average_\d+.dat', file) and 'all_data_points' not in root:
                    charge_density_average = Charge_density_average()
                    particle_densities.charge_density_average.append(charge_density_average)
                    device_length = []
                    electrons = []
                    holes = []
                    for i, line in enumerate(f):
                        columns = len(line.split())
                        break
                    f.seek(0)
                    value = np.zeros((3, columns))
                    for i, line in enumerate(f):
                        if i==0:
                            parts = line.split()
                            device_length = parts
                            
                        if i==1:
                            parts = line.split()
                            electrons = parts
                            
                        if i==2:
                            parts = line.split()
                            holes = parts
                    value[0] = device_length
                    value[1] = electrons
                    value[2] = holes  
                    charge_density_average.value = value
                '''
                if re.search(r'exciton_decay_density_average_\d+', file)  and 'all_data_points'  not in root:
                    sec_exciton_decay_density_average  = sec_particle_densities.m_create(Exciton_decay_density_average)
                    file_exciton_decay_density_average =  yaml.safe_load(f)
                                        
                    if 'CS' in file_exciton_decay_density_average['total']['annihilation']:
                        _cs = file_exciton_decay_density_average['total']['annihilation']['CS']
                        sec_exciton_decay_density_average.cs = _cs    
                    if 'EET' in file_exciton_decay_density_average['total']['annihilation']:
                        _eet = file_exciton_decay_density_average['total']['annihilation']['EET']
                        sec_exciton_decay_density_average.eet = _eet
                    if 'EPT' in file_exciton_decay_density_average['total']['annihilation']:
                        _ept = file_exciton_decay_density_average['total']['annihilation']['EPT']
                        sec_exciton_decay_density_average.ept = _ept                    
                    if 'PTQ' in file_exciton_decay_density_average['total']['annihilation']:
                        _ptq = file_exciton_decay_density_average['total']['annihilation']['PTQ']
                        sec_exciton_decay_density_average.ptq = _ptq
                    if 'SPQ' in file_exciton_decay_density_average['total']['annihilation']:
                        _spq = file_exciton_decay_density_average['total']['annihilation']['SPQ']
                        sec_exciton_decay_density_average.spq= _spq
                    if 'SSA' in file_exciton_decay_density_average['total']['annihilation']:
                        _ssa = file_exciton_decay_density_average['total']['annihilation']['SSA']
                        sec_exciton_decay_density_average.ssa = _ssa
                    if 'STA' in file_exciton_decay_density_average['total']['annihilation']:
                        _sta = file_exciton_decay_density_average['total']['annihilation']['STA']
                        sec_exciton_decay_density_average.sta = _sta
                    if 'TPQ' in file_exciton_decay_density_average['total']['annihilation']:
                        _tpq = file_exciton_decay_density_average['total']['annihilation']['TPQ']
                        sec_exciton_decay_density_average.tpq = _tpq
                    if 'TSA' in file_exciton_decay_density_average['total']['annihilation']:
                        _tsa = file_exciton_decay_density_average['total']['annihilation']['TSA']
                        sec_exciton_decay_density_average.tsa = _tsa 
                    if 'TTA' in file_exciton_decay_density_average['total']['annihilation']:
                        _tta = file_exciton_decay_density_average['total']['annihilation']['TTA']
                        sec_exciton_decay_density_average.tta = _tta
                    if 'TTF' in file_exciton_decay_density_average['total']['annihilation']:
                        _ttf = file_exciton_decay_density_average['total']['annihilation']['TTF']
                        sec_exciton_decay_density_average.ttf = _ttf
                    if 'radiative' in file_exciton_decay_density_average['total']['annihilation']:
                        _radiative = file_exciton_decay_density_average['total']['annihilation']['radiative']
                        sec_exciton_decay_density_average.radiative = _radiative
                    if 'thermal' in file_exciton_decay_density_average['total']['annihilation']:
                        _thermal = file_exciton_decay_density_average['total']['annihilation']['thermal']
                        sec_exciton_decay_density_average.thermal = _thermal
                    if 'photon' in file_exciton_decay_density_average['total']['creation']:
                        _photon = file_exciton_decay_density_average['total']['creation']['photon']
                        sec_exciton_decay_density_average.photon = _photon
                    if 'recombination' in file_exciton_decay_density_average['total']['creation']:
                        _recombination = file_exciton_decay_density_average['total']['creation']['recombination']
                        sec_exciton_decay_density_average.recombination = _recombination
                    if 'x_axis' in file_exciton_decay_density_average:
                        _exciton_decay_density_average_x_axis = file_exciton_decay_density_average['x_axis']
                        sec_exciton_decay_density_average.exciton_decay_density_average_x_axis = _exciton_decay_density_average_x_axis    
                
                if re.search(r'photon_creation_density_average_\d+', file) and 'all_data_points' not in root:
                    sec_photon_creation_density_average = sec_particle_densities.m_create(Photon_creation_density_average)
                    device_length = []
                    photons = []
                    excitons = []
                    for i, line in enumerate(f):
                        columns = len(line.split())
                        break
                    f.seek(0)
                    value = np.zeros((3, columns))
                    for i, line in enumerate(f):
                        if i == 0:
                            parts = line.split()
                            device_length = parts
                        if i == 1:
                            parts = line.split()
                            photons = parts
                        if i == 2:
                            parts = line.split()
                            excitons = parts
                    value[0] = device_length
                    value[1] = photons
                    value[2] = excitons
                    sec_photon_creation_density_average.value = value
                if re.search(r'quenching_density_average_\d+', file) and 'all_data_points' not in root:
                    sec_quenching_density_average = sec_particle_densities.m_create(Quenching_density_average)
                    device_length = []
                    excitons_quenched = []
                    electrons = []
                    holes = []
                    for i, line in enumerate(f):
                        columns = len(line.split())
                        break
                    f.seek(0)
                    value = np.zeros((4, columns))
                    for i, line in enumerate(f):
                        if i == 0:
                            parts = line.split()
                            device_length = parts
                        if i == 1:
                            parts = line.split()
                            excitons_quenched = parts
                        if i == 2:
                            parts = line.split()
                            electrons = parts
                        if i == 3:
                            parts = line.split()
                            holes = parts
                    value[0] = device_length     
                    value[1] = excitons_quenched
                    value[2] = electrons
                    value[3] = holes 
                    sec_quenching_density_average.value = value
                if 'exciton_molpairs' in root:
                    if not exciton_molpairs_hasrun:
                        sec_exciton_molpairs = sec_particle_densities.m_create(Exciton_molpairs)
                        exciton_molpairs_hasrun = True   
                    if re.search(r'host_emitter_transport_count', file) and 'all_data_points' not in root:
                        sec_host_emitter_transport_count = sec_exciton_molpairs.m_create(Host_emitter_transport_count) # might be wrong
                        file_host_emitter_transport_count = yaml.safe_load(f)
                        if 'dexter_S1S1' in file_host_emitter_transport_count:
                            _dexter_s1s1 = file_host_emitter_transport_count['dexter_S1S1']
                            sec_host_emitter_transport_count.dexter_s1s1 = _dexter_s1s1
                        if 'dexter_S1T1' in file_host_emitter_transport_count:
                            _dexter_s1t1 = file_host_emitter_transport_count['dexter_S1T1']
                            sec_host_emitter_transport_count.dexter_s1t1 = _dexter_s1t1
                        if 'dexter_T1S1' in file_host_emitter_transport_count:
                            _dexter_t1s1 = file_host_emitter_transport_count['dexter_T1S1']
                            sec_host_emitter_transport_count.dexter_t1s1 = _dexter_t1s1
                        if 'dexter_T1T1' in file_host_emitter_transport_count:
                            _dexter_t1t1 = file_host_emitter_transport_count['dexter_T1T1']
                            sec_host_emitter_transport_count.dexter_t1t1 = _dexter_t1t1
                        if 'foerster_S1S1' in file_host_emitter_transport_count:
                            _foerster_s1s1 = file_host_emitter_transport_count['foerster_S1S1']
                            sec_host_emitter_transport_count.foerster_s1s1 = _foerster_s1s1
                        if 'foerster_S1T1' in file_host_emitter_transport_count:
                            _foerster_s1t1 = file_host_emitter_transport_count['foerster_S1T1']
                            sec_host_emitter_transport_count.foerster_s1t1 = _foerster_s1t1
                        if 'foerster_T1S1' in file_host_emitter_transport_count:
                            _foerster_t1s1 = file_host_emitter_transport_count['foerster_T1S1']
                            sec_host_emitter_transport_count.foerster_t1s1 = _foerster_t1s1
                        if 'foerster_T1T1' in file_host_emitter_transport_count:
                            _foerster_t1t1 = file_host_emitter_transport_count['foerster_T1T1']
                            sec_host_emitter_transport_count.foerster_t1t1 = _foerster_t1t1
                        if 'x_axis' in file_host_emitter_transport_count:
                            _x_axis = file_host_emitter_transport_count['x_axis']
                            sec_host_emitter_transport_count.x_axis = _x_axis                                         
                    if re.search(r'emitter_emitter_transport_count.yml', file) and 'all_data_points' not in root:
                        sec_emitter_emitter_transport_count = sec_exciton_molpairs.m_create(Emitter_emitter_transport_count)
                        file_emitter_emitter_transport_count = yaml.safe_load(f)
                        if 'dexter_S1S1' in file_emitter_emitter_transport_count:
                            _dexter_s1s1 = file_emitter_emitter_transport_count['dexter_S1S1']
                            sec_emitter_emitter_transport_count.dexter_s1s1 = _dexter_s1s1
                        if 'dexter_S1T1' in file_emitter_emitter_transport_count:
                            _dexter_s1t1 = file_emitter_emitter_transport_count['dexter_S1T1']
                            sec_emitter_emitter_transport_count.dexter_s1t1 = _dexter_s1t1
                        if 'dexter_T1S1' in file_emitter_emitter_transport_count:
                            _dexter_t1s1 = file_emitter_emitter_transport_count['dexter_T1S1']
                            sec_emitter_emitter_transport_count.dexter_t1s1 = _dexter_t1s1
                        if 'dexter_T1T1' in file_emitter_emitter_transport_count:
                            _dexter_t1t1 = file_emitter_emitter_transport_count['dexter_T1T1']
                            sec_emitter_emitter_transport_count.dexter_t1t1 = _dexter_t1t1
                        if 'foerster_S1S1' in file_emitter_emitter_transport_count:
                            _foerster_s1s1 = file_emitter_emitter_transport_count['foerster_S1S1']
                            sec_emitter_emitter_transport_count.foerster_s1s1 = _foerster_s1s1
                        if 'foerster_S1T1' in file_emitter_emitter_transport_count:
                            _foerster_s1t1 = file_emitter_emitter_transport_count['foerster_S1T1']
                            sec_emitter_emitter_transport_count.foerster_s1t1 = _foerster_s1t1
                        if 'foerster_T1S1' in file_emitter_emitter_transport_count:
                            _foerster_t1s1 = file_emitter_emitter_transport_count['foerster_T1S1']
                            sec_emitter_emitter_transport_count.foerster_t1s1 = _foerster_t1s1
                        if 'foerster_T1T1' in file_emitter_emitter_transport_count:
                            _foerster_t1t1 = file_emitter_emitter_transport_count['foerster_T1T1']
                            sec_emitter_emitter_transport_count.foerster_t1t1 = _foerster_t1t1
                        if 'x_axis' in file_emitter_emitter_transport_count:
                            _x_axis = file_emitter_emitter_transport_count['x_axis']
                            sec_emitter_emitter_transport_count.x_axis = _x_axis
                    if re.search(r'host_host_transport_count.yml', file) and 'all_data_points' not in root:
                        sec_host_host_transport_count = sec_exciton_molpairs.m_create(Host_host_transport_count)
                        file_host_host_transport_count = yaml.safe_load(f)
                        if 'dexter_S1S1' in file_host_host_transport_count:
                            _dexter_s1s1 = file_host_host_transport_count['dexter_S1S1']
                            sec_host_host_transport_count.dexter_s1s1 = _dexter_s1s1 
                        if 'dexter_S1T1' in file_host_host_transport_count:
                            _dexter_s1t1 = file_host_host_transport_count['dexter_S1T1']
                            sec_host_host_transport_count.dexter_s1t1 = _dexter_s1t1
                        if 'dexter_T1S1' in file_host_host_transport_count:
                            _dexter_t1s1 = file_host_host_transport_count['dexter_T1S1']
                            sec_host_host_transport_count.dexter_t1s1 = _dexter_t1s1
                        if 'dexter_T1T1' in file_host_host_transport_count:
                            _dexter_t1t1 = file_host_host_transport_count['dexter_T1T1']
                            sec_host_host_transport_count.dexter_t1t1 = _dexter_t1t1
                        if 'foerster_S1S1' in file_host_host_transport_count:
                            _foerster_s1s1 = file_host_host_transport_count['foerster_S1S1']
                            sec_host_host_transport_count.foerster_s1s1 = _foerster_s1s1
                        if 'foerster_S1T1' in file_host_host_transport_count:
                            _foerster_s1t1 = file_host_host_transport_count['foerster_S1T1']
                            sec_host_host_transport_count.foerster_s1t1 = _foerster_s1t1
                        if 'foerster_T1S1' in file_host_host_transport_count:
                            _foerster_t1s1 = file_host_host_transport_count['foerster_T1S1']
                            sec_host_host_transport_count.foerster_t1s1 = _foerster_t1s1
                        if 'foerster_T1T1' in file_host_host_transport_count:
                            _foerster_t1t1 = file_host_host_transport_count['foerster_T1T1']
                            sec_host_host_transport_count.foerster_t1t1 = _foerster_t1t1
                        if 'x_axis' in file_host_host_transport_count:
                            _x_axis = file_host_host_transport_count['x_axis']
                            sec_host_host_transport_count.x_axis = _x_axis
                
                if 'runtime_analysis' in root:
                    
                    if not runtime_analysis_hasrun:    
                        sec_runtime_analysis = sec_experiments.m_create(Runtime_analysis)
                        runtime_analysis_hasrun = True
                    if re.search(r'event_counts_by_type_\d+', file) and not 'all_data_points' in root:
                        sec_event_counts_by_type = sec_runtime_analysis.m_create(Event_counts_by_type)
                        _spq = []
                        _sta = []
                        _tpq = []
                        _tta = []
                        _ttf = []
                        for i, line in enumerate(f):
                            line = line.lower()
                            parts = line.split(': ')
                            if 'dexter eeq' in line:
                                sec_event_counts_by_type.dexter_eeq = float(parts[1])
                            if 'dexter ept' in line:
                                sec_event_counts_by_type.dexter_ept = float(parts[1])
                            if 'spq' in line:
                                _spq.append(float(parts[1]))
                                sec_event_counts_by_type.spq = _spq       
                            if 'sta' in line:
                                _sta.append(float(parts[1]))
                                sec_event_counts_by_type.sta = _sta
                            if 'tpq' in line:
                                _tpq.append(float(parts[1]))
                                sec_event_counts_by_type.tpq = _tpq
                            if 'tta' in line:
                                _tta.append(float(parts[1]))
                                sec_event_counts_by_type.tta = _tta
                            if 'ttf' in line:
                                _ttf.append(float(parts[1]))
                                sec_event_counts_by_type.ttf = _ttf
                            if 'eject chg' in line:
                                sec_event_counts_by_type.eject_chg = float(parts[1])
                            if 'inject e' in line:
                                sec_event_counts_by_type.inject_e = float(parts[1])
                            if 'inject h' in line:
                                sec_event_counts_by_type.inject_h = float(parts[1])
                            if 'move chg' in line:
                                sec_event_counts_by_type.move_chg = float(parts[1])
                            if 'move exc dexter' in line:
                                sec_event_counts_by_type.move_exc_dexter = float(parts[1])
                            if 'move exc foerster' in line:
                                sec_event_counts_by_type.move_exc_foerster = float(parts[1])
                            if 'move+flip exc dexter' in line:
                                sec_event_counts_by_type.move_flip_exc_dexter = float(parts[1])
                            if 'move+flip exc foerster' in line:
                                sec_event_counts_by_type.move_flip_exc_foerster = float(parts[1])
                            if 'prtclrst' in line:
                                sec_event_counts_by_type.prtclRst = float(parts[1])
                            if 'rad decay' in line:
                                sec_event_counts_by_type.rad_decay = float(parts[1])
                            if 'recombination s1' in line:
                                sec_event_counts_by_type.recombination_s1 = float(parts[1])
                            if 'recombination t1' in line:
                                sec_event_counts_by_type.recombination_t1 = float(parts[1])
                            if 'seperate eh' in line:
                                sec_event_counts_by_type.seperate_eh = float(parts[1])
                            if 'seperate he' in line:
                                sec_event_counts_by_type.seperate_he = float(parts[1])
                            if 'setmult' in line:
                                sec_event_counts_by_type.setMult = float(parts[1])
                            if 'spin flip exc' in line:
                                sec_event_counts_by_type.spin_flip_exc = float(parts[1])
                            if 'thermal_decay' in line:
                                sec_event_counts_by_type.thermal_decay = float(parts[1])
    '''
    '''                                                                                         # ab hier ausklammern behalten 
                if 'device_data' in root:
                    if not device_data_hasrun:
                        sec_device_data = sec_material.m_create(Device_data)
                        device_data_hasrun = True
                    
                    if re.search(r'coord_\d+', file) and 'all_data_points' not in root:
                        sec_coordinates = sec_device_data.m_create(Coordinates)
                        
                        
                        for i, line in enumerate(f):
                            parts = line.split()
                            parts = [float(p) for p in parts]
                            _coordinates.append(parts)
                        
                        
                        if _coordinates[sum(coordinates_rows):] == _coordinates[sum(coordinates_rows) - coordinates_rows[-1]:sum(coordinates_rows)] and coordinates_counter >= 1:
                            sec_coordinates.text = 'This file has the same content as the previous file/previous repeating subsection.'
                        else: 
                            sec_coordinates.coordinates = _coordinates[sum(coordinates_rows):]
                        coordinates_counter += 1
                        coordinates_rows.append(i+1)
                    
                    if re.search(r'mol_types_\d+', file) and 'all_data_points' not in root:
                        sec_mol_types = sec_device_data.m_create(Mol_types)
                        _mol_types = []
                        for i, line in enumerate(f):
                            line = float(line)
                            _mol_types.append(line)
                        sec_mol_types.mol_types = _mol_types
                    
                    if re.search(r'site_energies_\d+', file) and 'all_data_points' not in root:
                        sec_site_energies = sec_device_data.m_create(Site_energies)
                        _site_energies = []
                        for i, line in enumerate(f):
                            parts = line.split()
                            parts = [float(p) for p in parts]
                            _site_energies.append(parts)
                        sec_site_energies.site_energies = _site_energies                        # bis hier
    '''                                                                                               
    '''            
                if 'Foerster' in root:
                    if not foerster_hasrun:
                        sec_foerster = sec_material.m_create(Foerster)
                        foerster_hasrun = True
                    if re.search(r'Dexter_\d+_', file) or re.search(r'[a-zA-Z]+\d[a-zA-Z]\d_', file):
                        sec_dexter_and_foerster = sec_foerster.m_create(Dexter_and_foerster)
                        _values = []
                        for i, line in enumerate(f):
                            parts = line.split()
                            parts = [float(x) for x in parts]
                            _values.append(parts)
                        sec_dexter_and_foerster.name = file
                        sec_dexter_and_foerster.values = _values
                    if re.search(r'foerster_expansion_errors', file) and 'all_data_points' not in root:
                        sec_foerster_expansion_errors = sec_foerster.m_create(Foerster_expansion_errors)
                        
                        for i, line in enumerate(f):
                            parts = line.split(': ')
                            if 'S1S1_0_0' in line:
                                sec_foerster_expansion_errors.s1s1_0_0 = float(parts[2])
                            if 'T1T1_0_0' in line:
                                sec_foerster_expansion_errors.t1t1_0_0 = float(parts[2])
                            if 'S1T1_0_0' in line:
                                sec_foerster_expansion_errors.s1t1_0_0 = float(parts[2])
                            if 'T1S1_0_0' in line:
                                sec_foerster_expansion_errors.t1s1_0_0 = float(parts[2])
                            if 'S1S1_0_1' in line:
                                sec_foerster_expansion_errors.s1s1_0_1 = float(parts[2])
                            if 'T1T1_0_1' in line:
                                sec_foerster_expansion_errors.t1t1_0_1 = float(parts[2])
                            if 'S1T1_0_1' in line:
                                sec_foerster_expansion_errors.s1t1_0_1 = float(parts[2])
                            if 'T1S1_0_1' in line:
                                sec_foerster_expansion_errors.t1s1_0_1 = float(parts[2])
                            if 'S1S1_1_0' in line:
                                sec_foerster_expansion_errors.s1s1_1_0 = float(parts[2])
                            if 'T1T1_1_0' in line:
                                sec_foerster_expansion_errors.t1t1_1_0 = float(parts[2])
                            if 'S1T1_1_0' in line:
                                sec_foerster_expansion_errors.s1t1_1_0 = float(parts[2])
                            if 'T1S1_1_0' in line:
                                sec_foerster_expansion_errors.t1s1_1_0 = float(parts[2])
                            if 'S1S1_1_1' in line:
                                sec_foerster_expansion_errors.s1s1_1_1 = float(parts[2])
                            if 'T1T1_1_1' in line:
                                sec_foerster_expansion_errors.t1t1_1_1 = float(parts[2])
                            if 'S1T1_1_1' in line:
                                sec_foerster_expansion_errors.s1t1_1_1 = float(parts[2])
                            if 'T1S1_1_1' in line:
                                sec_foerster_expansion_errors.t1s1_1_1 = float(parts[2]) 
                if 'material_data' in root:
                    if not material_data_hasrun:
                        sec_material_data = sec_lightforge_data.m_create(Material_data)
                        material_data_hasrun = True
                    if re.search(r'add_info_\d+', file):
                        sec_lf_add_info = sec_material_data.m_create(LF_add_info)
                        file_add_info = yaml.safe_load(f)
                        _length_layers = len(file_add_info['layers'])
                        
                        _lf_layer_id = []
                        _n_layer_sites = []
                        _sites_end_idx_in_device = []
                        _sites_start_idx_in_device = []
                        _add_info_thickness = []
                        _add_info_x_boundaries = []
                        for i in range(_length_layers):

                            if 'layer_id' in file_add_info['layers'][i]:
                                _lf_layer_id.append(int(file_add_info['layers'][i]['layer_id']))
                                
                            if 'n_layer_sites' in file_add_info['layers'][i]:
                                _n_layer_sites.append(int(file_add_info['layers'][i]['n_layer_sites']))
                                 
                            if 'sites_end_idx_in_device' in file_add_info['layers'][i]:
                                _sites_end_idx_in_device.append(file_add_info['layers'][i]['sites_end_idx_in_device'])
                                 
                            if 'sites_start_idx_in_device' in file_add_info['layers'][i]:
                                _sites_start_idx_in_device.append(file_add_info['layers'][i]['sites_start_idx_in_device'])
                                 
                            if 'thickness' in file_add_info['layers'][0]:
                                _add_info_thickness.append(file_add_info['layers'][i]['thickness'])
                                 
                            if 'x_boundaries' in file_add_info['layers'][0]:
                                _add_info_x_boundaries.append(file_add_info['layers'][i]['x_boundaries'])
                                 
                        sec_lf_add_info.lf_layer_id = _lf_layer_id
                        sec_lf_add_info.n_layer_sites = _n_layer_sites
                        sec_lf_add_info.sites_end_idx_in_device = _sites_end_idx_in_device
                        sec_lf_add_info.sites_start_idx_in_device = _sites_start_idx_in_device
                        sec_lf_add_info.add_info_thickness = _add_info_thickness
                        sec_lf_add_info.add_info_x_boundaries = _add_info_x_boundaries
                if 'run_lf.slr' in file:
                    sec_run_lf_slr = sec_input.m_create(Run_lf_slr)
                    for i, line in enumerate(f):
                        if 'nodes' in line:
                            parts = line.split('=')
                            sec_run_lf_slr.lf_nodes = int(parts[1])
                        if 'ntask' in line:
                            parts = line.split('=')
                            sec_run_lf_slr.lf_ntasks = int(parts[1])
                        if 'mem-per-cpu' in line:
                            parts = line.split('=')
                            sec_run_lf_slr.lf_mem_per_cpu = float(parts[1])
                if 'settings' in file:
                    sec_settings = sec_input.m_create(Settings)
                    materials_section = False   # counter for materials-section in settings-file
                    energies_section = False    # counter for energies-section under materials-section
                    layers_section = False
                    molecule_species_section = False
                    electrodes_section = False
                    pair_input_section = False
                    hole_transfer_integrals_section = False
                    electron_transfer_integrals_section = False
                    dexter_transfer_integrals_section = False
                    qp_output_files_section = False
                    _lf_energies = []
                    _lf_layers_materials = []
                    _lf_molecule_species_material = []
                    _lf_molecule_species_concentration = []
                    _lf_electrodes_workfunction = []
                    _lf_field_direction = []
                    for i, line in enumerate(f):
                        parts = line.split(':')
                        if re.search(r'^#', line):
                            continue
                        if 'pbc' in line:
                            sec_settings.lf_pbc = parts[1]
                            continue
                        if 'excitonics' in line:
                            sec_settings.lf_excitonics = parts[1]
                            continue
                        if 'connect_electrodes' in line:
                            sec_settings.connect_electrodes = parts[1]
                            continue
                        if 'coulomb_mesh' in line:
                            sec_settings.coulomb_mesh = parts[1]
                            continue
                        if re.search(r'\s+holes:', line):
                            sec_settings.particles_holes = parts[1]
                            continue
                        if re.search(r'\s+electrons:', line):
                            sec_settings.particles_electrons = parts[1]
                            continue                        
                        if re.search(r'\s+excitons:', line):
                            sec_settings.particles_excitons = parts[1]
                            continue
                        if 'morphology_width' in line:
                            sec_settings.morphology_width = parts[1]
                            continue
                        if 'materials' in line:
                            materials_section = True
                            continue
                        if 'name' in line and materials_section == True:
                            sec_settings_materials = sec_settings.m_create(Settings_materials)
                            sec_settings_materials.material_name = parts[1]
                            continue
                        if 'input_mode_transport' in line and materials_section == True:
                            _input_mode_transport = ''.join(parts[1:])
                            sec_settings_materials.input_mode_transport = _input_mode_transport
                            continue 
                        if 'exciton preset' in line and materials_section == True:
                            sec_settings_materials.lf_exciton_preset = parts[1]
                            continue
                        if 'molecule_pdb' in line and materials_section == True:
                            sec_settings_materials.lf_molecule_pdb = parts[1]
                            continue
                        if 'qp_output_sigma' in line.lower() and materials_section == True:
                            sec_settings_materials.lf_qp_output_sigma = parts[1]
                            continue
                        if 'qp_output_eaip' in line.lower() and materials_section == True:
                            sec_settings_materials.lf_qp_output_eaip = parts[1]
                            continue
                        if 'qp_output_lambda' in line.lower() and materials_section == True:
                            sec_settings_materials.lf_qp_output_lambda = parts[1]
                            continue
                        if 'energies' in line and materials_section == True:
                            _lf_energies = []
                            energies_section = True
                            continue
                        if re.search(r'\d+,\d+', line) and '-' in line and energies_section == True:
                            _lf_energies.append(line.replace('-', '').replace('[', '').replace(']', '').split(','))
                            sec_settings_materials.lf_energies = _lf_energies
                            continue
                        if ('[' not in line or '-' not in line) and energies_section == True:
                            energies_section = False
                        if re.search(r'^\w', line) and materials_section == True:
                            materials_section = False
                        if 'layers' in line:
                            layers_section = True
                            continue
                        if 'thickness' in line and layers_section == True:    
                            sec_settings_layers = sec_settings.m_create(Settings_layers)
                            sec_settings_layers.layer_thickness = parts[1]
                            continue
                        if 'morphology_input_mode' in line and layers_section == True:
                            sec_settings_layers.layer_morphology_input_mode = parts[1]
                            continue
                        if 'molecule_species' in line and layers_section == True:
                            sec_molecule_species = sec_settings_layers.m_create(Layer_molecule_species)
                            molecule_species_section = True
                            _lf_molecule_species_material = []
                            _lf_molecule_species_concentration = []
                            continue
                        if re.search(r'-\s*material', line) and molecule_species_section == True:    
                            _lf_molecule_species_material.append(parts[1])
                            sec_molecule_species.molecule_species_material = _lf_molecule_species_material
                            continue
                        if 'concentration' in line and molecule_species_section == True:
                            _lf_molecule_species_concentration.append(parts[1])
                            sec_molecule_species.molecule_species_concentration = (
                                _lf_molecule_species_concentration) 
                            continue
                        if (re.search(r'^\w', line) or re.search(r'^-', line) or len(parts)==1) and (
                                molecule_species_section == True):
                            
                            molecule_species_section = False
                        if re.search(r'\w', line) and layers_section == True:
                            layers_section = False
                        if 'neighbours' in line:
                            sec_settings.lf_neighbours = parts[1]
                            continue
                        if 'transfer_integral_source' in line:
                            sec_settings.transfer_integral_source = parts[1]
                            continue
                        if 'electrodes' in line:
                            electrodes_section = True
                            continue
                        if re.search(r'^-', line) and electrodes_section == True:
                            sec_settings_electrodes = sec_settings.m_create(Settings_electrodes)
                        if 'electrode_workfunction' in line and electrodes_section == True:    
                            sec_settings_electrodes.electrode_workfunction = parts[1]
                            continue
                        if 'coupling_model' in line and electrodes_section == True:
                            sec_settings_electrodes.electrode_coupling_model = parts[1]
                            continue
                        if 'electrode_wf_decay_length' in line and electrodes_section == True:
                            sec_settings_electrodes.electrode_wf_decay_length = parts[1]
                            continue
                        if 'electrode_coupling' in line and electrodes_section == True:
                            sec_settings_electrodes.electrode_coupling = parts[1]
                            continue
                        if re.search(r'^\w', line) and electrodes_section == True:
                            electrodes_section = False
                        if 'pair_input' in line:
                            pair_input_section = True
                            continue
                        if re.search(r'^-', line) and pair_input_section == True:
                            sec_settings_pair_input = sec_settings.m_create(Settings_pair_input)
                        if 'molecule 1' in line and pair_input_section == True:
                            sec_settings_pair_input.molecule_1_type = parts[1]
                            continue
                        if 'molecule 2' in line and pair_input_section == True:
                            sec_settings_pair_input.molecule_2_type = parts[1]
                            continue
                        if re.search(r'\sqp_output:', line.lower()) and pair_input_section == True:
                            sec_settings_pair_input.lf_qp_output = parts[1]
                            continue
                        if 'hole_transfer_integrals' in line and pair_input_section == True:
                            sec_hole_transfer_integrals = sec_settings_pair_input.m_create(
                                Settings_hole_transfer_integrals)
                            hole_transfer_integrals_section = True                            
                            continue
                        if 'wf_decay_length' in line and hole_transfer_integrals_section == True:
                            sec_hole_transfer_integrals.hole_transfer_integrals_wf_decay_length = parts[1]
                            continue
                        if 'maximum_ti' in line and hole_transfer_integrals_section == True:
                            sec_hole_transfer_integrals.hole_transfer_integrals_maximum_ti = parts[1]
                            continue
                        if ((re.search(r'electron', line) or re.search(r'dexter', line.lower()) or
                            re.search(r'^-', line) or re.search(r'^\w', line)) and 
                            hole_transfer_integrals_section == True):
                            hole_transfer_integrals_section = False  
                        if 'electron_transfer_integrals' in line and pair_input_section == True:
                            sec_electron_transfer_integrals = sec_settings_pair_input.m_create(
                                Settings_electron_transfer_integrals)
                            electron_transfer_integrals_section = True                            
                            continue
                        if 'wf_decay_length' in line and electron_transfer_integrals_section == True:
                            sec_electron_transfer_integrals.electron_transfer_integrals_wf_decay_length = (
                                parts[1])
                            continue
                        if 'maximum_ti' in line and electron_transfer_integrals_section == True:
                            sec_electron_transfer_integrals.electron_transfer_integrals_maximum_ti = parts[1]
                            continue
                        if ((re.search(r'hole', line) or re.search(r'dexter', line.lower()) or
                            re.search(r'^-', line) or re.search(r'^\w', line)) and 
                            electron_transfer_integrals_section == True):
                            electron_transfer_integrals_section = False
                        if 'dexter_transfer_integrals' in line.lower() and pair_input_section == True:
                            sec_dexter_transfer_integrals = sec_settings_pair_input.m_create(
                                Settings_dexter_transfer_integrals)
                            dexter_transfer_integrals_section = True                            
                            continue
                        if 'wf_decay_length' in line and dexter_transfer_integrals_section == True:
                            sec_dexter_transfer_integrals.dexter_transfer_integrals_wf_decay_length = parts[1]
                            continue
                        if 'maximum_ti' in line and dexter_transfer_integrals_section == True:
                            sec_dexter_transfer_integrals.dexter_transfer_integrals_maximum_ti = parts[1]
                            continue
                        if ((re.search(r'hole', line) or re.search(r'electron', line.lower()) or
                            re.search(r'^-', line) or re.search(r'^\w', line)) and 
                            dexter_transfer_integrals_section == True):
                            dexter_transfer_integrals_section = False
                        if re.search(r'^\w', line) and pair_input_section == True:
                            pair_input_section = False
                        if 'simulations' in line:
                            sec_settings.lf_simulations = int(parts[1])
                            continue
                        if 'measurement' in line:
                            sec_settings.lf_measurement = parts[1]
                            continue
                        if 'temperature' in line.lower():
                            sec_settings.lf_temperature = parts[1]
                            continue
                        if 'field_direction' in line:
                            _lf_field_direction = parts[1].replace('[', '').replace(']', '').replace(',', '').split()
                            sec_settings.lf_field_direction = _lf_field_direction
                            continue
                        if 'field_strength' in line:
                            _fields = parts[1].split()
                            sec_settings.lf_field_strength = _fields
                            continue
                        if 'initial_holes' in line:
                            sec_settings.lf_initial_holes = int(parts[1])
                            continue
                        if 'initial_electrons' in line:
                            sec_settings.lf_initial_electrons = int(parts[1])    
                            continue
                        if 'iv_fluctuation' in line:
                            sec_settings.lf_iv_fluctuation = parts[1]
                            continue
                        if 'max_iterations' in line:
                            sec_settings.lf_max_iterations = int(parts[1])
                            continue
                        if 'ti_prune' in line:
                            sec_settings.lf_ti_prune = parts[1]
                            continue
                        if 'noise_damping' in line:
                            sec_settings.lf_noise_damping = parts[1]
                            continue
                        if 'expansion_scheme' in line:
                            sec_settings.lf_expansion_scheme = parts[1]
                            continue
                        if 'qp_output_files' in line.lower():
                            qp_output_files_section = True
                            continue
                        if re.search(r'^-', line) and qp_output_files_section == True:
                            sec_qp_output_files = sec_settings.m_create(Settings_qp_output_files)
                        if 'name' in line and qp_output_files_section == True:
                            sec_qp_output_files.qp_output_files_name = parts[1]
                            continue
                        if 'qp_output.zip' in line.lower() and qp_output_files_section == True:
                            sec_qp_output_files.qp_output_files_output_zip = parts[1]
                            continue
                        if re.search(r'^\w', line) and qp_output_files_section == True:
                            qp_output_files_section = False
                        if re.search(r'^rates', line):
                            sec_settings.lf_rates = parts[1]
                            continue
                        if 'superexchange' in line:
                            sec_settings.lf_superexchange = parts[1]
                            continue
                        if 'epsilon_material' in line:
                            sec_settings.lf_epsilon_material = parts[1]
                            continue
                if re.search(r'molecule.pdb', file):
                    sec_molecule_pdb = sec_input.m_create(LF_molecule_pdb_file)
                    _lf_residue_class = []
                    _lf_atom_serial_number = []
                    _lf_atom_name = []
                    _lf_residue_type = []
                    _lf_chain_identifier = []
                    _lf_residue_sequence_number = []
                    _lf_molecule_pdb_coordinates = []
                    _lf_molecule_pdb_occupancy = []
                    _lf_molecule_pdb_temperature = []
                    _lf_molecule_pdb_element = []
                    _lf_molecule_pdb_charge = []
                    for i, line in enumerate(f):
                        if 'TER' in line:
                            continue
                        parts = line.split()
                        _lf_residue_class.append(parts[0])
                        _lf_atom_serial_number.append(int(parts[1]))
                        _lf_atom_name.append(parts[2])
                        _lf_residue_type.append(parts[3])
                        _lf_chain_identifier.append(parts[4])
                        _lf_residue_sequence_number.append(parts[5])
                        _lf_molecule_pdb_coordinates.append(parts[6:9])
                        _lf_molecule_pdb_occupancy.append(parts[9])
                        _lf_molecule_pdb_temperature.append(parts[10])
                        _lf_molecule_pdb_element.append(parts[11])
                        _lf_molecule_pdb_charge.append(parts[12])
                    sec_molecule_pdb.lf_residue_class = _lf_residue_class
                    sec_molecule_pdb.lf_atom_serial_number = _lf_atom_serial_number
                    sec_molecule_pdb.lf_atom_name = _lf_atom_name
                    sec_molecule_pdb.lf_residue_type = _lf_residue_type
                    sec_molecule_pdb.lf_chain_identifier = _lf_chain_identifier
                    sec_molecule_pdb.lf_residue_sequence_number = _lf_residue_sequence_number
                    sec_molecule_pdb.lf_molecule_pdb_coordinates = _lf_molecule_pdb_coordinates
                    sec_molecule_pdb.lf_molecule_pdb_occupancy = _lf_molecule_pdb_occupancy
                    sec_molecule_pdb.lf_molecule_pdb_temperature = _lf_molecule_pdb_temperature
                    sec_molecule_pdb.lf_molecule_pdb_element = _lf_molecule_pdb_element
                    sec_molecule_pdb.lf_molecule_pdb_charge = _lf_molecule_pdb_charge

                if 'files_for_kmc' in root:
                    if not files_for_kmc_hasrun:
                        sec_files_for_kmc = sec_input.m_create(Files_for_kmc)
                        files_for_kmc_hasrun = True
                    if 'COM' in file:
                        _COM = []
                        for i, line in enumerate(f):
                            _a = [float(x) for x in list(line.split())]
                            _COM.append(_a)
                        sec_files_for_kmc.lf_COM = _COM
                    if 'inner_idxs' in file:
                        _lf_inner_idxs = []
                        for i, line in enumerate(f):
                            _lf_inner_idxs.append(float(line))
                        sec_files_for_kmc.lf_inner_idxs = _lf_inner_idxs
                    if 'ip_ea' in file.lower():
                        _lf_ip_ea = []
                        for i, line in enumerate(f):
                            _a = [float(x) for x in list(line.split())]
                            _lf_ip_ea.append(_a)
                        sec_files_for_kmc.lf_ip_ea = _lf_ip_ea
                    if re.search(r'js_homo_mol_pairs_\d+', file.lower()):
                        sec_js_homo_mol_pairs = sec_files_for_kmc.m_create(Js_homo_mol_pairs)
                        _js_homo_mol_pairs_value = []
                        for i, line in enumerate(f):
                            _a = [float(x) for x in list(line.split())]
                            _js_homo_mol_pairs_value.append(_a)
                        sec_js_homo_mol_pairs.js_homo_mol_pairs_value = _js_homo_mol_pairs_value
                    
                    if re.search(r'js_lumo_mol_pairs_\d+', file.lower()):
                        sec_js_lumo_mol_pairs = sec_files_for_kmc.m_create(Js_lumo_mol_pairs)
                        _js_lumo_mol_pairs_value = []
                        for i, line in enumerate(f):
                            _a = [float(x) for x in list(line.split())]
                            _js_lumo_mol_pairs_value.append(_a)
                        sec_js_lumo_mol_pairs.js_lumo_mol_pairs_value = _js_lumo_mol_pairs_value
                    if re.search(r'js_dexter_mol_pairs_\d+', file.lower()):
                        sec_js_dexter_mol_pairs = sec_files_for_kmc.m_create(Js_dexter_mol_pairs)
                        _js_dexter_mol_pairs_value = []
                        for i, line in enumerate(f):
                            _a = [float(x) for x in list(line.split())]
                            _js_dexter_mol_pairs_value.append(_a)
                        sec_js_dexter_mol_pairs.js_dexter_mol_pairs_value = _js_dexter_mol_pairs_value
                    
                    if re.search(r'$sigma_mol_pairs_\d+', file.lower()):
                        sec_sigma_mol_pairs = sec_files_for_kmc.m_create(Sigma_mol_pairs)
                        _sigma_mol_pairs_value = []
                        for i, line in enumerate(f):
                            _sigma_mol_pairs_value.append(float(line))
                        sec_sigma_mol_pairs.sigma_mol_pairs_value = _sigma_mol_pairs_value
                
                if 'vacuum_lambda' in root:
                    sec_lf_vacuum_lambda = sec_input.m_create(LF_vacuum_lambda)
                    if 'lambda' in file:
                        for i, line in enumerate(f):
                            if 'hole:' in line.lower():
                                parts = line.split(':')
                                sec_lf_vacuum_lambda.lf_vacuum_lambda_hole = float(parts[1]) 
                            if 'electron:' in line.lower():
                                parts = line.split(':')
                                sec_lf_vacuum_lambda.lf_vacuum_lambda_electron = float(parts[1])
                
                if 'runtime_data' in root:
                    if not runtime_data_hasrun:
                        sec_runtime_data = sec_lightforge_data.m_create(Runtime_data)
                        runtime_data_hasrun = True
                    if re.search(r'experiment_inventory_\d+', file):
                        sec_experiment_inventory = sec_runtime_data.m_create(LF_experiment_inventory)
                        experiment_inventory_field_direction_section = False
                        experiment_inventory_initial_charges_section = False
                        _lf_experiment_inventory_field_direction = []
                        _lf_experiment_inventory_initial_charges = []
                        for i, line in enumerate(f):
                            if 'field_direction' in line:
                                experiment_inventory_field_direction_section = True
                                continue
                            if re.search(r'^-', line) and experiment_inventory_field_direction_section == True:
                                parts = line.replace(' ', '').split('-')
                                _lf_experiment_inventory_field_direction.append(parts[1])
                                continue
                            if re.search(r'^\w', line) and experiment_inventory_field_direction_section == True:
                                sec_experiment_inventory.lf_experiment_inventory_field_direction = _lf_experiment_inventory_field_direction
                                experiment_inventory_field_direction_section = False
                            if 'field_strength' in line:
                                parts = line.split(':')
                                sec_experiment_inventory.lf_experiment_inventory_field_strength = parts[1]
                                continue
                            if 'initial_charges' in line:
                                experiment_inventory_initial_charges_section = True
                                continue
                            if re.search(r'^-', line) and experiment_inventory_initial_charges_section == True:
                                parts = line.replace(' ', '').split('-')
                                _lf_experiment_inventory_initial_charges.append(parts[1])
                                continue
                            if re.search(r'^\w', line) and experiment_inventory_initial_charges_section == True:
                                sec_experiment_inventory.lf_experiment_inventory_initial_charges = _lf_experiment_inventory_initial_charges
                                experiment_inventory_initial_charges_section = False
                            if 'job_id' in line:
                                parts = line.split(':')
                                sec_experiment_inventory.lf_experiment_inventory_job_id = int(parts[1])
                                continue
                            if 'mat_id' in line:
                                parts = line.split(':')
                                sec_experiment_inventory.lf_experiment_inventory_mat_id = int(parts[1])
                                continue
                            if 'settings_id' in line:
                                parts = line.split(':')
                                sec_experiment_inventory.lf_experiment_inventory_settings_id = int(parts[1])
                                continue
                    if re.search(r'particle_positions_\d+', file):
                        sec_particle_positions = sec_runtime_data.m_create(LF_particle_positions)
                        _lf_particle_positions_value = []
                        for i, line in enumerate(f):
                            parts = line.split()
                            _lf_particle_positions_value.append(parts)
                        sec_particle_positions.lf_particle_positions_value = _lf_particle_positions_value 
    '''
class NewParser(MatchingParser):

    def parse(
        self, 
        mainfile: str, 
        archive: 'EntryArchive',
        logger: 'BoundLogger', 
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None: 
        logger.info('NewParser.parse', parameter=configuration.parameter)
        
        sec_program = archive.m_setdefault('run.program')
        sec_program.name = "Lightforge"
        sec_workflow = archive.m_setdefault('workflow')
        sec_workflow.type = 'single_point'

        mainfile = Path(mainfile)
#        archive.workflow2 = Workflow(name='test')
#        sec_program = archive.m_setdefault('run.program')
#        sec_program.name = "Lightforge test test"
        
        DetailedParser(mainfile, archive)