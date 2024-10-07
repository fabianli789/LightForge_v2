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

#from lightforge_v2.schema_packages.schema_package import NewSchemaPackage
from lightforge_v2.schema_packages.schema_package import Currents, LightforgeCalculation, LightforgeRun


configuration = config.get_plugin_entry_point(
    'lightforge_v2.parsers:parser_entry_point'
)


class NewParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None:
        logger.info('NewParser.parse', parameter=configuration.parameter)

        archive.workflow2 = Workflow(name='test')
        sec_program = archive.m_setdefault('run.program')
        sec_program.name = "Lightforge test test"
              
#        sec_run = archive.m_create(Run)
#        sec_calc = sec_run.m_create(LightforgeCalculation) 
#        sec_currents = sec_calc.m_create(Currents)

#        sec_currents.current_density = 6
#        sec_calc.pressure = 5

        run = Run()
        archive.run.append(run)
        currents = Currents()
        
        calculation = LightforgeCalculation(currents=currents)
        currents.current_density = 9
        run.calculation.append(calculation)