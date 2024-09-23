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

#from nomad.datamodel.metainfo.simulation.run import Run                        # no import error but error at sec_run = archive.m_create(Run)
#from nomad.datamodel.metainfo.simulation.calculation import Calculation

#from runschema.run import Run, Program                                         # import error: no module named 'runschema'
#from nomad_simulations.schema_packages.general import Program, Simulation      # import error: no module 'nomad_simulations'

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
#        sec_program = archive.m_setdefault('run.program')
#        sec_program.name = "Lightforge test test"
    
#        sec_run = archive.m_create(Run)
#        sec_calc = Calculation()
#        sec_calc.pressure = 5