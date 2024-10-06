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
from nomad.datamodel.metainfo.annotations import ELNAnnotation, ELNComponentEnum

from nomad.metainfo import Quantity, SchemaPackage, Section, MSection, SubSection   # changed
from runschema.run import Run
from runschema.calculation import Calculation 

configuration = config.get_plugin_entry_point(
    'lightforge_v2.schema_packages:schema_package_entry_point'
)

m_package = SchemaPackage()


class Currents(MSection):
    current_density = Quantity(type=np.float64)

class LightforgeCalculation(Calculation):
    currents = SubSection(sub_section=Currents, repeats=False)

class LightforgeRun(Run):
    lightforge_calculation = SubSection(sub_section=LightforgeCalculation, repeats=False)


class NewSchemaPackage(Schema):
    name = Quantity(
        type=str, a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity)
    )
    message = Quantity(type=str)

    def normalize(self, archive: 'EntryArchive', logger: 'BoundLogger') -> None:
        super().normalize(archive, logger)

        logger.info('NewSchema.normalize', parameter=configuration.parameter)
        self.message = f'Hello {self.name}!'


m_package.__init_metainfo__()
