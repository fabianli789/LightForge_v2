from nomad.config.models.plugins import NormalizerEntryPoint
from pydantic import Field


class NewNormalizerEntryPoint(NormalizerEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from lightforge_v2.normalizers.normalizer import NewNormalizer

        return NewNormalizer(**self.dict())


normalizer_entry_point = NewNormalizerEntryPoint(
    name='NewNormalizer',
    description='New normalizer entry point configuration.',
)
