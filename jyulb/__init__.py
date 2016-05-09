from simphony.engine import ABCEngineExtension
from simphony.engine import EngineInterface
from simphony.engine.decorators import register

from .internal.isothermal import JYULBEngine as InternalWrapper
from .fileio.isothermal import JYULBEngine as FileIOWrapper
from .cuba_extension import CUBAExtension

__all__ = ["InternalWrapper", "FileioWrapper", "CUBAExtension"]


@register
class JYULBExtension(ABCEngineExtension):
    """JYU-LB extension.

    This extension provides support for JYU-LB engine.
    """

    def get_supported_engines(self):
        """Get metadata about supported engines.

        Returns
        -------
        list: a list of EngineMetadata objects
        """
        # TODO: Add proper features as soon as the metadata classes are ready.
        # lammps_features =\
        #     self.create_engine_metadata_feature(MolecularDynamics(),
        #                                         [Verlet()])
        jyulb_features = None

        jyulb = self.create_engine_metadata('JYULB',
                                            jyulb_features,
                                            [EngineInterface.Internal,
                                             EngineInterface.FileIO])

        return [jyulb]

    def create_wrapper(self, cuds, engine_name, engine_interface):
        """Creates a wrapper to the requested engine.

        Parameters
        ----------
        cuds: CUDS
          CUDS computational model data
        engine_name: str
          name of the engine, must be supported by this extension
        engine_interface: EngineInterface
          the interface to interact with engine

        Returns
        -------
        ABCEngineExtension: A wrapper configured with cuds and ready to run
        """
        use_internal_interface = False
        if engine_interface == EngineInterface.Internal:
            use_internal_interface = True

        if engine_name == 'JYULB':
            if use_internal_interface:
                return InternalWrapper(cuds=cuds)
            else:
                return FileIOWrapper(cuds=cuds)
        else:
            raise Exception('Only JYULB engine is supported. '
                            'Unsupported engine: %s', engine_name)
