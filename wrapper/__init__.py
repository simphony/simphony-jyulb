from simphony.engine import ABCEngineExtension
from simphony.engine import EngineInterface
from simphony.engine.decorators import register

from .jyulb_engine import Wrapper as InternalWrapper

__all__ = ["InternalWrapper"]


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
        # Flow type, relaxation model etc.
        # jyulb_features =\
        #     self.create_engine_metadata_feature(LAMINAR_FLOW, TRT)

        jyulb_features = None

        jyulb = self.create_engine_metadata('JYU-LB',
                                            jyulb_features,
                                            [EngineInterface.Internal])
        return [jyulb]

    def create_wrapper(self, cuds, engine_name, engine_interface):
        """Create a wrapper to the requested engine.

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
        use_internal_interface = True
        if engine_interface == EngineInterface.FileIO:
            use_internal_interface = False

        if engine_name == 'JYU-LB':
            if use_internal_interface:
                return InternalWrapper(cuds=cuds)
            else:
                raise Exception('File-IO wrapper for the JYU-LB engine '
                                'is not supported.')
        else:
            raise Exception('Only JYU-LB engine is supported. '
                            'Unsupported engine: %s', engine_name)
