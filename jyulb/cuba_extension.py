from enum import Enum, unique


@unique
class CUBAExtension(Enum):
    """ Provisional CUBA keywords specific for JYU-LB
    These are additional CUBA-Keywords that are not included
    in simphony-common yet. The proposed description for
    these new CUBA keywords is:
    - description: Collision operator of a model Boltzmann equation
    domain: [LBM]
    key: COLLISION_OPERATOR
    name: CollisionOperator
    number: 200
    shape: [1]
    type: integer
    - description: Reference value for density
    domain: [LBM]
    key: REFERENCE_DENSITY
    name: ReferenceDensity
    number: 201
    shape: [1]
    type: double
    - description: Toggle for an external force
    domain: [LBM]
    key: EXTERNAL_FORCING
    name: ExternalForcing
    number: 202
    shape: [1]
    type: boolean
    - description: Classification of a fluid flow
    domain: [LBM]
    key: FLOW_TYPE
    name: FlowType
    number: 203
    shape: [1]
    type: integer
    - description: Acceleration due to gravity
    domain: [LBM]
    key: GRAVITY
    name: Gravity
    number: 204
    shape: [3]
    type: double
    """

    COLLISION_OPERATOR = "COLLISION_OPERATOR"
    REFERENCE_DENSITY = "REFERENCE_DENSITY"
    EXTERNAL_FORCING = "EXTERNAL_FORCING"
    FLOW_TYPE = "FLOW_TYPE"
    GRAVITY = "GRAVITY"