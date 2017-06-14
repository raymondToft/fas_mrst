Models and examples that implement dual-porosity models in the AD framework.

Description of folders and files:

- models:
--- DualPorosityReservoirModel: Base class for all reservoirs that show dual-porosity behaviour
--- ThreePhaseBlackOilDPModel:  Three-phase black-oil model. This serves as a base class for the two phase model. 
							    However, the equations for the 3ph case are still not implemented and will come 
								in a future release. *****Do not instantiate this model directly, as it will give 
								an error because the equations function does not exist.*****
--- TwoPhaseOilWaterDPModel: Two-phase compressible model. This is ready to be used and the basic usage is shown
							 in an example.							
--- SequentialPressureTransportDPModel: A sequential model for the two-phase compressible model. This is a
										dual-porosity implementation of the SequentialPressureTransportModel in the
										blackoil-sequential module.
--- TransportOilWaterDPModel: A model that represents the transport equation in the sequential model.
--- PressureOilWaterWaterDPModel: A model that represents the pressure equation in the sequential model.

- models\equations:
--- equationsOilWaterDP: Equation for the TwoPhaseOilWaterDPModel.
--- transportEquationOilWaterDP: Transport equation for the TransportOilWaterDPModel.
--- pressureEquationOilWaterDP: Pressure equation for the PressureOilWaterWaterDPModel.

- utils:
--- getSequentialDPModelFromFI: useful function to transform a fully coupled model (TwoPhaseOilWaterDPModel) into
								a sequential one (SequentialPressureTransportDPModel).
								
- examples:
---	quarterSpotDP: a simple quarter-spot model using TwoPhaseOilWaterDPModel.
--- quarterSpotDPSequential: a simple quarter-spot model using SequentialPressureTransportDPModel.

References:
[1] Kazemi et al. Numerical Simulation of Water-Oil Flow in Naturally Fractured Reservoirs, SPE Journal, 1976
[2] Quandalle and Sabathier. Typical Features of a Multipurpose Reservoir Simulator, SPE Journal, 1989
[3] Lu et al. General Transfer Functions for Multiphase Flow in Fractured Reservoirs, SPE Journal, 2008

