%% Create adigator files
daq = generateInitialDaq();
setup = daq.header.setup;
setup.functions.continuous        = @fourwheelContinuous;
setup.functions.endpoint          = @fourwheelEndpoint;
adigatorGenFiles4gpops2(setup);