%% Create adigator files
daq = generateInitialDaq('loadGuess',false);
setup = daq.header.setup;
setup.functions.continuous        = @fourwheelContinuous;
setup.functions.endpoint          = @fourwheelEndpoint;
adigatorGenFiles4gpops2(setup);