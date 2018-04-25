%Class definition used for unit conversion in JATec
%
%Creation: 8 Jan 2013 - Jeff Anderson
%NOTE!! When adding conversion I parse the original units via the 2, so
%that will mess up if you have mps^2 for instance..so use just double s's.  Example mpss2g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef myConstants
    properties (Constant = true)
        g = 9.80665;
        kph2mps = 0.277778;
        kph2mph = 1/1.60934;
        mps2kph = 3.6;
        mph2mps = 0.44704;
        mph2kph = 1.60934;
        deg2rad = pi/180;
        g2mps2  = 9.80665;
        mpss2g  = 1/9.80665;  
        rad2deg = 180/pi;
        lbf2n = 4.44822;
        n2lbf = 1/4.44822;
        pa2mpa = 1/10^6;
        mpa2pa = 10^6;
        kpa2bar = 1/100;
        kpa2psi = 1/(6894.75729/1000);
        bar2kpa = 100;
        psi2pa = 6894.75729;
        pa2psi = 1/6894.75729;
        psi2kpa = 6894.75729/1000;
        psi2bar = 0.0689475729;
        bar2psi = 1/0.0689475729;
        ft2m = 0.3048;
        ftlbf2nm = 0.3048*4.44822;   %WRONG
        nm2ftlbf = 1/(0.3048*4.44822); %WRONG
        rpm2radps = 2*pi/60;
        radps2rpm = 1/(2*pi/60);
        radps2hz = 1/(2*pi);
        hz2radps = 2*pi;
        kg2lb = 2.20462262;
        in2m = 0.0254;
        m2in = 1/0.0254;
        mm2in = 1/25.4;
        m2mi = 0.000621371;
        mi2m = 1/0.000621371;
        in2mm = 25.4;
        kg2N = 9.81;
        cm2m = 0.01;
        mm2m = 0.001;
        npdeg2nprad = 180/pi;
        nmpdeg2nmprad = 180/pi;
		mmps2kph = 0.0036;
        km2mi = 0.621371;
        mi2km = 1/0.621371;
        w2hp = 0.00134102;
        kW2hp = 1.34102;
        hp2w = 1/0.00134102;
        hp2kW = 1/1.34102;
    end
    

end