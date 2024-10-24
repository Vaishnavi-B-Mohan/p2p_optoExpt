classdef ionic_channel_hemond
    methods(Static)
        function params = SetParams(ion)
            %% Defining conductances for different ionic channels in a neuron based on Hemond et al., 2008
            gNa = 22;
            gKdr = 10;
            gH = 0.01;
            gCaL = 0.01;
            gKa = 20;
            gKM = 0.5;
            gL = 0.04;
            EH = -30;
            ENa = 55;
            EK = -90;
            EL = -70;
            tauCa = 100;
            IDC = 0;
            Cm = 1.41;

        end
    end
end