classdef opsin_photocurrent
    properties
        dt = 0.1;
    end
    methods(Static)

        function params = SetParams(opsin)
            % Defining the parameters of the model for different opins
            params.name = opsin;
            switch opsin
                case 'ChR2'
                    params.Gd1 = 0.09;
                    params.Gd2 = 0.01;
                    params.Gr = 0.5e-3;
                    params.g0 = 0.12e-2; % this varies for different neurons: hippocampal, RGNs
                    params.g0_ph = 5.9;
                    params.phim = 4e16;
                    params.k1 = 3;
                    params.k2 = 0.18;
                    params.Gf0 = 0.015;
                    params.Gb0 = 0.005;
                    params.kf = 0.03;
                    params.kb = 0.003;
                    params.gamma = 0.05;
                    params.p = 1;
                    params.q = 1;
                    params.E = 0;

                case 'ReaChR'
                    params.Gd1 = 7.7e-3;
                    params.Gd2 = 1.25e-3;
                    params.Gr = 3.33e-5;
                    params.g0 = 0.28e-2; % this varies for different neurons: hippocampal, RGNs
                    params.g0_ph = 14.28;
                    params.phim = 5e17; %5e21;
                    params.k1 = 1.2;
                    params.k2 = 0.01;
                    params.Gf0 = 0.0005;
                    params.Gb0 = 0.0005;
                    params.kf = 0.012;
                    params.kb = 0.001;
                    params.gamma = 0.05;
                    params.p = 1;
                    params.q = 1;
                    params.E = 7;

                case 'ChrimsonR'
                    params.Gd1 = 0.067;
                    params.Gd2 = 0.01;
                    params.Gr = 0.5e-3;
                    params.g0 = 0.24e-2; % this varies for different neurons: hippocampal, RGNs
                    params.g0_ph = 12.25;
                    params.phim = 20e17; %20e20; 
                    params.k1 = 6;
                    params.k2 = 0.1;
                    params.Gf0 = 0.02;
                    params.Gb0 = 0.05;
                    params.kf = 0.1;
                    params.kb = 0.001;
                    params.gamma = 0.05;
                    params.p = 0.6;
                    params.q = 1;
                    params.E = 0;

                case 'CsChrimson'
                    params.Gd1 = 0.033;
                    params.Gd2 = 0.017;
                    params.Gr = 5e-6;
                    params.g0 = 0.37e-2; % this varies for different neurons: hippocampal, RGNs
                    params.g0_ph = 18.48;
                    params.phim = 6e16;
                    params.k1 = 3;
                    params.k2 = 0.04;
                    params.Gf0 = 0.005;
                    params.Gb0 = 0.01;
                    params.kf = 0.01;
                    params.kb = 0.6;
                    params.gamma = 0.05;
                    params.p = 1;
                    params.q = 1;
                    params.E = -10;

                case 'bReaChES'
                    params.Gd1 = 0.025;
                    params.Gd2 = 0.01;
                    params.Gr = 3.3e-5;
                    params.g0 = 0.73e-2; % this varies for different neurons: hippocampal, RGNs
                    params.g0_ph = 36.5;
                    params.phim = 6e15; % 6e17; 
                    params.k1 = 0.4;
                    params.k2 = 0.01;
                    params.Gf0 = 0.002;
                    params.Gb0 = 0.002;
                    params.kf = 0.01;
                    params.kb = 0.04;
                    params.gamma = 0.05;
                    params.p = 1;
                    params.q = 1;
                    params.E = 10;

                case 'ChRmine'
                    params.Gd1 = 0.02; % in ms-1
                    params.Gd2 = 0.013; % in ms-1
                    params.Gr = 5.9e-4; % in ms-1
                    params.g0 = 2.2e-2; % this varies for different neurons: hippocampal, RGNs in nS/sq .cm
                    params.g0_ph = 110; % in nS
                    params.phim = 2.1e15;%2.1e15; % in photons/sq. mm/s
                    params.k1 = 0.2; % in ms-1
                    params.k2 = 0.01; % in ms-1
                    params.Gf0 = 0.0027; % in ms-1
                    params.Gb0 = 0.0005; % in ms-1
                    params.kf = 0.001; % in ms-1
                    params.kb = 0; % in ms-1 
                    params.gamma = 0.05;
                    params.p = 0.8;
                    params.q = 1;
                    params.E = 5.64; % in mV    
            end
        end

    end
    methods
        function [I_opsin] = get_opsin_current(this, opsin, V, lambda, Irr)
            %% Opsin model equations - four state photocurrent model 
            % dC1/dt = Gd1 * O1 + Gr * C2 - Ga1(phi) * C1 
            % dO1/dt = Ga1(phi) * C1 + Gb(phi) * O2 - (Gd1 + Gr(phi)) * O1
            % dO2/dt  = Ga2(phi) * C2 + Gf(phi) * O1 - (Gd2 + Gb(phi)) * O2
            % dC2/dt = Gd2 * O2 - (Gr + Ga2(phi)) * C2
            % where, C1 + O1 + O2 + C2 = 1


            if ~exist('lambda', 'var')
                    lambda = 590; % in nm
            end
            if ~exist('Irr', 'var')
                Irr = [ones(1, 1500), zeros(1, 500)]; % irradiance is determined by the background lighting condition
            end
            if strcmp(opsin, 'ChR2')
                lambda = 460;
            end
%             time = 0:this.dt:2000;
            NumPts = size(Irr,1);
            params = this.SetParams(opsin);
            h = 6.62607015e-17; % Planck's constant after cancelling some decimal pts with labmda and c
            c = 2.99792458; % speed of light (after cancelling some decimal points with Planck's)
            epsilon = 1;
            % defining the initial conditions: setting the proportion of closed state C1 to 1 and assuming all open states are 0
            pts = ones(size(Irr));
            C1 = 1*pts;
            C2 = 0*pts;
            O1 = 0*pts;
            O2 = 0*pts;
            f_phi(1,:) = O1(1,:) + params.gamma * O2(1,:);

            for i = 1:NumPts-1

                phi = (lambda*Irr(i,:))/(h*c);
                Ga1 = (epsilon * params.k1 * (phi.^params.p)) ./ (phi.^params.p + params.phim.^params.p);
                Ga2 = (epsilon * params.k2 * (phi.^params.p)) ./ (phi.^params.p + params.phim.^params.p);
                Gf = params.Gf0 + (epsilon * params.kf *(phi.^params.q)) ./ (phi.^params.q + params.phim.^params.q);
                Gb = params.Gb0 + (epsilon * params.kb *(phi.^params.q)) ./ (phi.^params.q + params.phim.^params.q);

                C1(i+1,:) = (params.Gd1 .* O1(i,:) + params.Gr .* C2(i,:) - Ga1 .* C1(i,:))*this.dt + C1(i,:);
                O1(i+1,:) = (Ga1 .* C1(i,:) + Gb .* O2(i,:) - (params.Gd1 + Gf) .* O1(i,:))*this.dt + O1(i,:);
                O2(i+1,:) = (Ga2 .* C2(i) + Gf .* O1(i) - (params.Gd2 + Gb) .* O2(i))*this.dt + O2(i,:);
                C2(i+1,:) = (params.Gd2 .* O2(i,:) - (params.Gr + Ga2) .* C2(i,:))*this.dt + C2(i,:);
                f_phi(i+1,:) = O1(i+1,:) + params.gamma .* O2(i+1,:);
            end
%                 phi(end+1) = 0;
                g_opsin = params.g0*f_phi;
%                 figure(3)
%                 plot(time, f_phi);
%                 title('Photon flux with time');
%                 xlabel('Time in ms');
%                 grid on;
%                 hold on;
                I_opsin = g_opsin * (V - params.E);
        end

    end
end