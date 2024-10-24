classdef ionic_channel_fm
    properties
        dt = 0.001;
        T = 295.15; % Temperature in K or 22 degrees C;
    end

    methods(Static)
        function params = SetParams(ion, V)
            params.name = ion;
            switch ion
                case 'Na'
                    params.g = 50;
                    params.E = 35;
                    params.p = 3;
                    params.alpha_on = (-0.6*(V+30)) / (exp(-(V+30)/10) - 1);
                    params.beta_on = 20 * exp((-V - 55) / (18));
                    params.q = 1;
                    params.alpha_off = 0.4 * exp((-V - 50) / (20));
                    params.beta_off = 6 / (exp(-(V+20)/10) + 1);

                case 'K'
                    params.g = 12;
                    params.E = -75; % This is given in the paper explicitly. But this is the inactivating potassium channel current.
                    params.p = 4;
                    params.alpha_on = (-0.02*(V+40)) / (exp(-(V+40)/10) - 1);
                    params.beta_on = 0.4 * exp((-V - 50) / (80));
                    params.q = 0;
                    params.alpha_off = 1;
                    params.beta_off = 1;

                case 'Ca'
                    params.g = 2.2;
                    params.p = 3;
                    params.alpha_on = (-0.3*(V+13)) / (exp(-(V+13)/10) - 1);
                    params.beta_on = 10 * exp((-V - 38) / (18));
                    params.q = 0;
                    params.alpha_off = 1;
                    params.beta_off = 1;

                case 'Ka'
                    params.g = 36;
                    params.E = -75; % This is not given in the paper explicitly. But this is the activating potassium channel current.
                    params.p = 3;
                    params.alpha_on = (-0.006*(V+90)) / (exp(-(V+90)/10) - 1);
                    params.beta_on = 0.1 * exp((-V - 30) / (10));
                    params.q = 1;
                    params.alpha_off = 0.04 * exp((-V - 70) / (20));
                    params.beta_off = 0.6 / (exp(-(V+40)/10) + 1);

                case 'KCa'
                    params.g = 0.05;
                    params.p = 0;
                    params.q = 0;
                    params.alpha_on = 1; params.alpha_off = 1;
                    params.beta_on = 1; params.beta_off = 1;

                case 'L'
                    params.g = 0.147;
                    params.E = -61;
                    params.p = 0;
                    params.q = 0;
                    params.alpha_on = 1; params.alpha_off = 1;
                    params.beta_on = 1; params.beta_off = 1;
            end
        end
    end

    methods
        function [gv_on, gv_off] = get_gating_variable(this, params)
            time = [0: this.dt: 1];

            if params.p == 0
                gv_on = 1;
            else
                gv_on(1) = 0;
            end
            if params.q == 0
                gv_off = 1;
            else
                gv_off(1) = 0;
            end

            
            

            if params.p ~= 0
                for i = 1:length(time)-1
                    gv_on(i+1) = (1 - this.dt*(params.alpha_on + params.beta_on))*gv_on(i) + this.dt*params.alpha_on;
                    if params.q ~= 0
                        
                        gv_off(i+1) = (1 - this.dt*(params.alpha_off + params.beta_off))*gv_off(i) + this.dt*params.alpha_off;
                    end
                end
            end
        end


        function [ICa, ECa, CCa_int] = get_calcium_current(this, V)
            % Calcium current follows the same equation as other ion
            % channels, i.e., If = gf * m^p * h^q * (V-Ef)
            % But, ECa follows the Nernst equation:
            % ECa = (RT/2F)*ln((Ca2+)ext / (Ca2+)int (t))
            % where,
            % R: is the Universal Gas Constant
            % T: Temperature in K
            % F: Faraday's constant which is the charge held in 1M of
            % electrons
            % (Ca2+)ext: extracellular calcium ion concentration
            % (Ca2+)int (t): intracellular calcium ion concentration which
            % is a function of time

            % (Ca2+)int (t) follows the first order diff. eqn:
            % d/dt ((Ca2+)int) = -(3/2Fr) ICa -  ((Ca2+)int - (Ca2+)res)/Tau_Ca
            % where,
            % 3/2Fr = 0.000015
            % (Ca2+)res = 1e-4 mM
            %

            params = this.SetParams('Ca', V);

            % Defining constants
            R = 8.3144626; % Universal gas constant in J/K/mol
            F = 9.6485332; % Faraday constant in  C/mol
            c1 = R*this.T/(2*F);
            time = [0: this.dt: 1];
            nt = length(time);

            CCa_ext = 1.8; % extracellular Ca2+ concentration = 1.8 mM - this is a constant?
            CCa_int = 1e-4 * ones(1,nt); % initial intracellular Ca2+ concentration = 1e-4 mM - this varies with calcium current

            % Initialization
            ICa = zeros(1,nt); m = zeros(1,nt);
            ECa = (c1 * log(CCa_ext/CCa_int(1))) * ones(1, nt);

            for i = 2: nt
                CCa_int(i) = 0.015e-3 * this.dt * ICa(i-1) + (1 - this.dt/0.05) * CCa_int(i-1) + (this.dt/0.05 * 1e-4); % Intracellular Ca2+ concentration is a first-order kinetic equation which is a function of Calcium ion current
                if CCa_int(i) < 0 
                    print("There's something wrong here!")
                end
                ECa(i) = c1 * log(CCa_ext/CCa_int(i)); % Resting potential follows the Nernst equation
                m(i) = (1 - this.dt*(params.alpha_on + params.beta_on))*m(i-1) + this.dt*params.alpha_on; % Gating variable - first order kinetic eqn. like the other ions
                ICa(i) = params.g * (m(i-1).^params.p) * (V - ECa(i)); % Calcium current depends on the resting potential which changes with time per Nernst eqn.
            end
        end

        function [I] = ion_current(this, V, ion)
            % Ionic current for ion channel f is given by the equation:
            % If = gf * m^p * h^q * (V-Ef)
            % where, gf: maximum conductance for that ion channel
            % m: activation variable
            % p: activation exponent
            % h: inactivation variable
            % q: inactivation exponent
            % Ef: depolarization potential

            params = this.SetParams(ion, V); % Constants for Na, K, Ka, L. Resting potential (E) and conductance (for KCa) vary as a function of intracellular Ca2+ concentration and hence are a function of time for 'Ca' and 'KCa'
            [m, h] = this.get_gating_variable(params);

            if ion == 'Ca'
                [I, ~, ~] = this.get_calcium_current(V);
            elseif ion == 'KCa'
                [~, ECa, CCa_int] = this.get_calcium_current(V);
                CCa_diss = 1e-3; % Dissociation constant of Ca2+
                k = ((CCa_int./CCa_diss).^2);
                gCa = params.g * (k) ./ (1 + k); % Ligand-gated conductance of KCa varies with intracellular calcium concentration, hence is a function of time and calcium current
                I = gCa .* (m.^params.p) .* (h.^params.q) .* (V - ECa); % Resting potential is calculated by ECa which follows the Nernst eqn.
            else
                I = params.g * (m.^params.p) .* (h.^params.q) * (V - params.E);
            end
        end

        function [I] = get_ionic_current(this, V)
            ions = ["Na", "K", "Ka", "Ca", "KCa", "L"];
            time = [0: this.dt: 1];
            I = zeros(1,length(time));
            figure()
            for i = 1: length(ions)
                I_ion = ion_current(this, V, ions(i));
                I = I + I_ion;
                subplot(3,2,i)
                plot(time, I_ion);
                title(ions(i));
                grid on;
            end
        end
    end
end