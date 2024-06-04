using Plots
using ModelingToolkit
using DifferentialEquations

@connector function GasFlowPort(; name)
    # Parameters:
    # - `name`: Name of the connector
    states = @variables begin
        pressure(t) = 1.013e5    # Pressure in Pascals
        flow_rate(t) = 1.0       # Flow rate in g/s
        temperature(t) = 300.0   # Temperature in Kelvin
    end
    ODESystem(Equation[], t, states, []; name=name)
end

function ElectricalGasPort(; name)
    # Parameters:
    # - `name`: Name of the connector
    @named positive_pin = Pin()
    @named negative_pin = Pin()
    @named gas_port = GasFlowPort()
    states = @variables begin
        voltage(t) = 1.0          # Voltage in volts
        current(t) = 1.0          # Current in amperes
        temperature(t) = 573      # Temperature in Kelvin
        flow_rate(t) = 0          # Flow rate in g/s
    end
    parameters = @parameters begin
        pressure = 1.013e5       # Pressure in Pascals
    end
    equations = [
        voltage ~ positive_pin.voltage - negative_pin.voltage,
        0 ~ positive_pin.current + negative_pin.current,
        temperature ~ gas_port.temperature,
        flow_rate ~ gas_port.flow_rate,
        current ~ positive_pin.current
    ]
    return compose(ODESystem(equations, t, states, []; name=Symbol(name)), positive_pin, negative_pin, gas_port)
end

function HydrogenTank(; name, initial_pressure=4e6)
    # Parameters:
    # - `name`: Name of the tank
    # - `initial_pressure`: Initial pressure in Pascals
    @named gas_in = GasFlowPort()
    @named gas_out = GasFlowPort()
    states = @variables begin
        pressure(t) = initial_pressure # Pressure in Pascals
        temperature(t) = 293           # Temperature in Kelvin
        gas_density(t) = initial_pressure / 293 / 4124 # Density of gas in kg/m3
        metal_density(t) = 8400        # Density of metal hydride in kg/m3
        absorption_rate(t) = gas_in.flow_rate # Absorption rate in kg/m3/s
        equilibrium_pressure(t) = 0    # Equilibrium pressure in Pascals
    end
    parameters = @parameters(
        Cp_metal = 419,        # Specific heat of metal hydride in J/kg/K
        Cp_gas = 14890,        # Specific heat of hydrogen gas in J/kg/K
        reaction_heat = 35620, # Heat of reaction in J/mol
        initial_metal_density = 8400,  # Density of the metal without hydrogen in kg/m3
        saturated_metal_density = 8517, # Density of the MH at steady state in kg/m3
        A = 10.570,            # Constant for equilibrium pressure calculation
        B = 3704.6,            # Constant for equilibrium pressure calculation
        pressure_ref = 100000, # Reference pressure in Pascals
        void_ratio = 0.5,      # Void ratio of the hydrogen storage tank
        absorption_constant = 59.187,  # Constant for hydrogen absorption rate in 1/s
        activation_energy = 21170,     # Activation energy in J/mol
        R = 8.314,            # Ideal gas constant in J/mol/K
        molar_mass_H2 = 2.016e-3, # Molar mass of H2 in kg/mol
        gas_constant = 8.314 / 2.016e-3 # Specific gas constant for hydrogen in J/kg/K
    )
    equations = [
        equilibrium_pressure ~ pressure_ref * exp(A - B / temperature),
        absorption_rate ~ absorption_constant * exp(-activation_energy / temperature / R) * log(pressure / equilibrium_pressure) * (saturated_metal_density - metal_density), # Absorption rate in kg/m3/s
        ∂(gas_density) ~ -absorption_rate / void_ratio,
        ∂(metal_density) ~ absorption_rate / (1 - void_ratio),
        ∂(temperature) ~ (absorption_rate * reaction_heat / gas_constant + absorption_rate * temperature * (Cp_gas - Cp_metal)) / (void_ratio * gas_density * Cp_gas + (1 - void_ratio) * metal_density * Cp_metal),
        pressure ~ gas_density * gas_constant * temperature
    ]
    return compose(ODESystem(equations, t, states, parameters; name=Symbol(name)), gas_in, gas_out)
end

function HydrogenDesorptionTank(; name, initial_pressure=101325, reference_temperature=313.5, reference_flow_rate=10, max_void_ratio=0.7, initial_temperature=300, saturation_degree=0.99, reference_heat_transfer=300)
    # Parameters:
    # - `name`: Name of the tank
    # - `initial_pressure`: Initial pressure in Pascals
    # - `reference_temperature`: Reference temperature in Kelvin
    # - `reference_flow_rate`: Reference flow rate in liters per minute
    # - `max_void_ratio`: Maximum void ratio
    # - `initial_temperature`: Initial temperature in Kelvin
    # - `saturation_degree`: Degree of saturation
    # - `reference_heat_transfer`: Reference heat transfer coefficient in W/m2
    @named gas_in = GasFlowPort()
    @named gas_out = GasFlowPort()
    states = @variables begin
        tank_pressure(t) = initial_pressure       # Pressure in Pascals
        tank_temperature(t) = initial_temperature # Temperature in Kelvin
        hydrogen_mass_rate(t) = 0                 # Absorption or desorption rate in kg/m3/s
        metal_density(t) = 8517 * saturation_degree # Density of the metal hydride in kg/m3
        gas_density(t) = initial_pressure / initial_temperature / 4124 # Density of gas in kg/m3
        gas_flow_out(t) = 0                       # Gas flow rate out of the tank in kg/m3/s
        equilibrium_pressure(t) = 0               # Equilibrium pressure in Pascals
        water_temperature_out(t) = 0              # Circulating water outflow temperature in Kelvin
        hydride_proportion(t) = 0                 # Proportion of hydride
        heat_release(t) = 0                       # Heat release of circulating water in Joules
    end
    parameters = @parameters(
        tank_diameter = 0.3,         # Tank diameter in meters
        tank_length = 0.5,           # Tank length in meters
        tank_volume = pi * 0.3^2 / 4 * 0.5, # Tank volume in cubic meters
        void_ratio_max = max_void_ratio, # Maximum void ratio
        metal_hydride_volume = pi * 0.3^2 / 4 * 0.5 * max_void_ratio, # Volume of metal hydride in cubic meters
        surface_area = pi * 0.3 * 0.5, # Surface area in square meters
        water_flow_rate = reference_flow_rate / 60, # Water flow rate in kg/s
        heat_transfer_coefficient = reference_heat_transfer, # Heat transfer coefficient in W/m2
        initial_metal_density = 8400,   # Density of the metal without hydrogen in kg/m3
        saturated_metal_density = 8517, # Density of the MH at steady state in kg/m3
        specific_heat_metal = 419,      # Specific heat of metal hydride in J/kg/K
        specific_heat_water = 1860,     # Specific heat of water in J/kg/K
        specific_heat_gas = 14890,      # Specific heat of hydrogen gas in J/kg/K
        void_ratio = 0.5,               # Void ratio of the hydrogen storage tank
        molar_mass_H2 = 2.016e-3,       # Molar mass of H2 in kg/mol
        gas_constant = 8.314,           # Ideal gas constant in J/mol/K
        specific_gas_constant = 8.314 / 2.016e-3, # Specific gas constant for hydrogen in J/kg/K
        reaction_heat = 1.54e7,         # Heat of reaction in J/kgH2
        activation_energy = 21179.6,    # Activation energy in J/mol
        desorption_constant = 59.187,   # Constant for hydrogen desorption rate in 1/s
        overall_heat_transfer_coefficient = reference_heat_transfer, # Overall heat transfer coefficient in W/m2
        alpha = reference_heat_transfer * pi * 0.3 * 0.5 / (reference_flow_rate / 60 * 1860), # Alpha coefficient for heat transfer
        alpha1 = 1,                    # Alpha1 coefficient for equilibrium pressure calculation
        alpha2 = 0.5,                  # Alpha2 coefficient for equilibrium pressure calculation
        beta = 0.137,                  # Beta coefficient for equilibrium pressure calculation
        constant_a = 13.44,            # Constant A for equilibrium pressure calculation
        constant_b = 3780,             # Constant B for equilibrium pressure calculation
        constant_phi = 0.038,          # Constant phi for equilibrium pressure calculation
        outlet_pressure = 101325,      # Outlet pressure in Pascals
        inlet_water_temperature = reference_temperature, # Inlet water temperature in Kelvin
        outlet_diameter = 0.005,       # Outlet diameter of the tank in meters
        throttle_area = pi * 0.005^2 / 4, # Throttle area at the tank exit in square meters
        specific_heat_ratio = 1.409,   # Specific heat ratio for hydrogen
        reference_pressure = 101325,   # Reference pressure in Pascals
        discharge_coefficient = 9.57,  # Discharge coefficient for the outlet
        energy_discharge = 16473       # Energy discharge constant in J/mol
    )
    equations = [
        tank_pressure ~ gas_density * specific_gas_constant * tank_temperature,
        hydride_proportion ~ (metal_density - initial_metal_density) / (saturated_metal_density - initial_metal_density),
        equilibrium_pressure ~ reference_pressure * exp(constant_a - constant_b / tank_temperature + constant_phi * tan(alpha1 * pi * (hydride_proportion - alpha2)) + beta / 2),
        hydrogen_mass_rate ~ -desorption_constant * exp(-activation_energy / specific_gas_constant / tank_temperature) * log(tank_pressure / equilibrium_pressure) * (saturated_metal_density - metal_density), # Mass rate of hydrogen desorption in kg/m3/s
        ∂(gas_density) ~ (hydrogen_mass_rate - gas_flow_out) / (void_ratio_max - 1 + void_ratio),
        ∂(metal_density) ~ -hydrogen_mass_rate / (1 - void_ratio),
        ∂(tank_temperature) ~ (-hydrogen_mass_rate * reaction_heat + heat_release / metal_hydride_volume) / ((void_ratio_max - 1 + void_ratio) * specific_heat_gas * gas_density + (1 - void_ratio) * specific_heat_metal * metal_density),
        water_temperature_out ~ tank_temperature + (inlet_water_temperature - tank_temperature) * exp(-alpha),
        heat_release ~ water_flow_rate * specific_heat_water * (inlet_water_temperature - tank_temperature) * (1 - exp(-alpha)),
        gas_flow_out ~ tank_pressure / metal_hydride_volume / sqrt(specific_gas_constant * tank_temperature) * throttle_area * sqrt(specific_heat_ratio) * 0.578, # Gas flow rate in kg/m3/s
        gas_out.temperature ~ tank_temperature,
        gas_out.flow_rate ~ gas_flow_out
    ]
    return compose(ODESystem(equations, t, states, parameters; name=Symbol(name)), gas_in, gas_out)
end

function PEMElectrolyzer(; name, initial_temperature=313.15, initial_pressure=1.023, number_of_cells=1)
    # Parameters:
    # - `name`: Name of the electrolyzer
    # - `initial_temperature`: Initial temperature in Kelvin
    # - `initial_pressure`: Initial pressure in atm
    # - `number_of_cells`: Number of electrolyzer cells in series
    @named electrical_gas_port = ElectricalGasPort()
    @unpack voltage, current = electrical_gas_port
    states = @variables begin
        hydrogen_molar_rate(t)  # Molar rate of hydrogen production in mol/s
        faraday_efficiency(t)   # Faraday efficiency, dimensionless
        hydrogen_mass_rate(t) = 0 # Mass rate of hydrogen production in g/s
        hydrogen_mass_rate_derivative(t) = 0 # Derivative of hydrogen mass rate
        voltage_ohmic(t) = 0    # Ohmic voltage drop
        voltage_activation(t) = 0 # Activation voltage drop
        heat_transfer_ohmic(t) = 0 # Ohmic heat transfer
    end
    parameters = @parameters(
        electrode_area = 20 * 0.0001, # Electrode area in m2
        reversible_voltage = 1.229,   # Reversible voltage in V
        electrolyzer_temperature = initial_temperature, # Electrolyzer operating temperature in K
        electrolyzer_pressure = initial_pressure, # Electrolyzer operating pressure in atm
        number_of_cells = number_of_cells, # Number of electrolyzer cells in series
        resistance_coefficient1 = 4.45153e-5, # Ohmic resistance parameter 1
        resistance_coefficient2 = 6.88874e-9, # Ohmic resistance parameter 2
        ohmic_coefficient1 = -3.12996e-6,     # Ohmic resistance parameter 1
        ohmic_coefficient2 = 4.47137e-7,      # Ohmic resistance parameter 2
        activation_coefficient = 0.33824,     # Activation overpotential coefficient
        activation_coefficient1 = -0.01539,   # Activation overpotential coefficient 1
        activation_coefficient2 = 2.00181,    # Activation overpotential coefficient 2
        activation_coefficient3 = 15.24178,   # Activation overpotential coefficient 3
        faraday_coefficient1 = 478645.74,     # Faraday efficiency coefficient 1
        faraday_coefficient2 = -953.15,       # Faraday efficiency coefficient 2
        specific_faraday_coefficient1 = 478645.74 + (-953.15 * initial_temperature), # Specific Faraday efficiency coefficient 1
        specific_faraday_coefficient2 = 1.03960 + (-0.00104 * initial_temperature), # Specific Faraday efficiency coefficient 2
        resonant_frequency = 17,              # Resonant frequency in Hz
        damping_ratio = 0.7,                  # Damping ratio
        faraday_constant = 96485,             # Faraday constant in C/mol
        number_of_electrons = 2,              # Number of electrons transferred
        hydrogen_molar_mass = 2.016,          # Molar mass of hydrogen in g/mol
        overpotential_coefficient1 = 0.09901, # Overpotential coefficient 1
        overpotential_coefficient2 = -0.00207,# Overpotential coefficient 2
        overpotential_coefficient3 = 1.31064e-5, # Overpotential coefficient 3
        overpotential_coefficient4 = -0.08483,# Overpotential coefficient 4
        overpotential_coefficient5 = 0.00179, # Overpotential coefficient 5
        overpotential_coefficient6 = -1.13390e-5, # Overpotential coefficient 6
        overpotential_coefficient7 = 1481.45, # Overpotential coefficient 7
        overpotential_coefficient8 = -23.60345,# Overpotential coefficient 8
        overpotential_coefficient9 = -0.25774, # Overpotential coefficient 9
        activation_energy_coefficient1 = 3.71417, # Activation energy coefficient 1
        activation_energy_coefficient2 = -0.93063,# Activation energy coefficient 2
        activation_energy_coefficient3 = 0.05817, # Activation energy coefficient 3
        activation_energy_coefficient4 = -3.72068,# Activation energy coefficient 4
        activation_energy_coefficient5 = 0.93219, # Activation energy coefficient 5
        activation_energy_coefficient6 = -0.05826,# Activation energy coefficient 6
        activation_energy_coefficient7 = -18.38215,# Activation energy coefficient 7
        activation_energy_coefficient8 = 5.87316, # Activation energy coefficient 8
        activation_energy_coefficient9 = -0.46425 # Activation energy coefficient 9
    )
    equations = [
        voltage ~ reversible_voltage + (resistance_coefficient1 + resistance_coefficient2 * (electrolyzer_temperature - 273.15) + ohmic_coefficient1 + ohmic_coefficient2 * electrolyzer_pressure) * (current / electrode_area) + activation_coefficient * log10((activation_coefficient1 + activation_coefficient2 / (electrolyzer_temperature - 273.15) + activation_coefficient3 / (electrolyzer_temperature - 273.15)^2) * abs(current / electrode_area) + 1),
        voltage_ohmic ~ (resistance_coefficient1 + resistance_coefficient2 * (electrolyzer_temperature - 273.15) + ohmic_coefficient1 + ohmic_coefficient2 * electrolyzer_pressure) * (current / electrode_area),
        voltage_activation ~ activation_coefficient * log10((activation_coefficient1 + activation_coefficient2 / (electrolyzer_temperature - 273.15) + activation_coefficient3 / (electrolyzer_temperature - 273.15)^2) * abs(current / electrode_area) + 1),
        faraday_efficiency ~ specific_faraday_coefficient2 * ((current / electrode_area)^2 / ((current / electrode_area)^2 + specific_faraday_coefficient1)),
        hydrogen_molar_rate ~ (current / electrode_area) * number_of_cells * faraday_efficiency / (number_of_electrons * faraday_constant),
        ∂(hydrogen_mass_rate) ~ hydrogen_mass_rate_derivative,
        ∂(hydrogen_mass_rate_derivative) ~ 2 * damping_ratio * resonant_frequency * hydrogen_mass_rate_derivative - resonant_frequency^2 * hydrogen_mass_rate + resonant_frequency^2 * hydrogen_molar_mass * hydrogen_molar_rate,
        heat_transfer_ohmic ~ overpotential_coefficient1 + overpotential_coefficient2 * electrolyzer_temperature + overpotential_coefficient3 * electrolyzer_temperature^2 + (overpotential_coefficient4 + overpotential_coefficient5 * electrolyzer_temperature + overpotential_coefficient6 * electrolyzer_temperature^2) * exp((overpotential_coefficient7 + overpotential_coefficient8 * electrolyzer_temperature + overpotential_coefficient9 * electrolyzer_temperature^2) / (current / electrode_area)) + activation_energy_coefficient1 + activation_energy_coefficient2 * electrolyzer_pressure + activation_energy_coefficient3 * electrolyzer_pressure^2 + (activation_energy_coefficient4 + activation_energy_coefficient5 * electrolyzer_pressure + activation_energy_coefficient6 * electrolyzer_pressure^2) * exp((activation_energy_coefficient7 + activation_energy_coefficient8 * electrolyzer_pressure + activation_energy_coefficient9 * electrolyzer_pressure^2) / (current / electrode_area))
    ]
    extend(ODESystem(equations, t, states, parameters; name=Symbol(name)), electrical_gas_port)
end

function PistonCompressor(; name, outlet_pressure=350000)
    # Parameters:
    # - `name`: Name of the compressor
    # - `outlet_pressure`: Outlet pressure in Pascals
    @named gas_in = GasFlowPort()
    @named gas_out = GasFlowPort()
    parameters = @parameters(
        polytropic_exponent = 1.35,            # Polytropic exponent for compression
        inlet_pressure = gas_in.pressure,      # Inlet pressure in Pascals
        outlet_pressure = outlet_pressure,     # Outlet pressure in Pascals
        inlet_temperature = 298.15,            # Inlet temperature in Kelvin
        outlet_temperature = 298.15 * (outlet_pressure / gas_in.pressure)^((1.35 - 1) / 1.35) # Outlet temperature in Kelvin
    )
    equations = [
        gas_out.temperature ~ outlet_temperature,
        gas_out.pressure ~ outlet_pressure,
        gas_out.flow_rate ~ gas_in.flow_rate
    ]
    return compose(ODESystem(equations, t, [], parameters; name=name), gas_in, gas_out)
end

function pv_battery(reference_temperature, reference_irradiance, reference_short_circuit_current, reference_open_circuit_voltage, reference_maximum_power_current, reference_maximum_power_voltage, actual_temperature, actual_irradiance)
    # Parameters:
    # - `reference_temperature`: Reference temperature in Kelvin
    # - `reference_irradiance`: Reference irradiance in W/m^2
    # - `reference_short_circuit_current`: Reference short circuit current in A
    # - `reference_open_circuit_voltage`: Reference open circuit voltage in V
    # - `reference_maximum_power_current`: Reference maximum power current in A
    # - `reference_maximum_power_voltage`: Reference maximum power voltage in V
    # - `actual_temperature`: Actual temperature in Kelvin
    # - `actual_irradiance`: Actual irradiance in W/m^2
    temperature_coefficient_current = 0.025
    temperature_coefficient_voltage = 0.0005
    irradiance_coefficient_voltage = 0.00288

    delta_temperature = actual_temperature - reference_temperature
    delta_irradiance = actual_irradiance - reference_irradiance

    short_circuit_current = reference_short_circuit_current * actual_irradiance / reference_irradiance * (1 + temperature_coefficient_current * delta_temperature)
    open_circuit_voltage = reference_open_circuit_voltage * (1 - irradiance_coefficient_voltage * delta_temperature) * log(exp(1) + temperature_coefficient_voltage * delta_irradiance)
    maximum_power_current = reference_maximum_power_current * actual_irradiance / reference_irradiance * (1 + temperature_coefficient_current * delta_temperature)
    maximum_power_voltage = reference_maximum_power_voltage * (1 - irradiance_coefficient_voltage * delta_temperature) * log(exp(1) + temperature_coefficient_voltage * delta_irradiance)

    constant_C2 = (maximum_power_voltage / open_circuit_voltage - 1) * (log(1 - maximum_power_current / short_circuit_current))^-1
    constant_C1 = (1 - maximum_power_current / short_circuit_current) * exp(-maximum_power_voltage / (constant_C2 * open_circuit_voltage))

    return short_circuit_current, open_circuit_voltage, maximum_power_current, maximum_power_voltage, constant_C1, constant_C2
end
