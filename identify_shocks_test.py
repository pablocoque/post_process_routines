import numpy as np

class shock_integration:
    def __init__(self, time_series, mass_evo, rh_evo, component):
        self.time = time_series
        self.mass_evo = mass_evo
        self.rh_evo = rh_evo
        self.component = component
        self.integrating_flag = np.zeros_like(time_series, dtype=bool)
        self.shock_events = 0
        self.tidal_heating_events = []
        self.tidal_heating_cumulative = 0.
        self.single_shock_integral = 0.
        self.single_shock_duration = 0.

    def __identify_shocks_method_1__(self):
        # Dummy implementation for shock identification
        # In a real scenario, this would involve complex logic
        max_value = 0.
        for i in range(1, len(self.time)):
            dt = self.time[i] - self.time[i-1]
            if self.integrating_flag[i-1]:
                if(self.component[i-1] * self.component[i] < 0):
                    self.single_shock_integral += 0.5 * 0.5 * dt * np.abs(self.component[i-1])
                    self.single_shock_duration += 0.5 * dt
                    self.__calculate_heating_parameter__(i)
                    self.__restart_single_shock__()
                elif((np.abs(self.component[i]) < 0.88 * np.abs(self.component[i-1])) or
                     (np.abs(self.component[i]) > np.abs(self.component[i-1]))):
                    self.__calculate_heating_parameter__(i)
                    self.__restart_single_shock__()
                else:
                    self.single_shock_integral += 0.5 * dt * np.abs(self.component[i] + self.component[i-1])
                    self.single_shock_duration += dt
                    self.integrating_flag[i] = True
            else:
                if(self.component[i] * self.component[i-1] < 0):
                    self.single_shock_integral += 0.5 * 0.5 * dt * np.abs(self.component[i])
                    self.single_shock_duration += 0.5 * dt
                    self.integrating_flag[i] = True
                    max_value = np.abs(self.component[i])
                    # self.__calculate_heating_parameter__(i)
                    # self.__restart_single_shock__()
                elif(np.abs(self.component[i]) < np.abs(self.component[i-1])):
                    self.single_shock_integral += 0.5 * dt * np.abs(self.component[i] + self.component[i-1])
                    self.single_shock_duration += dt
                    self.integrating_flag[i] = True
                    max_value = np.abs(self.component[i])
                else:
                    continue

    def __identify_shocks_method_2__(self):
        # Dummy implementation for shock identification
        # In a real scenario, this would involve complex logic
        max_value = 0.
        behavior = 0
        for i in range(1, len(self.time)):
            dt = self.time[i] - self.time[i-1]
            if self.integrating_flag[i-1]:
                if(self.component[i-1] * self.component[i] < 0):
                    self.single_shock_integral += 0.5 * 0.5 * dt * np.abs(self.component[i-1])
                    self.single_shock_duration += 0.5 * dt
                    self.__calculate_heating_parameter__(i)
                    self.__restart_single_shock__()
                    max_value = 0.
                    behavior = 0
                elif((np.abs(self.component[i]) < 0.88 * max_value) or
                     ((np.abs(self.component[i]) > np.abs(self.component[i-1]))
                       and (behavior < 0))):
                    self.__calculate_heating_parameter__(i)
                    self.__restart_single_shock__()
                    max_value = 0.
                    behavior = 0
                else:
                    self.single_shock_integral += 0.5 * dt * np.abs(self.component[i] + self.component[i-1])
                    self.single_shock_duration += dt
                    self.integrating_flag[i] = True
                    if(np.abs(self.component[i]) > max_value):
                        max_value = np.abs(self.component[i])
                    if(np.abs(self.component[i])<np.abs(self.component[i-1])):
                        behavior = -1
            else:
                if(self.component[i] * self.component[i-1] < 0):
                    self.single_shock_integral += 0.5 * 0.5 * dt * np.abs(self.component[i])
                    self.single_shock_duration += 0.5 * dt
                    max_value = np.abs(self.component[i])
                    self.integrating_flag[i] = True
                    behavior = 1
                    # self.__calculate_heating_parameter__(i)
                    # self.__restart_single_shock__()
                elif(np.abs( 0.88 * self.component[i]) > np.abs(self.component[i-1])):
                    self.single_shock_integral += 0.5 * dt * np.abs(self.component[i] + self.component[i-1])
                    self.single_shock_duration += dt
                    max_value = np.abs(self.component[i])
                    self.integrating_flag[i] = True
                    behavior = 1
                

    def __calculate_heating_parameter__(self, index):
        Itid = self.single_shock_integral**2 * self.calculate_adiabatic_correction(index)
        self.tidal_heating_events.append([self.time[index], Itid])
        self.tidal_heating_cumulative += Itid
    
    def calculate_adiabatic_correction(self, index):
        G = 4.3e-3 * (3.15e16/3.09e13)**2 # Gravitational constant pc3 msun-1 gyr-2
        tdyn = 8 * np.pi * self.rh_evo[index]**3 / (3 * self.mass_evo[index] * G)
        return (1 + self.single_shock_duration**2 / tdyn**2)**(-1.5)
        
    def __restart_single_shock__(self):
        self.shock_events += 1
        self.single_shock_integral = 0.
        self.single_shock_duration = 0.

    def restart_all_shocks(self):
        self.shock_events = 0
        self.tidal_heating_events = []
        self.tidal_heating_cumulative = 0.
        self.integrating_flag = np.zeros_like(self.time, dtype=bool)

    def get_shock_events_method_1(self):
        self.__identify_shocks_method_1__()

    def get_shock_events_method_2(self):
        self.__identify_shocks_method_2__()