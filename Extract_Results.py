import numpy as np
import pandas as pd
from pythermalcomfort import pmv
import psychrolib
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import os

def CreateCustomPoly(coordinates, n, color, p, mode, linestyle, label):
    
    if mode == 'rel_hum':
        # New list to store interpolated coordinates
        interpolated_coordinates = []
        
        # New list to store interpolated coordinates
        interpolated_coordinates = []
        # Interpolate points with the same y-coordinate
        for i in range(len(coordinates) - 1):
            x1, y1 = coordinates[i]
            x2, y2 = coordinates[i + 1]
            
            interpolated_coordinates.append(coordinates[i])
            
            if y1 == y2:
                interpolated_x = np.linspace(x1, x2, n + 2)[1:-1]  # Exclude first and last points
                for x in interpolated_x:
                    interpolated_coordinates.append((x, y1))
                    
        
        # Add the last point
        interpolated_coordinates.append(coordinates[-1])
        
        #Create Polygons
        coordinates_t_humrat=[]
        for i in interpolated_coordinates:
            coordinates_t_humrat.append((i[0],psychrolib.GetHumRatioFromRelHum(i[0], i[1], p)))
        
    elif mode == 'hum_rat':
        coordinates_t_humrat=coordinates
        
    elif mode =='vap_pres':
        coordinates_t_humrat=[]
        for i in coordinates:
            pw=i[1]*133.322 #mmHg to Pa
            #https://www.built-envi.com/wp-content/uploads/cae331_513_lecture12_psychrometric-equations.pdf
            coordinates_t_humrat.append((i[0],0.622*(pw/(p-pw))))
        
    poly1=Polygon(coordinates_t_humrat)
    
    #plot

    x,y = poly1.exterior.xy
    #plot, = plt.plot(x,y,color=color, linestyle=linestyle, label=label)
    plot="N/A"
    
    return plot, poly1



def F_to_C(coordinates):
    converted_coordinates = []
    for x, y in coordinates:
        x_celsius = (x - 32) * 5 / 9
        converted_coordinates.append((x_celsius, y))
    return converted_coordinates

def extract_epw_metadata_and_data(file_path,header_only):
    """
    Extract metadata (including latitude and longitude) and specific weather data from an EPW file.

    Parameters:
    file_path (str): Path to the EPW file.

    Returns:
    tuple: A tuple containing a dictionary with metadata and a DataFrame with the extracted weather data.
    """

    # Read the first 8 lines for metadata
    with open(file_path, 'r') as file:
        header_lines = [next(file) for _ in range(8)]

    # Extracting latitude and longitude from the metadata
    location_data = header_lines[0].split(',')
    metadata = {
        'Location': location_data[1].strip(),
        'Format': location_data[4].strip().split(' ')[0],
        'WMO_no': str(location_data[5]),
        'Latitude': float(location_data[6]),
        'Longitude': float(location_data[7]),
        'Time Zone': float(location_data[8]),
        'Elevation': float(location_data[9])
    }
    
    extracted_data=[]
    if not(header_only):
        # EPW files have a fixed column structure, these are the columns of interest for weather data
        columns_of_interest = {
            'Dry Bulb Temperature': 6,  # in degrees Celsius
            'Relative Humidity': 8,  # in percent
            'Global Horizontal Radiation': 13,  # Wh/m2
            'Wind Speed': 21,  # m/s
        }
    
        # Read the EPW file for weather data, skipping the header
        data = pd.read_csv(file_path, skiprows=8, header=None, encoding='utf-8')
    
        # Extract only the columns of interest
        extracted_data = data[[columns_of_interest['Dry Bulb Temperature'],
                               columns_of_interest['Relative Humidity'],
                               columns_of_interest['Global Horizontal Radiation'],
                               columns_of_interest['Wind Speed']]]
    
        # Rename the columns for clarity
        extracted_data.columns = ['Dry Bulb Temperature (°C)', 'Relative Humidity (%)',
                                  'Global Horizontal Radiation (Wh/m2)', 'Wind Speed (m/s)']

    return metadata, extracted_data

def indoor_wind(u_met):
    #Ashrae fundamentals; airflow round buildings 
    #Urban/suburban areas
    U_indoor=u_met*(270/10)**(0.14)*(10/370)**(0.22)
    return U_indoor

def CliCo_Comfort_Model(Y_airspeed,epw_folder,df):
    strategy_str=df.index.get_level_values('strategy').unique() #strategy
    strategy=strategy_str[0]
    epw_file_str=df.index.get_level_values('epw_filename').unique()
    epw=epw_file_str[0]

    metadata, extracted_data = extract_epw_metadata_and_data(os.path.join(epw_folder,epw),False)
    extracted_data['Wind_indoor']=np.vectorize(indoor_wind)(extracted_data['Wind Speed (m/s)'])
    
    if strategy in ["AllVent","AllVentMass"]: 
        ## Altitude in meters to pressure in Pascals
        ## Reference  https://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html 
        ### Altitude of location in Meters 
        AltitudeFromSeaLevel = 10
        p = (101325)*(1-(2.25577*(10)**(-5))*(AltitudeFromSeaLevel))**5.25588
        #______________________________________________________________________________
        #Start draw comfort polygons - Climate Consultant
        #______________________________________________________________________________
        
        n = 1  # Number of points to interpolate
        
        T_comf_low_left=21
        T_comf_up_left=19
        T_comf_low_right=27.5
        T_comf_up_right=26.5
        
        coordinates_ClimateConsulatant_all=[(T_comf_up_left,12/1000),(T_comf_low_left,0),(T_comf_low_right,0),(T_comf_up_right,12/1000)]
        poly_ClimateConsulatant_all,poly1=CreateCustomPoly(coordinates_ClimateConsulatant_all, n, "green", p, "hum_rat",'-','Comfort')
        
        #for ventilation max valocity is 1.5 m/s (3.7 °C Reduction)
        Red=3.7 #°C
        coordinates_Ventilation=[(T_comf_low_right,0),(T_comf_up_right,12/1000),(T_comf_up_right+Red,12/1000),(T_comf_low_right+Red,0)]
        poly_Ventilation,poly2=CreateCustomPoly(coordinates_Ventilation, n, "cyan", p, "hum_rat",'-','Ventilation')
        
        #Heating
        Red=0 #°C
        coordinates_heating=[(18.7-Red,13.5/1000),(T_comf_low_left-Red,0),(-50,0),(-50,13.5/1000)]
        poly_heating,poly5=CreateCustomPoly(coordinates_heating, n, "red", p, "hum_rat",'-','Heating')
        
        #for fan, max velocity is 0.8 m/s (3 °C Reduction)
        Red=2.7 #°C
        coordinates_Fan=[(T_comf_low_right,0),(T_comf_up_right,12/1000),(T_comf_up_right+Red,12/1000),(T_comf_low_right+Red,0)]
        poly_Fan,poly6=CreateCustomPoly(coordinates_Fan, n, "lime", p, "hum_rat",'-','Fan')
        
        
        #______________________________________________________________________________
        #End draw comfort polygons - Climate Consultant
        #______________________________________________________________________________

        #______________________________________________________________________________
        #START Polygon checking
        Simulation_results = pd.DataFrame(columns=['Heating', 'Comfort','Ventilation_poly','Fan','Comfort_with_Fan_NatVent'])
        #Generate the points for the polygon
        points={}
        points['Sim']=np.vectorize(Point)(df["Zone Mean Air Temperature"], df['Zone Mean Air Humidity Ratio'])
        #check heating
        Simulation_results['Heating'] =np.vectorize(poly5.contains)(points['Sim'])
        #check comfort
        Simulation_results['Comfort'] =np.vectorize(poly1.contains)(points['Sim'])

        #__________________________________________________________________________
        #NatVent
        #__________________________________________________________________________
        #Nat_Vent calculate wind speed

        #check if inside the ventilation polygon
        inside_vent_poly = np.vectorize(poly2.contains)(points['Sim'])

        if Y_airspeed=="ASHRAE_fundamentals":
            air_speed=extracted_data['Wind_indoor']

        if Y_airspeed=="Sim_half_window":
            rho_air=1.204
            area=0.5
            air_speed_sim=df['Zone Ventilation Mass Flow Rate'] / rho_air / area
            air_speed = air_speed_sim.reset_index(drop=True)
        
    
        #is airspeed enough
        wind_effective = (air_speed> 0.2) #& (extracted_data['Wind_indoor'] <= 1.6)

        window_open_vent= (df['Zone Ventilation Mass Flow Rate'] > 0)
        window_open_vent_reset = window_open_vent.reset_index(drop=True)
        Simulation_results['Ventilation_poly']=inside_vent_poly*wind_effective*window_open_vent_reset 


        #__________________________________________________________________________

        #Fan
        Simulation_results['Fan'] = np.vectorize(poly6.contains)(points['Sim'])
        
        #Fan and Vent
        columns_to_check = ['Fan', 'Comfort', 'Heating', 'Ventilation_poly'] #'thermal_mass_night_vent'
        Simulation_results['Comfort_with_Fan_NatVent'] = Simulation_results[columns_to_check].any(axis=1).astype(int)
        
        columns_to_check = ['Comfort', 'Heating', 'Ventilation_poly'] #'thermal_mass_night_vent'
        Simulation_results['Comfort_with_NatVent'] = Simulation_results[columns_to_check].any(axis=1).astype(int)
        
        columns_to_check = ['Comfort', 'Heating', 'Fan'] #'thermal_mass_night_vent'
        Simulation_results['Comfort_with_Fan'] = Simulation_results[columns_to_check].any(axis=1).astype(int)
        
        columns_to_check = ['Comfort', 'Heating'] #'thermal_mass_night_vent'
        Simulation_results['Comfort_Baseline'] = Simulation_results[columns_to_check].any(axis=1).astype(int)
        
        #occupancy check above for everything
        disc_h_Sim_CliCoPoly_vent=8760-(Simulation_results['Comfort_with_NatVent'].sum())
        disc_h_Sim_CliCoPoly_fan=8760-(Simulation_results['Comfort_with_Fan'].sum())
        disc_h_Sim_CliCoPoly_fan_vent=8760-(Simulation_results['Comfort_with_Fan_NatVent'].sum())
        disc_h_Baseline=8760-(Simulation_results['Comfort_Baseline'].sum())
        #how many hours in comfort + heating (as baseline?)

        occ_profile_str=df['People Occupant Count']
        occ_profile = occ_profile_str.reset_index(drop=True)
        occ_bin=occ_profile>0

        disc_h_Sim_CliCoPoly_vent_occ_only=((1-Simulation_results['Comfort_with_NatVent'])*occ_bin).sum()
        disc_h_Sim_CliCoPoly_fan_occ_only=((1-Simulation_results['Comfort_with_Fan'])*occ_bin).sum()
        disc_h_Sim_CliCoPoly_fan_vent_occ_only=((1-Simulation_results['Comfort_with_Fan_NatVent'])*occ_bin).sum()
        disc_h_Baseline_occ_only=((1-Simulation_results['Comfort_Baseline'])*occ_bin).sum()

        
        
    else:
        disc_h_Sim_CliCoPoly_vent="N/A"
        disc_h_Sim_CliCoPoly_fan="N/A"
        disc_h_Sim_CliCoPoly_fan_vent="N/A"
        disc_h_Baseline="N/A"

        disc_h_Sim_CliCoPoly_vent_occ_only="N/A"
        disc_h_Sim_CliCoPoly_fan_occ_only="N/A"
        disc_h_Sim_CliCoPoly_fan_vent_occ_only="N/A"
        disc_h_Baseline_occ_only="N/A"

    #make data frame for the results
    data = {'Disc_h_Sim_CliCoPoly_vent': [disc_h_Sim_CliCoPoly_vent],
                            'Disc_h_Sim_CliCoPoly_fan': [disc_h_Sim_CliCoPoly_fan],
                            'Disc_h_Sim_CliCoPoly_fan_vent': [disc_h_Sim_CliCoPoly_fan_vent],
                            'Disc_h_Baseline': [disc_h_Baseline],
                            'Disc_h_Sim_CliCoPoly_vent_occ_only': [disc_h_Sim_CliCoPoly_vent_occ_only],
                            'Disc_h_Sim_CliCoPoly_fan_occ_only': [disc_h_Sim_CliCoPoly_fan_occ_only],
                            'Disc_h_Sim_CliCoPoly_fan_vent_occ_only': [disc_h_Sim_CliCoPoly_fan_vent_occ_only],
                            'Disc_h_Baseline_occ_only': [disc_h_Baseline_occ_only]}
    
    if Y_airspeed=="ASHRAE_fundamentals":
        index=['airspeed_from_wind_ashrae']
    if Y_airspeed=="Sim_half_window":
        index=['airspeed_from_sim']
    results = pd.DataFrame(data, index=index)    
    
    return results, extracted_data['Wind_indoor']

def calculate_PMV_related_metrics(PMV_df,index,Tair_indoor,Mrt_indoor,RH_indoor,vr, met, clo, PMV_too_hot_thres,hum_Rat_indoor,hum_Rat_thres,occ):
    #below is vectorized
    hourly_PMV=np.vectorize(pmv)(tdb=Tair_indoor, 
                                       tr=Mrt_indoor, 
                                       vr=vr, 
                                       rh=RH_indoor,
                                       met=met, 
                                       clo=clo, 
                                       units='SI',
                                       limit_inputs=False,
                                       standard='ASHRAE')   #ASHRAE uses SET model for elevated airspeed
    hourly_PMV_too_hot = hourly_PMV>PMV_too_hot_thres
    
    hourly_PMV_comf=(hourly_PMV > -0.5) & (hourly_PMV < PMV_too_hot_thres)
    #impose a hard 0.012 kgw/kgda upper limit for humidity
    hum_rat_too_high=hum_Rat_indoor>hum_Rat_thres
    uncomf_bc_hum_rat_too_high=hum_rat_too_high*hourly_PMV_comf
    
    
    uncomf_hours_with_humlim=(hourly_PMV_too_hot+uncomf_bc_hum_rat_too_high)!=0
    
    #PMV_hours above 0.5
    PMV_hours_annual=(hourly_PMV-PMV_too_hot_thres)*hourly_PMV_too_hot
    #overheating PMV hours
    
    #Sums
    #sums that do not have the upper hum limit imposed
    sum_overheating=hourly_PMV_too_hot.sum()
    sum_PMV_hours=(PMV_hours_annual).sum()
    #sums that have the upper humlimit imposed
    sum_overheating_humlim=(uncomf_hours_with_humlim).sum()
    
    #only when occupied sums
    occ_binary=occ>0
    sum_overheating_occ=(hourly_PMV_too_hot*occ_binary).sum()
    sum_PMV_hours_occ=(PMV_hours_annual*occ_binary).sum()
    #sums that have the upper humlimit imposed
    sum_overheating_humlim_occ=(uncomf_hours_with_humlim*occ_binary).sum()
    
    new_row_df = pd.DataFrame({'Sum Overheating': [sum_overheating], 
                               'Sum PMV Hours': [sum_PMV_hours], 
                               'Sum Overheating humlim': [sum_overheating_humlim],
                               'occ Sum Overheating': [sum_overheating_occ], 
                               'occ Sum PMV Hours': [sum_PMV_hours_occ], 
                               'occ Sum Overheating humlim': [sum_overheating_humlim_occ]}, 
                              index=[index])
    
    PMV_df = pd.concat([PMV_df, new_row_df])

    return PMV_df



def Simulation_Post_Proc(df,epw_folder):
    #Climate Consulatant comofrt model start
    #________________________________________________________________________________ 
    Y_airspeed="ASHRAE_fundamentals" #Either  or "Sim_half_window"
    psychrolib.SetUnitSystem(psychrolib.SI) #important!!
    results1, airspeed_from_Ashrae_wind = CliCo_Comfort_Model(Y_airspeed,epw_folder,df)
    Y_airspeed="Sim_half_window"
    results2, airspeed_from_Ashrae_wind = CliCo_Comfort_Model(Y_airspeed,epw_folder,df)    
    results_CliCo_comf_model_sim=pd.concat([results1, results2], axis=0)

    #PMV_START
    #________________________________________________________________________________    
    # Create an empty DataFrame with specified column names
    columns = ['Sum Overheating', 'Sum PMV Hours', 'Sum Overheating humlim','occ Sum Overheating', 'occ Sum PMV Hours', 'occ Sum Overheating humlim']
    PMV_df = pd.DataFrame(columns=columns)
    
    #still air
    PMV_too_hot_thres=0.5
    clo=0.5
    met=1.1
    still_air=0.1
    hum_Rat_thres=0.012
    PMV_df=calculate_PMV_related_metrics(PMV_df,"still_air",df["Zone Mean Air Temperature"],df['Zone Mean Radiant Temperature'],df['Zone Air Relative Humidity'],still_air,met, clo, PMV_too_hot_thres,df['Zone Mean Air Humidity Ratio'],hum_Rat_thres,df['People Occupant Count'])
    
    #nat vent air speed
    rho_air=1.204 #kg/m3
    area=0.5 #m2 half the window area
    #mass flow is in kg/s
    velocity_vent=(df['Zone Ventilation Mass Flow Rate'] / rho_air / area).clip(upper=1.6,lower=still_air) #set 1.5 m/s as upper limit; window only
    #ONLY TAKE OVER 0.2 m/s into account?? from Climate Consultant help
    velocity_vent_newmin = velocity_vent.mask(velocity_vent < 0.2, still_air)
    PMV_df=calculate_PMV_related_metrics(PMV_df,"velocity_from_window",df["Zone Mean Air Temperature"],df['Zone Mean Radiant Temperature'],df['Zone Air Relative Humidity'],velocity_vent_newmin,met, clo, PMV_too_hot_thres,df['Zone Mean Air Humidity Ratio'],hum_Rat_thres,df['People Occupant Count'])

    #fan and nat vent air speed
    vr_fan=0.8
    velocity_vent_fan=velocity_vent.clip(lower=vr_fan)#ensure always at least 0.8; fan and window
    PMV_df=calculate_PMV_related_metrics(PMV_df,"velocity_from_window_and_fan",df["Zone Mean Air Temperature"],df['Zone Mean Radiant Temperature'],df['Zone Air Relative Humidity'],velocity_vent_fan,met, clo, PMV_too_hot_thres,df['Zone Mean Air Humidity Ratio'],hum_Rat_thres,df['People Occupant Count'])
    
    #fan only
    PMV_df=calculate_PMV_related_metrics(PMV_df,"fan_only",df["Zone Mean Air Temperature"],df['Zone Mean Radiant Temperature'],df['Zone Air Relative Humidity'],vr_fan,met, clo, PMV_too_hot_thres,df['Zone Mean Air Humidity Ratio'],hum_Rat_thres,df['People Occupant Count'])
    
    #NatVent_with "indoor wind" from ASHRAE
    SET_with_Ashrae_wind=(airspeed_from_Ashrae_wind).clip(upper=1.6,lower=still_air)
    SET_with_Ashrae_wind_newmin = SET_with_Ashrae_wind.mask(SET_with_Ashrae_wind < 0.2, still_air)
    
    window_open=df['Zone Ventilation Mass Flow Rate']
    window_open = window_open.reset_index(drop=True)
    
    SET_with_Ashrae_wind_newmin=SET_with_Ashrae_wind_newmin*(window_open > 0)
    PMV_df=calculate_PMV_related_metrics(PMV_df,"window_with_Ashrae_fund_wind",df["Zone Mean Air Temperature"],df['Zone Mean Radiant Temperature'],df['Zone Air Relative Humidity'],SET_with_Ashrae_wind_newmin,met, clo, PMV_too_hot_thres,df['Zone Mean Air Humidity Ratio'],hum_Rat_thres,df['People Occupant Count'])

    #PMV_END
    #________________________________________________________________________________
    
    #HEATING_START
    #________________________________________________________________________________
    #Heating
    area=59.55
    heating_kWh=df['Zone Ideal Loads Zone Total Heating Energy']/3.6e6
    heating_kWh_tol=(heating_kWh/3.6e6)>0.1 #0.1 kWh is the tolerance for counting as a heating hour
    heating_hours=heating_kWh_tol.sum()
    heating_demand=(heating_kWh/area).sum()
    
    #HEATING_END
    #________________________________________________________________________________
    #shading hours
    shading_hours=df['Surface Shading Device Is On Time Fraction'].sum()
    
    Result_dict={}
    Result_dict["PMV"]=PMV_df
    Result_dict["CliCo_comf_model_sim"]=results_CliCo_comf_model_sim
    Result_dict["Heating_hours"]=heating_hours
    Result_dict["Heating_Demand_kWH_per_m2"]=heating_demand
    Result_dict["Shading_down_hours"]=shading_hours
    
    return Result_dict

df = pd.read_hdf("ARE_Abu.Dhabi.412170_IWEC_office_20_East_AllVent.hdf")
df = df.stack(level=-1) # convert wide form to long form dataframe
df = df.droplevel([0,1],axis=1) # drop the unnecessary energyplus categorization of different meters)



epw_folder="IWEC"
Result_dict=Simulation_Post_Proc(df,epw_folder)







