from eppy.modeleditor import IDF
import os
import copy
import shutil
from eppy.runner.run_functions import runIDFs
import pandas as pd
import time





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


def Create_Schedule_Seasonal_Shading_Ctrl(idf1,epw_path,output_dir,fname1):
    #this function creates a seasonal shading availabilty schedule based on the last 24h rolling average of temperature of the epw file
    name_idf_shd_sched=fname1.split('.')[0]+'with_shd_sched'
    fname_shading_schedule=name_idf_shd_sched+'.idf'
    idf1.idfname = name_idf_shd_sched
    
    #read epw
    metadata, extracted_data = extract_epw_metadata_and_data(epw_path,False)
    WMO_ID=metadata['WMO_no'] #extract WMO ID and weather data format
    weather_format=metadata['Format']
    #__________________________________________________________________________________
    #rolling average last 24 hrs
    avg_window=24
    rolling_avg = extracted_data['Dry Bulb Temperature (°C)'].rolling(window=avg_window, min_periods=1).mean()
    #check if above
    threshold_sun=12
    Y_rolling_avg=rolling_avg>threshold_sun
    Y_rolling_avg = Y_rolling_avg.astype(int)
    Y_rolling_avg_24 = [Y_rolling_avg[i:i+24] for i in range(0, len(Y_rolling_avg), 24)]
    
    #Daily
    week_days=['Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday']
    special_days=['Holiday', 'SummerDesignDay', 'WinterDesignDay', 'CustomDay1','CustomDay2']
    w=1
    
    #add annual schedule
    #use GlobalCO2
    sched_annual_vorlage=copy.deepcopy(idf1.idfobjects['Schedule:Year'][0])
    idf1.idfobjects['Schedule:Year'].append(sched_annual_vorlage)
    sched_annual=idf1.idfobjects['Schedule:Year'][-1]
    shading_sched_name='Y_shading_yearly'
    sched_annual.Name=shading_sched_name
    sched_annual.Schedule_Type_Limits_Name='On/Off'
    
    for n in range (0,365):
        sched=idf1.newidfobject('Schedule:Day:Hourly')
        sched_name='Y_sun_d_'+str(n+1)
        sched.Name=sched_name
        sched.Schedule_Type_Limits_Name='On/Off'
        for i,item in enumerate(Y_rolling_avg_24[n],start=1):
            setattr(sched, f'Hour_{i}', item)
            
        #caclulate day of the week
        dow=(n%7)+1
        if dow == 1:
            sched_week=idf1.newidfobject('Schedule:Week:Daily')
            week_sched_name='Y_sun_w_'+str(w)
            sched_week.Name=week_sched_name
            #set special days (needed for syntax, i dont think is actually used!, just set to sunday)
            for sp_day in special_days:
                setattr(sched_week, f'{sp_day}_ScheduleDay_Name', sched_name)
            
            #yearly schedule rename
            setattr(sched_annual, f'ScheduleWeek_Name_{w}', week_sched_name)
            w=w+1
            
        dow_eng=week_days[dow-1]
        #create weekly schedule
        setattr(sched_week, f'{dow_eng}_ScheduleDay_Name', sched_name)
    
    #set remaining days in last week of year (only for syntax)
    for item in week_days[dow:]:
        setattr(sched_week, f'{item}_ScheduleDay_Name', sched_name)
        
    #idf1.saveas(os.path.join(output_dir,fname_shading_schedule)) #save as base
    
    return shading_sched_name, idf1, WMO_ID, weather_format



def gen_idf(epw_dir,epw_filename,template_dir,template_name,output_dir,WWR,window_orientation,strategy):
    #__________________________________________________________________________________
    #Template Selection

    #Select proper Template Office or Residential (Source for Gains SIA Merkblatt 2024)
    if template_name == "residential":
        fname1 = "Residential_totalBaseline.idf"
        fname1_path=os.path.join(template_dir,fname1)
        #drop file extension
        basecase_pre=fname1.split(".")[0]
    elif template_name == "office":
        fname1 = "Office_totalBaseline.idf"
        fname1_path=os.path.join(template_dir,fname1)
        basecase_pre=fname1.split(".")[0]

    iddfile="C:\EnergyPlusV22-2-0\Energy+.idd"
    IDF.setiddname(iddfile)
    epw_path=os.path.join(epw_dir,epw_filename)
    idf1 = IDF(fname1_path, epw_path)

    #__________________________________________________________________________________
    #Create shading schedule (last 24h rolling average above 12deg C, shading on)
    shading_sched_name, idf1, WMO_ID, weather_format = Create_Schedule_Seasonal_Shading_Ctrl(idf1,epw_path,output_dir,fname1)
    #__________________________________________________________________________________
    

    #Change the Orientation of the building
    building = idf1.idfobjects['BUILDING'][0]
    
    if window_orientation == "South":
        building.North_Axis = 0.0 #in template window is on south facade and 0.0 is the default
    elif window_orientation == "East":
        building.North_Axis = 270.0
    elif window_orientation == "West":
        building.North_Axis = 90.0
    elif window_orientation == "North":
        building.North_Axis = 180.0

    #__________________________________________________________________________________
    #Set WWR Ratio
    window_idf1=idf1.idfobjects['FenestrationSurface:Detailed'][0]
    upper_z_wind=window_idf1.Vertex_3_Zcoordinate
    height_wall=4 #4m
    new_z=upper_z_wind-height_wall*(WWR/100)
    window_idf1.Vertex_1_Zcoordinate=new_z
    window_idf1.Vertex_2_Zcoordinate=new_z
    
    #Ventilation Day Start
    #__________________________________________________________________________________
    
    idf_ventday=copy.deepcopy(idf1) #make a deep copy
    ventilation=idf_ventday.idfobjects['ZoneVentilation:WindandStackOpenArea'][0]
    ventilation.Opening_Area_Fraction_Schedule_Name="NatVent_day" #only avail between 7 am and 6pm
    ventilation.Opening_Area=1 #1m2
    ventilation.Height_Difference=1 #1m
    ventilation.Minimum_Indoor_Temperature=22 #22deg C

    ventday_pre=template_name+"_DayVent_"+str(WWR)+"_"+window_orientation+"_"+str(WMO_ID)+"_"+ weather_format
    fname_ventday=ventday_pre+".idf"
    idf_ventday.idfname = fname_ventday
    if strategy == "DayVent":
        idf_ventday.saveas(os.path.join(output_dir,fname_ventday))
        output_filename = fname_ventday

    #Shading start
    #__________________________________________________________________________________
    idf_shading=copy.deepcopy(idf_ventday) #make a deep copy
    #Add Shading
    shading = idf_shading.idfobjects['WindowShadingControl'][0]
    shading.Shading_Control_Type='OnIfHighSolarOnWindow'
    shading.Setpoint = 130 #setpoint for shading control
    #add seasonal schedule
    shading.Shading_Control_Is_Scheduled='Yes'
    shading.Schedule_Name=shading_sched_name

    shading_pre=template_name+"_Shading_"+str(WWR)+"_"+window_orientation+"_"+str(WMO_ID)+"_"+ weather_format
    fname_shading=shading_pre+".idf"
    idf_shading.idfname = fname_shading

    if strategy == "Shading":
        idf_shading.saveas(os.path.join(output_dir,fname_shading))
        output_filename = fname_shading
    #__________________________________________________________________________________

    #Ventilation Day/Night Start
    #__________________________________________________________________________________
    idf_ventall=copy.deepcopy(idf_shading) #make a deep copy

    ventilation=idf_ventall.idfobjects['ZoneVentilation:WindandStackOpenArea'][0]
    ventilation.Opening_Area_Fraction_Schedule_Name="AllOn" #shading_sched_name #"AllOn" #always available (at night too)


    ventall_pre=template_name+"_AllVent_"+str(WWR)+"_"+window_orientation+"_"+str(WMO_ID)+"_"+ weather_format
    fname_ventall=ventall_pre+".idf"
    idf_ventall.idfname = fname_ventall
    if strategy == "AllVent":
        idf_ventall.saveas(os.path.join(output_dir,fname_ventall))
        output_filename=fname_ventall
  

    #Thermal Mass Start
    #__________________________________________________________________________________
    idf_thermalmass=copy.deepcopy(idf_ventall) #make a deep copy
    detailed_surfaces= idf_thermalmass.idfobjects['BuildingSurface:Detailed']

    #set wall material
    walls=[sf for sf in detailed_surfaces if sf.Surface_Type=='Wall']
    for wall in walls:
        wall.Construction_Name = "UVal_0.2_Mass" # make sure such a construction exists in the model

    roofs=[sf for sf in detailed_surfaces if sf.Surface_Type=='Roof']
    for roof in roofs:
        roof.Construction_Name = "UVal_0.2_Mass" # make sure such a construction exists in the model

    floors=[sf for sf in detailed_surfaces if sf.Surface_Type=='Floor']
    for floor in floors:
        floor.Construction_Name = "Slab_Mass_FLIPPED" # make sure such a construction exists in the model

    ventall_mass_pre=template_name+"_AllVentMass_"+str(WWR)+"_"+window_orientation+"_"+str(WMO_ID)+"_"+ weather_format
    fname_thmass=ventall_mass_pre+".idf"
    idf_thermalmass.idfname = fname_thmass
    if strategy == "AllVentMass": #only save if this strategy is selected
        idf_thermalmass.saveas(os.path.join(output_dir,fname_thmass))
        output_filename = fname_thmass
    
    return output_filename, epw_path


start_time = time.time()

"""
#Single Test below:
#__________________________________________________________________________________
epw_dir="IWEC_test"
epw_filename="ARE_Abu.Dhabi.412170_IWEC.epw"
template_dir="Templates"
template_name="office"
output_dir="IDF_dir_test"
WWR=20
window_orientation="South"
strategy="AllVentMass"

output_filename=gen_idf(epw_dir,epw_filename,template_dir,template_name,output_dir,WWR,window_orientation,strategy)
print(output_filename)
#output_filename is: template_name +"_" + strategy + "_" + WWR + "_" + window_orientation + "_" + WMO_ID + "_" + weather_format + ".idf"
"""


#Iterations
WWRs = [20,30,40,50,60,70,80] #20,30,40,50,
orientations=["South","East","West","North"]
building_types = ["office","residential"]
strategies=["DayVent","Shading","AllVent","AllVentMass"]

epw_dir="IWEC_test"
template_dir="Templates"
output_dir="IDF_dir_test"

Simulations=[] #list of tuples (idf,epw)

for epw_filename in os.listdir(epw_dir):
    for template_name in building_types:
        for WWR in WWRs:
            for window_orientation in orientations:
                for strategy in strategies:
                    output_filename, corresponding_epw=gen_idf(epw_dir,epw_filename,template_dir,template_name,output_dir,WWR,window_orientation,strategy)
                    Simulations.append((output_filename,corresponding_epw)) #generate simulation list


print("Process finished --- %s seconds ---" % (time.time() - start_time))