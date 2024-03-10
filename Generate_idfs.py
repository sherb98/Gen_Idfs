#from eppy import modeleditor
from eppy.modeleditor import IDF
import os
import copy
import shutil
from eppy.runner.run_functions import runIDFs
import pandas as pd
import time


start_time = time.time()


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


def Create_Schedule_Seasonal_Shading_Ctrl(idf1,epwfile,Sim_dir,fname1):
    #this function creates a seasonal shading availabilty schedule based on the last 24h rolling average of temperature of the epw file
    name_idf_shd_sched=fname1.split('.')[0]+'with_shd_sched'
    fname_shading_schedule=name_idf_shd_sched+'.idf'
    idf1.idfname = name_idf_shd_sched
    
    #read epw
    metadata, extracted_data = extract_epw_metadata_and_data(epwfile,False)
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
        
    idf1.saveas(os.path.join(Sim_dir,fname_shading_schedule)) #save as base
    
    return shading_sched_name, fname_shading_schedule, idf1



def gen_idf_run_sim(Sim_dir,building_type,epwfile,WWRs,name_epw_current,orientations):
    #Generates the idf files needed per epw file
    #__________________________________________________________________________________
    #Set Simulation Directory

    weather_Folder="Weather"

    #__________________________________________________________________________________
    #Template Selection

    #Select proper Template Office or Residential (Source for Gains SIA Merkblatt 2024)
    if building_type == "residential":
        fname1 = "Residential_totalBaseline.idf"
        #drop file extension
        basecase_pre=fname1.split(".")[0]
    elif building_type == "office":
        fname1 = "Office_totalBaseline.idf"
        basecase_pre=fname1.split(".")[0]

    iddfile="C:\EnergyPlusV22-2-0\Energy+.idd"
    IDF.setiddname(iddfile)
    idf1 = IDF(fname1, epwfile)

    #__________________________________________________________________________________
    #Create shading schedule (last 24h rolling average above 12deg C, shading on)
    shading_sched_name, fname_shading_schedule, idf1 = Create_Schedule_Seasonal_Shading_Ctrl(idf1,epwfile,Sim_dir,fname1)
    #__________________________________________________________________________________
    
    csv_prefixes=[]
    idf_list=[]

    for window_orientation in orientations:

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

        for n, WWR in enumerate(WWRs):
            
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
            ventilation.Opening_Area_Fraction_Schedule_Name="NatVentAvail_Zone_0" #only avail between 7 am and 6pm
            ventilation.Opening_Area=1 #1m2
            ventilation.Height_Difference=1 #1m
            ventilation.Minimum_Indoor_Temperature=22 #22deg C
        
            ventday_pre=building_type+"_DayVent_"+str(WWR)+"_"+window_orientation+"_"+name_epw_current
            fname_ventday=ventday_pre+".idf"
            idf_ventday.idfname = fname_ventday
            idf_ventday.saveas(os.path.join(Sim_dir,fname_ventday))
        
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
        
            shading_pre=building_type+"_Shading_"+str(WWR)+"_"+window_orientation+"_"+name_epw_current
            fname_shading=shading_pre+".idf"
            idf_shading.idfname = fname_shading
            idf_shading.saveas(os.path.join(Sim_dir,fname_shading))
            #__________________________________________________________________________________
        
            #Ventilation Day/Night Start
            #__________________________________________________________________________________
            idf_ventall=copy.deepcopy(idf_shading) #make a deep copy
        
            ventilation=idf_ventall.idfobjects['ZoneVentilation:WindandStackOpenArea'][0]
            ventilation.Opening_Area_Fraction_Schedule_Name="AllOn" #shading_sched_name #"AllOn" #always available (at night too)

        
            ventall_pre=building_type+"_AllVent_"+str(WWR)+"_"+window_orientation+"_"+name_epw_current
            fname_ventall=ventall_pre+".idf"
            idf_ventall.idfname = fname_ventall
            idf_ventall.saveas(os.path.join(Sim_dir,fname_ventall))
        
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
        
            ventall_mass_pre=building_type+"_AllVent_mass_"+str(WWR)+"_"+window_orientation+"_"+name_epw_current
            fname_thmass=ventall_mass_pre+".idf"
            idf_thermalmass.idfname = fname_thmass
            idf_thermalmass.saveas(os.path.join(Sim_dir,fname_thmass))
        
        
        
            idf_shading.idfname = fname_shading
            idf_ventday.idfname = fname_ventday
            idf_ventall.idfname = fname_ventall
            idf_thermalmass.idfname = fname_thmass
            
            csv_prefixes.append([ventday_pre, shading_pre, ventall_pre, ventall_mass_pre])
            idf_list.append([idf_shading.idfname,idf_ventday.idfname,idf_ventall.idfname,idf_thermalmass.idfname])

    return csv_prefixes, idf_list, iddfile




if __name__ == '__main__':

    #Metadata 1 Start
    #__________________________________________________________________________________
    building_types = ["office","residential"]
    WWRs = [20,30,40,50,60,70,80] #20,30,40,50,
    orientations=["South","East","West","North"]
    #this is only FYI; it is hard-coded in the function
    #Variants=["_DayVent_", "_Shading_", "_AllVent_", "_AllVent_mass_"]
    #__________________________________________________________________________________
    #Metadata 1 End


    epw_folder="IWEC_test"

    #epws=[]
    epws_corr_ent=[]
    names_epws=[]
    
    csv_prefixes=[]
    idf_list=[]

    Sim_dir="IDF_dir"#,building_type) #create a directory for the simulations
    if os.path.exists(Sim_dir):
        shutil.rmtree(Sim_dir)
    os.makedirs(Sim_dir)

    for file in os.listdir(epw_folder):
        epwfile=os.path.abspath(os.path.join(epw_folder,file))
        data,dummy= extract_epw_metadata_and_data(epwfile,True)

        #Metadata 2 start
        #__________________________________________________________________________________
        name_epw_current=str(data['WMO_no'])+"_"+ data['Format'] #extract WMO ID and weather data format
        #__________________________________________________________________________________
        #Metadata 2 end

        names_epws.append(name_epw_current)

        for building_type in building_types:
            #__________________________________________________________________________________
            csv_prefixes_cur,idf_list_cur,iddfile=gen_idf_run_sim(Sim_dir,building_type,epwfile,WWRs,name_epw_current,orientations)
            #idfs are named like so: building_type_variant_WWR_orientation_epwname.idf ; epw name: WMOnumber_format
            #__________________________________________________________________________________

            idf_list_cur_flat=[item for sublist in idf_list_cur for item in sublist] #flatten list
            csv_list_cur_flat=[item+'.csv' for sublist in csv_prefixes_cur for item in sublist] #flatten list

            for iter in range (0,len(csv_list_cur_flat)): #append epw file the same number of times as the number of simulations
                epws_corr_ent.append(epwfile) #this is the corresponding epw list for the simulations

            csv_prefixes.append(csv_list_cur_flat)
            idf_list.append(idf_list_cur_flat)

    idf_list_flat=[item for sublist in idf_list for item in sublist] #flatten list
    csv_list_flat=[item for sublist in csv_prefixes for item in sublist] #flatten list

    main_dir=os.path.dirname(os.path.realpath(__file__))

    print(csv_list_flat)
    print(idf_list_flat)

    print("Process finished --- %s seconds ---" % (time.time() - start_time))
    








