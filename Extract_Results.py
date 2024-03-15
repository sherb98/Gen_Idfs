import numpy as np
import pandas as pd
from pythermalcomfort import pmv



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

def Simulation_Post_Proc(df):

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
    velocity_vent=(df['Zone Ventilation Mass Flow Rate'] / rho_air / area).clip(upper=1.5,lower=still_air) #set 1.5 m/s as upper limit; window only
    PMV_df=calculate_PMV_related_metrics(PMV_df,"velocity_from_window",df["Zone Mean Air Temperature"],df['Zone Mean Radiant Temperature'],df['Zone Air Relative Humidity'],velocity_vent,met, clo, PMV_too_hot_thres,df['Zone Mean Air Humidity Ratio'],hum_Rat_thres,df['People Occupant Count'])
    
    #fan and nat vent air speed
    vr_fan=0.8
    velocity_vent_fan=velocity_vent.clip(lower=vr_fan)#ensure always at least 0.8; fan and window
    PMV_df=calculate_PMV_related_metrics(PMV_df,"velocity_from_window_and_fan",df["Zone Mean Air Temperature"],df['Zone Mean Radiant Temperature'],df['Zone Air Relative Humidity'],velocity_vent_fan,met, clo, PMV_too_hot_thres,df['Zone Mean Air Humidity Ratio'],hum_Rat_thres,df['People Occupant Count'])
    
    
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
    Result_dict["Heating_hours"]=heating_hours
    Result_dict["Heating_Demand_kWH_per_m2"]=heating_demand
    Result_dict["Shading_down_hours"]=shading_hours
    
    return Result_dict

df = pd.read_hdf("ARE_Abu.Dhabi.412170_IWEC_office_20_East_AllVent.hdf")
df = df.stack(level=-1) # convert wide form to long form dataframe
df = df.droplevel([0,1],axis=1) # drop the unnecessary energyplus categorization of different meters)

Result_dict=Simulation_Post_Proc(df)







