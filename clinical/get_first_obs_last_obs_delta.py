import pandas as pd
data=pd.read_csv("APPiPSC_clinicalMRIdata_031522.csv",header=0,sep='\t')
cols_to_keep=["subj_id",
              "dm_status",
              "visit",
              "ados_severity",
              "ados_sa_severity",
              "ados_rrb_severity",
              "gi_symptoms",
              "gi_severity",
              "vitals_hc_cm",
              "vitals_wt_lb",
              "vitals_bmi",
              "mori_total_l1_cerebellum",
              "mori_total_l1_hemisphere_r",
              "mori_total_l1_hemisphere_l",
              "mori_total_l1_total_volume",
              "mori_total_l1_csf",
              "mori_total_l1_brainstem",
              "mori_total_l3_frontal_l",
              "mori_total_l3_frontal_r",
              "mori_total_l3_parietal_l",
              "mori_total_l3_parietal_r",
              "mori_total_l3_temporal_l",
              "mori_total_l3_temporal_r",
              "mori_total_l3_limbic_l",
              "mori_total_l3_limbic_r",
              "mori_total_l3_occipital_l",
              "mori_total_l3_occipital_r",
              "mori_total_l3_lateralventricle_l",
              "mori_total_l3_lateralventricle_r",
              "mori_total_l3_iii_ventricle",
              "mori_total_l3_iv_ventricle",
              "iq_viq",
              "iq_nviq",
              "iq_fsiq",
              "vine_communication_ss",
              "vine_living_ss",
              "vine_social_ss",
              "vine_motor_ss",
              "vine_abc_ss"]
data=data[cols_to_keep]
colnames=cols_to_keep[3::]
subject_to_condition={}
all_vals={} 
first_vals={}
last_vals={} 
delta_vals={}

for index,row in data.iterrows():
    subject=row['subj_id']
    condition=row['dm_status']
    if subject not in subject_to_condition:
        subject_to_condition[subject]=condition
    if subject not in all_vals:
        all_vals[subject]={}
    for key in colnames:
        cur_value=row[key]
        if str(cur_value) != "nan":
            if key not in all_vals[subject]:
                all_vals[subject][key]=[cur_value]
            else:
                all_vals[subject][key].append(cur_value)

for subject in all_vals:
    first_vals[subject]={}
    last_vals[subject]={}
    delta_vals[subject]={}
    
    for key in colnames:
        if key in all_vals[subject]:
            first=all_vals[subject][key][0]
            last=all_vals[subject][key][-1]
            first_vals[subject][key]=first
            last_vals[subject][key]=last 
            num_samples=len(all_vals[subject][key])
            if num_samples > 1:
                delta=last-first
                delta_vals[subject][key]=delta

            
first_val_df=pd.DataFrame.from_dict(first_vals,orient='index')
first_val_df.to_csv("first_vals.txt",header=True,index=True,sep='\t')

last_val_df=pd.DataFrame.from_dict(last_vals,orient='index')
last_val_df.to_csv("last_vals.txt",header=True,index=True,sep='\t')

delta_val_df=pd.DataFrame.from_dict(delta_vals,orient='index')
delta_val_df.to_csv("delta_vals.txt",header=True,index=True,sep='\t')

