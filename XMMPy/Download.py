from pysas.wrapper import Wrapper as w
import os
import shutil
import re
import glob
import warnings
import pandas as pd
from astropy.units import UnitsWarning
from astroquery.esa.xmm_newton import XMMNewton
import Utilities as utils
import numpy as np

def retrieve_obsids(source, working_dir):
    """Creates summary CSV file of all Obs. IDs for provided source

    Args:
        source (str): Valid format of source; could be catalog name, J2000 name, IAU name, etc.
        working_dir (str): Absolute path to folder where source folder will be created. Obs. ID file will be stored in source folder
    """

    work_dir = os.path.join(working_dir, source)
    os.makedirs(work_dir, exist_ok = True)
    os.chdir(work_dir)
    global obs_ids
    obs_ids = os.path.join(work_dir, f"{source}.csv")

    #Retrieve the Obs. IDs and save them to a file
    if not os.path.exists(obs_ids):
        dfs = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category = UnitsWarning)
            epic_source, cat_4xmm, stack_4xmm, slew_source = XMMNewton.get_epic_metadata(target_name=source)
            tables = [epic_source, cat_4xmm, stack_4xmm]
            table1 = [epic_source, cat_4xmm]
            for table in tables:
                if table in table1:
                    column_data = table["observation_id"]
                    data_values = column_data.data
                    df = pd.DataFrame(data_values, columns = ["Observation ID"])
                    df["Observation ID"] = df["Observation ID"].astype(str)
                    df["Observation ID"] = df["Observation ID"].apply(utils.format_observation_id)
                    df = df.drop_duplicates().reset_index(drop = True)
                    dfs.append(df)
                else:
                    column_data = table["obs_id"]
                    data_values = column_data.data
                    df = pd.DataFrame(data_values, columns = ["Observation ID"])
                    df["Observation ID"] = df["Observation ID"].astype(str)
                    df["Observation ID"] = df["Observation ID"].apply(utils.format_observation_id)
                    df = df.drop_duplicates().reset_index(drop = True)
                    dfs.append(df)
                            
            #Combine all unique Obs. IDs and remove duplicates and duds
            final_df = pd.concat(dfs).drop_duplicates().reset_index(drop = True)
            final_df = final_df[~final_df["Observation ID"].isna() & (final_df["Observation ID"] != 0)]
            final_df["Observation ID"] = df["Observation ID"].astype(str)
            final_df["Observation ID"] = final_df["Observation ID"].apply(utils.format_observation_id)
            final_df = final_df[~final_df["Observation ID"].isin(["0000000000", f"0000000{np.nan}"])]
            final_df.to_csv(obs_ids, index = False)

def download_and_reprocess_obsid(obs_id, work_dir):
    """Downloads and reprocesses data for provided Obs. ID

    Args:
        obs_id (str): Obs. ID to be downloaded and processed
        work_dir (str): Absolute path to folder where all Obs. ID data will be stored. Another folder with name of Obs. ID will be made 
        within it, and all data will be stored here
    """

    #Starting SAS Session
    inargs = [f"odfid={obs_id}",f"workdir={work_dir}"]
    w("startsas", inargs).run()

    #Moving files
    obs_dir = os.path.join(work_dir, obs_id)
    os.makedirs(obs_dir, exist_ok = True)
    shutil.move(os.path.join(work_dir, "ccf.cif"), obs_dir)
    for file in glob.glob(os.path.join(work_dir, "*SUM.SAS*")):
        shutil.move(file, obs_dir)

    #Add Environment Variables
    os.environ["SAS_CCF"] = os.path.join(obs_dir, "ccf.cif")
    matching_files = glob.glob(os.path.join(obs_dir, "*SUM.SAS"))
    for file in matching_files:
        odf = os.path.basename(file)
    os.environ["SAS_ODF"] = os.path.join(obs_dir, odf)
    os.chdir(obs_dir)

    #Reprocessing files
    w("epproc", []).run()
    w("emproc", []).run()

    #Renaming files
    p1 = r".*EPN.*ImagingEvts.*"
    p2 = r".*EMOS1.*ImagingEvts.*"
    p3 = r".*EMOS2.*ImagingEvts.*"

    for file in os.listdir(obs_dir):
        if re.search(p1, file):
            pn = os.path.join(obs_dir, file) 
        elif re.search(p2, file):
            mos1 = os.path.join(obs_dir, file)
        elif re.search(p3, file):  
            mos2 = os.path.join(obs_dir, file)

    os.rename(pn, f"{obs_id}_PN.ds")
    os.rename(mos1, f"{obs_id}_MOS1.ds")
    os.rename(mos2, f"{obs_id}_MOS2.ds")

    pn = f"{obs_id}_PN.ds"
    mos1  = f"{obs_id}_MOS1.ds"
    mos2 = f"{obs_id}_MOS2.ds"

    #Barycentric Corrections
    inargs = [f"table={pn}:EVENTS"]
    w("barycen", inargs).run()

    inargs = [f"table={mos1}:EVENTS"]
    w("barycen", inargs).run()

    inargs = [f"table={mos2}:EVENTS"]
    w("barycen", inargs).run()

    #Image Extraction
    xbin = 20
    ybin = 20
    xcoord = "X"
    ycoord = "Y"

    out_image_pn = os.path.join(obs_dir, f"{obs_id}_PN_Image.fits")
    inargs = [f"table={f"{obs_id}_PN.ds"}", f"imageset={out_image_pn}","withimageset=yes",f"xcolumn={xcoord}",f"ycolumn={ycoord}",f"ximagebinsize={xbin}",f"yimagebinsize={ybin}"]
    w("evselect", inargs).run()

    out_image_mos1 = os.path.join(obs_dir, f"{obs_id}_MOS1_Image.fits")
    inargs = [f"table={f"{obs_id}_MOS1.ds"}", f"imageset={out_image_mos1}","withimageset=yes",f"xcolumn={xcoord}",f"ycolumn={ycoord}",f"ximagebinsize={xbin}",f"yimagebinsize={ybin}"]
    w("evselect", inargs).run()

    out_image_mos2 = os.path.join(obs_dir, f"{obs_id}_MOS2_Image.fits")
    inargs = [f"table={f"{obs_id}_MOS2.ds"}", f"imageset={out_image_mos2}","withimageset=yes",f"xcolumn={xcoord}",f"ycolumn={ycoord}",f"ximagebinsize={xbin}",f"yimagebinsize={ybin}"]
    w("evselect", inargs).run()