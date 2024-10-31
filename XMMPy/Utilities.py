import os
import shutil
import io
from contextlib import redirect_stdout
from pysas.wrapper import Wrapper as w
from datetime import datetime

def format_observation_id(obs_id):
    """Add "0"s to beginning if length of Obs. ID < 10 (First 0 signifies the telescope is XMM Newton)

    Args:
        obs_id (str): The Obs. ID to be form
    Returns:
        str: Formatted Obs. ID
    """

    if len(str(obs_id)) < 10:
        num_of_missing = 10 - len(str(obs_id))
        return "0" * num_of_missing + str(obs_id)
    else:
        return str(obs_id)

def copy_over_files(work_dir, obs_dir):
    """Copies dataset, calibration index, and image files to the destination folder

    Args:
        work_dir (str): Absolute path of the folder to move the files from
        obs_dir (str): Absolute path of folder to move the files to
    """
    
    obs_id = obs_dir.split("/")[-1]

    #Source Files
    pn_old = os.path.join(work_dir, f"{obs_id}_PN.ds")
    mos1_old = os.path.join(work_dir, f"{obs_id}_MOS1.ds")
    mos2_old = os.path.join(work_dir, f"{obs_id}_MOS2.ds")
    pn_image_old = os.path.join(work_dir, f"{obs_id}_PN_Image.fits")
    mos1_image_old = os.path.join(work_dir, f"{obs_id}_MOS1_Image.fits")
    mos2_image_old = os.path.join(work_dir, f"{obs_id}_MOS2_Image.fits")
    ccf_old = os.path.join(work_dir, "ccf.cif")

    #Destination Files
    pn_new = os.path.join(obs_dir, f"{obs_id}_PN.ds")
    mos1_new = os.path.join(obs_dir, f"{obs_id}_MOS1.ds")
    mos2_new = os.path.join(obs_dir, f"{obs_id}_MOS2.ds")
    pn_image_new = os.path.join(obs_dir, f"{obs_id}_PN_Image.fits")
    mos1_image_new = os.path.join(obs_dir, f"{obs_id}_MOS1_Image.fits")
    mos2_image_new = os.path.join(obs_dir, f"{obs_id}_MOS2_Image.fits")
    ccf_new = os.path.join(obs_dir, "ccf.cif")

    #Copy them over
    shutil.copy(pn_old, pn_new)
    shutil.copy(mos1_old, mos1_new)
    shutil.copy(mos2_old, mos2_new)
    shutil.copy(pn_image_old, pn_image_new)
    shutil.copy(mos1_image_old, mos1_image_new)
    shutil.copy(mos2_image_old, mos2_image_new)
    shutil.copy(ccf_old, ccf_new)

def time_conv(input_time):
    """Converts mission time from MRT to Calendar format

    Args:
        input_time (str): Time to be converted from MRT

    Returns:
        output: Formatted date and time
    """

    output_capture = io.StringIO()
    inargs=[f"time={input_time}", "format=MRT"]
    with redirect_stdout(output_capture):
        w("xmmtimeconv", inargs).run()

    captured_output = output_capture.getvalue()
    selected_line = [line for line in captured_output.splitlines() if line.startswith("FITS")]
    string = (selected_line[0]).split(": ")[1]
    date = (datetime.strptime(string.split("T")[0], "%Y-%m-%d")).strftime("%B %d, %Y")
    time = (datetime.strptime(string.split("T")[1], "%H:%M:%S.%f")).strftime("%I:%M:%S %p")
    output = f"{date}; {time}"

    return output