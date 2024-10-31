from pysas.wrapper import Wrapper as w
import os
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import Plotting as plot
import Utilities as utils
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import warnings
from astropy.utils.exceptions import AstropyWarning
from XMMPy import values 

def save_source_region(obs_dir, source_name, instr, radius = 0.0055):
    """Save Region file for Source with Given Coordinates

    Args:
        obs_dir (str): Absolute path to working folder where region file is to be generated. Folder name should be of the Obs. ID
        source_name (str): Valid format of source; could be catalog name, J2000 name, IAU name, etc.
        instr (str): Chooses instrument to save regions for, with options being "pn", "mos1", "mos2" for PN, MOS1, and MOS2 respectively
        radius (float, optional): Radius of annulus to be saved in degrees. Defaults to 0.004.
    """

    warnings.filterwarnings("ignore", category = AstropyWarning)
    obs_id = obs_dir.split("/")[-1]

    #Open image and get header info
    fits_file = os.path.join(obs_dir, f"{obs_id}_{instr.upper()}_Image.fits")
    with fits.open(fits_file) as hdul:
        header = hdul[0].header
        wcs = WCS(header)

    coords = source_name.split("J")[1]

    #Extract coordinates
    if "+" in coords:
        sign_pos = coords.index("+")
    elif "-" in coords:
        sign_pos = coords.index("-")

    ra_raw = coords[:sign_pos]  
    dec_raw = coords[sign_pos:] 
    ra = f"{ra_raw[:2]}:{ra_raw[2:4]}:{ra_raw[4:]}"
    dec = f"{dec_raw[:3]}:{dec_raw[3:5]}:{dec_raw[5:]}"

    #Convert RA and Dec from sexagesimal to decimal degrees
    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    ra_deg = coord.ra.deg
    dec_deg = coord.dec.deg

    #Convert Decimal Degree WCS to Pixel coordinates
    pixel_coords = wcs.all_world2pix(ra_deg, dec_deg, 1)
    x_pix, y_pix = pixel_coords

    #Calculate the pixel scale in both directions
    scale_x = wcs.wcs.cdelt[0] 
    scale_y = wcs.wcs.cdelt[1]
    average_scale = (abs(scale_x) + abs(scale_y)) / 2 

    #Convert the radius from degrees to pixels using the average scale
    radius_pix = radius / average_scale

    #Apply the LTM/LTV transformations 
    ltm = header.get("LTM1_1", 1)
    ltm2 = header.get("LTM2_2", 1)
    ltv1 = header.get("LTV1", 0)   
    ltv2 = header.get("LTV2", 0)  

    #Calculate final values
    x_phys = (x_pix - ltv1) / ltm if ltm else x_pix
    y_phys = (y_pix - ltv2) / ltm2 if ltm2 else y_pix
    physical_radius = radius_pix / ltm if ltm else radius_pix

    #Output strings
    region = f"circle({x_phys:.3f}, {y_phys:.3f}, {physical_radius:.3f}) # color=white"
    region_content = '''# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
physical
'''

    #Make region file
    with open(os.path.join(obs_dir, f"{instr.lower()}_regions.reg"), "w") as reg_file:
        reg_file.write(region_content.strip() + "\n")
        reg_file.write(region + "\n")

def epic_lightcurve(obs_dir, energy_min, energy_max, binsize, band, instr, bkg_sub):
    """EPIC light curve generation function

    Args:
        obs_dir (str): Absolute path to working folder where region file is to be generated. Folder name should be of the Obs. ID
        energy_min (str): Minimum energy of band in eV
        energy_max (str): Maximum energy of band in eV
        binsize (int): Binsize (s)
        band (str): Name of energy band
        instr (str): Chooses instrument to save regions for, with options being "pn", "mos1", "mos2" for PN, MOS1, and MOS2 respectively
        bkg_sub (str): Choose to background subtract, "yes" or "no"
    """

    obs_id = obs_dir.split("/")[-1]
    os.chdir(obs_dir)

    #Light curve Parameters
    if instr.lower() == "pn":
        q_flag = "#XMMEA_EP" 
        pattern = 4   
        epic = f"{obs_id}_PN.ds"
    else:
        q_flag = "#XMMEA_EM" 
        pattern = 12
        if instr.lower() == "mos1":
            epic = f"{obs_id}_MOS1.ds"
        else:
            epic = f"{obs_id}_MOS2.ds"

    if bkg_sub == "yes":
        with open(os.path.join(obs_dir, f"{instr.lower()}_regions.reg"), "r") as file:
            lines = file.readlines()
        region1 = lines[3].strip()
        region2 = lines[4].strip()
        c1 = region1.partition("color")[2].replace("=", "").replace("\n", "")

        if(c1.split(" ")[0] == "white"):
            src_region = region1.split(" # ")[0]
            bkg_region = region2.split(" # ")[0]
        else:
            src_region = region2.split(" # ")[0]
            bkg_region = region1.split(" # ")[0]

        src_lc = f"EPIC_{instr.upper()}_Source_{band}_Lightcurve.lc"
        expression = f"{q_flag}&&(PATTERN<={pattern})&& ((X,Y) IN {src_region})&&(PI in [{energy_min}:{energy_max}])"
        inargs = [f"table={epic}", "energycolumn=PI", "withrateset=yes", f"rateset={src_lc}", f"timebinsize={binsize}", "maketimecolumn=yes", "makeratecolumn=yes", f"expression={expression}"]
        w("evselect", inargs).run()

        bkg_lc = f"EPIC_{instr.upper()}_Background_{band}_Lightcurve.lc"
        expression = f"{q_flag}&&(PATTERN<={pattern})&& ((X,Y) IN {bkg_region})&&(PI in [{energy_min}:{energy_max}])"
        inargs = [f"table={epic}", "energycolumn=PI", "withrateset=yes", f"rateset={bkg_lc}", f"timebinsize={binsize}", "maketimecolumn=yes", "makeratecolumn=yes", f"expression={expression}"]
        w("evselect", inargs).run()

        sub_lc = f"EPIC_{instr.upper()}_Corrected_{band}_Lightcurve.lc"
        inargs = [f"eventlist={epic}", f"srctslist={src_lc}", f"outset={sub_lc}", f"bkgtslist={bkg_lc}", "withbkgset=yes", "applyabsolutecorrections=yes"]
        w("epiclccorr", inargs).run()   

        os.remove(src_lc)
        os.remove(bkg_lc)

    else:
        with open(os.path.join(obs_dir, f"{instr.lower()}_regions.reg"), "r") as file:
            lines = file.readlines()
        region1 = lines[3].strip()
        region = region1.split(" # ")[0]

        #Count Rate
        src_rate_lc = f"EPIC_{instr.upper()}_Source_Rate_{band}_Lightcurve.lc"
        expression = f"{q_flag}&&(PATTERN<={pattern})&& ((X,Y) IN {region})&&(PI in [{energy_min}:{energy_max}])"
        inargs = [f"table={epic}", "energycolumn=PI", "withrateset=yes", f"rateset={src_rate_lc}", f"timebinsize={binsize}", "maketimecolumn=yes", "makeratecolumn=yes", f"expression={expression}"]
        w("evselect", inargs).run()
            
        #Counts
        src_cts_lc = f"EPIC_{instr.upper()}_Source_Counts_{band}_Lightcurve.lc"
        expression = f"{q_flag}&&(PATTERN<={pattern})&& ((X,Y) IN {region})&&(PI in [{energy_min}:{energy_max}])"
        inargs = [f"table={epic}", "energycolumn=PI", "withrateset=yes", f"rateset={src_cts_lc}", f"timebinsize={binsize}", "maketimecolumn=yes", "makeratecolumn=no", f"expression={expression}"]
        w("evselect", inargs).run()

def generate_epic_lightcurve(obs_dir, binsize, source_name, instr, bkg_sub):
    """Create all the EPIC light curves for a source

    Args:
        obs_dir (str): Absolute path to working folder where region file is to be generated. Folder name should be of the Obs. ID
        binsize (int): Binsize (s)
        source_name (str): Valid format of source; could be catalog name, J2000 name, IAU name, etc.
        instr (str): Chooses instrument to save regions for, with options being "pn", "mos1", "mos2" for PN, MOS1, and MOS2 respectively
        bkg_sub (str): Choose to background subtract, "yes" or "no"
    """

    os.chdir(obs_dir)
    binsize = int(binsize)
    obs_id = obs_dir.split("/")[-1]
    
    if bkg_sub == "yes":
        for _, row in values.iterrows():
            band = row["Band"]
            energy_min = row["Energy Min"]
            energy_max = row["Energy Max"]

            epic_lightcurve(obs_dir, energy_min, energy_max, binsize, band, instr, bkg_sub = "yes")

        rate_dict = {}
        error_list = []
        time_list = []
        bin_list = []

        with fits.open(os.path.join(obs_dir, f"EPIC_{instr.upper()}_Corrected_Broadband_Lightcurve.lc")) as hdul:
            data = hdul[1].data
            for entry in data["ERROR"]:
                error_list.append(entry)

            for entry in data["TIME"]:
                time_list.append(entry)
                bin_list.append((entry - time_list[0]) / 500 + 1)

        for _, row in values.iterrows():
            band = row["Band"]
            rate_file = os.path.join(obs_dir, f"EPIC_{instr.upper()}_Corrected_{band}_Lightcurve.lc")

            with fits.open(rate_file) as hdul:
                data = hdul[1].data
                rate = []
                for value in data["RATE"]:
                    rate.append(value)
                rate_dict[band] = rate

            os.remove(rate_file)

        time = pd.Series(time_list)
        errors = pd.Series(error_list)
        rates = pd.DataFrame(rate_dict)
        rates.insert(0, "Time", time)
        rates.insert(1, "Bin", pd.Series(map(int, bin_list)))
        rates["Error"] = errors

        rates.to_csv(os.path.join(obs_dir, f"{source_name}_{obs_id}_{instr.upper()}_Count_Rates.csv"))
        rates_csv = os.path.join(obs_dir, f"{source_name}_{obs_id}_{instr.upper()}_Count_Rates.csv")

        if rates["Broadband"].sum() == 0:
            os.remove(rates_csv)
        else:
            #Plotting
            fig, axs = plt.subplots(3, 1, figsize = (10, 15), dpi = 300)

            plt.sca(axs[0])
            axs[0].set_title("Broadband Count Rate", fontsize = 12)
            axs[0].text(0.9935, 1.05, f"Start: {utils.time_conv(rates["Time"].iloc[0])}", transform = axs[0].transAxes, fontsize = 8, verticalalignment = "top", horizontalalignment = "right", bbox = dict(facecolor = "white", edgecolor = "black", boxstyle = "square,pad=0.5", linewidth = 0.8))
            plot.rate_plotter(plt, file = rates_csv, band = "Broadband", error = "Yes")
            plt.grid(True, linestyle = "--")

            plt.sca(axs[1])
            axs[1].set_title("Separated Energy Band Count Rates", fontsize = 12, y = 1.08)
            x_pos = 0
            for _, row in values.iterrows():
                band = row["Band"]
                color = row["Color"]
                text = row["Text"]

                # Plot the rate for each energy band
                if band == "Broadband":
                    plot.rate_plotter(plt, band, file = rates_csv, c = color, dashed = "Yes")
                    axs[1].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2, linestyle = "dotted"))
                else:
                    plot.rate_plotter(plt, band, file = rates_csv, c = color)
                    axs[1].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2))
                axs[1].text(x_pos + 0.04, 1.04, text, fontsize = 9, transform = axs[1].transAxes, va = "center", ha = "left")

                text_width = len(text) * 0.008
                x_pos += 0.0325 + text_width
            plt.grid(True, linestyle = "--")

            plt.sca(axs[2])
            axs[2].set_title("Hardness Ratios", fontsize = 12, y = 1.08)
            plot.hr_plotter(plt, file = rates_csv)
            x_pos = 0
            for _, row in plot.hr_values.iterrows():  
                formula = row["Formula"]    
                color = row["Color"]

                axs[2].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2))
                
                axs[2].text(x_pos + 0.04, 1.04, formula, fontsize = 9, transform = axs[2].transAxes, va = "center", ha = "left")       
                text_width = len(text) * 0.011
                x_pos += 0.01805 + text_width
            plt.grid(True, linestyle = "--")

            fig.suptitle(f"{instr.upper()} Corrected Light Curve for Obs. ID {obs_id} with bin size of {binsize}s", fontsize = 15, y = 0.995) 
            fig.tight_layout(pad = 1)
            fig.savefig(f"{source_name}_{obs_id}_{instr.upper()}.png")

            plt.close(fig)
    else:
        for _, row in values.iterrows():
            band = row["Band"]
            energy_min = row["Energy Min"]
            energy_max = row["Energy Max"]

            epic_lightcurve(obs_dir, energy_min, energy_max, binsize, band, instr, bkg_sub = "no")

        counts_dict = {}
        rate_dict = {}
        error_list = []
        time_list = []
        bin_list = []

        with fits.open(os.path.join(obs_dir, f"EPIC_{instr.upper()}_Source_Rate_Broadband_Lightcurve.lc")) as hdul:
            data = hdul[1].data
            for entry in data["ERROR"]:
                error_list.append(entry)

            for entry in data["TIME"]:
                time_list.append(entry)
                bin_list.append((entry - time_list[0]) / binsize + 1)

        for _, row in values.iterrows():
            band = row["Band"]
            counts_file = os.path.join(obs_dir, f"EPIC_{instr.upper()}_Source_Counts_{band}_Lightcurve.lc")
            rate_file = os.path.join(obs_dir, f"EPIC_{instr.upper()}_Source_Rate_{band}_Lightcurve.lc")

            with fits.open(counts_file) as hdul:
                data = hdul[1].data
                counts = []
                for value in data["COUNTS"]:
                    counts.append(value)
                counts_dict[band] = counts

            with fits.open(rate_file) as hdul:
                data = hdul[1].data
                rate = []
                for value in data["RATE"]:
                    rate.append(value)
                rate_dict[band] = rate

            os.remove(counts_file)
            os.remove(rate_file)
            
        time = pd.Series(time_list)
        errors = pd.Series(error_list)
        rates = pd.DataFrame(rate_dict)
        counts = pd.DataFrame(counts_dict)
        rates.insert(0, "Time", time)
        rates.insert(1, "Bin", pd.Series(map(int, bin_list)))
        rates["Error"] = errors
        counts.insert(0, "Time", time)
        counts.insert(1, "Bin", pd.Series(map(int, bin_list)))

        rates.to_csv(os.path.join(obs_dir, f"{source_name}_{obs_id}_{instr.upper()}_Count_Rates.csv"))
        counts.to_csv(os.path.join(obs_dir, f"{source_name}_{obs_id}_{instr.upper()}_Counts.csv"))

        rates_csv = os.path.join(obs_dir, f"{source_name}_{obs_id}_{instr.upper()}_Count_Rates.csv")
        counts_csv = os.path.join(obs_dir, f"{source_name}_{obs_id}_{instr.upper()}_Counts.csv")

        if rates["Broadband"].sum() == 0:
            os.remove(rates_csv)
            os.remove(counts_csv)
        else:
            fig, axs = plt.subplots(6, 1, figsize = (10, 30), dpi = 300)

            plt.sca(axs[0])
            axs[0].set_title("Broadband Count Rate", fontsize = 12)
            axs[0].text(0.9935, 1.05, f"Start: {utils.time_conv(rates["Time"].iloc[0])}", transform = axs[0].transAxes, fontsize = 8, verticalalignment = "top", horizontalalignment = "right", bbox = dict(facecolor = "white", edgecolor = "black", boxstyle = "square,pad=0.5", linewidth = 0.8))
            plot.rate_plotter(plt, file = rates_csv, band = "Broadband", error = "Yes")
            plt.grid(True, linestyle = "--")

            plt.sca(axs[1])
            axs[1].set_title("Separated Energy Band Count Rates", fontsize = 12, y = 1.08)
            x_pos = 0
            for _, row in values.iterrows():
                band = row["Band"]
                color = row["Color"]
                text = row["Text"]

                # Plot the rate for each energy band
                if band == "Broadband":
                    plot.rate_plotter(plt, band = band, file = rates_csv, c = color, dashed = "Yes")
                    axs[1].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2, linestyle = "dotted"))
                else:
                    plot.rate_plotter(plt, band = band, file = rates_csv, c = color)
                    axs[1].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2))
                    
                axs[1].text(x_pos + 0.04, 1.04, text, fontsize = 9, transform = axs[1].transAxes, va = "center", ha = "left")
                text_width = len(text) * 0.008
                x_pos += 0.0325 + text_width
            plt.grid(True, linestyle = "--")

            plt.sca(axs[2])
            axs[2].set_title("Hardness Ratios", fontsize = 12, y = 1.08)
            plot.hr_plotter(plt, file = rates_csv)
            x_pos = 0
            for _, row in plot.hr_values.iterrows():  
                formula = row["Formula"]    
                color = row["Color"]

                axs[2].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2))
                
                axs[2].text(x_pos + 0.04, 1.04, formula, fontsize = 9, transform = axs[2].transAxes, va = "center", ha = "left")       
                text_width = len(text) * 0.011
                x_pos += 0.01805 + text_width
            plt.grid(True, linestyle = "--")

            plt.sca(axs[3])
            axs[3].set_title("Broadband Counts", fontsize = 12)
            plot.counts_plotter(plt, file = counts_csv, band = "Broadband", c = "purple")
            plt.grid(True, linestyle = "--")

            plt.sca(axs[4])
            axs[4].set_title("Cumulative Counts", fontsize = 12)
            plot.cumulative_counts_plotter(plt, file = counts_csv)
            plt.grid(True, linestyle = "--")

            plt.sca(axs[5])
            axs[5].set_title("Separated Energy Band Counts", fontsize = 12, y = 1.08)
            x_pos = 0
            for _, row in values.iterrows():
                band = row["Band"]
                color = row["Color"]
                text = row["Text"]

                if band == "Broadband":
                    plot.counts_plotter(plt, band = band, file = counts_csv, c = color, dashed = "Yes")
                    axs[5].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2, linestyle = "dotted"))
                else:
                    plot.counts_plotter(plt, band = band, file = counts_csv, c = color)
                    axs[5].annotate(
                        "", xy = (x_pos + 0.03, 1.04), xytext = (x_pos, 1.04), 
                        xycoords = "axes fraction", textcoords = "axes fraction",
                        arrowprops = dict(arrowstyle = "-", color = color, lw = 2))
                    
                axs[5].text(x_pos + 0.04, 1.04, text, fontsize = 9, transform = axs[5].transAxes, va = "center", ha = "left")       
                text_width = len(text) * 0.008
                x_pos += 0.0325 + text_width
            plt.grid(True, linestyle = "--")

            fig.suptitle(f"{instr.upper()} Light Curve for Obs. ID {obs_id} with bin size of {binsize}s", fontsize = 15, y = 0.995) 
            fig.tight_layout(pad = 1)
            fig.savefig(f"{source_name}_{obs_id}_{instr.upper()}.png")

            plt.close(fig)