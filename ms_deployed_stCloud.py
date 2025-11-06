#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 15:38:25 2025

@author: danfeldheim
"""

# Imports
import pandas as pd
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib.ticker import MaxNLocator
# from scipy.signal import savgol_filter
# from scipy.stats import linregress
# from scipy.stats import t
# import scipy.stats as stats
# from statsmodels.tsa.stattools import acf
# from scipy.stats import norm
# import statsmodels.api as sm
# from statsmodels.graphics.tsaplots import plot_acf
from scipy.signal import find_peaks
from scipy.integrate import simpson
import os
# from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
# from st_aggrid.shared import JsCode
import gc
# import psutil
# import heapq
from PIL import Image
from pyteomics import mzml
import io
import time
from io import StringIO, BytesIO


# Clean up any leftover figures or memory at the start of each run
plt.close("all")
gc.collect()

class Flow_Control():
    """This class makes all of the calls to other classes and methods."""
    
    def __init__(self):
        
        # Get the path relative to the current file (inside Docker container)
        BASE_DIR = os.path.dirname(__file__)
         
        
    def all_calls(self):
        """This is the main logic workflow. All calls to other functions are here."""
        
        #---------------------------------------------------------------------------------------------
       
        state_variables_dict = {
                                "snow":False,
                                }
            
        Utilities.add_to_state(state_variables_dict)
     
        # Set up the header, login form, and check user credentials
        # Create Setup instance
        setup = Setup()
        
        # nav_bar = setup.navigation()
        
        # Render the header
        header = setup.header()
        
        if not st.session_state["snow"]:
            st.snow()
            st.session_state["snow"] = True
            
        import_data = Load_Data()
        files = import_data.upload()
        
        if files:
            select_pfas = import_data.select_pfas()
            
            if select_pfas:
                data_dict = import_data.import_to_df(files, select_pfas)
                
                analyze = Analyze_Spectra()
                
                # Plot spectra and return dict of adjusted spectra
                adjusted_spectra_dict, plot_spectra_dict = analyze.plot_mass_spectra(data_dict)
                
                # Find peaks and integrate
                # Returns a dict of filename and m/z, height, and integration
                # for all files uploaded
                peak_search = analyze.peak_finder(adjusted_spectra_dict)
                
                # Combine all peak height, integration data into a single df
                all_results_df = analyze.combine_results(peak_search)
                
                # Create a table for all results
                analyze.download_results(all_results_df, plot_spectra_dict)
                
                
            
class Setup():
    """Class that lays out the app header and sidebar."""
    
    def __init__(self):
        
        pass
    
    def header(self):
        
        # Draw line across the page
        st.divider()
        
        # Add a logo and title
        col1, col2 = st.columns([0.75,5])

        with col1:
        
            st.image(st.session_state['logo'])
            
        with col2:

            # st.write('')
            st.write('')
            st.write('')
            # st.write('')
            
            st.markdown(f"<p style='color: Blue; \
                          font-size: 32px; \
                          margin: 0;'>Biota Mass Spectrometry Analysis Tool</p>",
                          unsafe_allow_html=True
                        )
        st.divider()
        
class Load_Data():
    """Creates a file uploader and imports data as dataframe."""
    
    def __init__(self):
        
        # Dictionary of PFAS with retention time and m/z windows
       self.pfas_dict = {
                         "PFOA": {"rt_min": 0.22, "rt_max": 0.50, "mz_min": 412.60, "mz_max": 413.10},
                         "PFOS": {"rt_min": 0.65, "rt_max": 0.90, "mz_min": 498.90, "mz_max": 499.70},
                         "PFNA": {"rt_min": 0.55, "rt_max": 0.80, "mz_min": 462.30, "mz_max": 463.00},
                       }
    
    def upload(self):
        
        # File uploader that allows multiple files
        uploaded_files = st.file_uploader(
                                          "Choose files", 
                                          type=["mzML"],  
                                          accept_multiple_files=True, 
                                          )
        
        return uploaded_files
    
    def select_pfas(self):
        
        """Allows user to select a PFAS target from dropdown."""
        
        st.write('')
        st.write('')
        
        col1, col2 = st.columns([1,2])
        
        with col1:
        
            st.markdown(f"<p style='color: DarkRed; \
                  font-size: 18px; \
                  margin: 0;'>Select a PFAS for Analysis</p>",
                  unsafe_allow_html=True)
            
            pfas_name = st.selectbox("", ['Choose'] + list(self.pfas_dict.keys()), label_visibility="collapsed")
            
        if pfas_name == 'Choose':
            return None
        

        return pfas_name, self.pfas_dict[pfas_name]
    
    def import_to_df(self, uploaded_files, select_pfas):
        
        if not uploaded_files:
            return None
    
        # Unpack PFAS parameters
        pfas_name, params = select_pfas
        rt_min, rt_max = params["rt_min"], params["rt_max"]
        mz_min, mz_max = params["mz_min"], params["mz_max"]
    
        data_dict = {}
    
        for file in uploaded_files:
            
            file_name = file.name.replace(".mzML", "")
            content = io.BytesIO(file.read())
    
            with mzml.MzML(content) as reader:
                
                spectra_data = []
    
                for spectrum in reader:
                    
                    rt = spectrum.get('scanList', {}).get('scan', [{}])[0].get('scan start time', None)
                    
                    if rt is None or not (rt_min <= rt <= rt_max):
                        continue
    
                    mz_array = spectrum.get('m/z array', [])
                    intensity_array = spectrum.get('intensity array', [])
    
                    for m, intensity in zip(mz_array, intensity_array):
                        if mz_min <= m <= mz_max:
                            spectra_data.append({
                                                'retention_time': rt,
                                                'mz': m,
                                                'intensity': intensity
                                                })
    
                df = pd.DataFrame(spectra_data)
                
                data_dict[file_name] = df
    
    
        return data_dict
    
class Analyze_Spectra():
    
    def __init__(self):
        
        pass
    
    def plot_mass_spectra(self, data_dict):
        """
        Loops through a dictionary of DataFrames and plots mass spectra for each file with:
          - A progress bar for all plots
          - A two-way slider for adjusting m/z axis window
          - Proper clearing of the progress bar upon completion
          - Returns a dict of spectra after adjusting the x-axis window
        """
        if not data_dict:
            st.warning("No data to plot.")
            return
    
        total_files = len(data_dict)
    
        # Use an st.empty() placeholder instead of container for dynamic clearing
        progress_placeholder = st.empty()
        progress_bar = progress_placeholder.progress(0)
        
        # Create a dict for filtered spectra
        adjusted_spectra_dict = {}
        
        # Create a dict for plots so they can be exported as pdfs
        plot_dict = {}
    
        # Loop through dict
        for idx, (file_name, df) in enumerate(data_dict.items(), start=1):
            
            # Print file name
            st.markdown(f"<h3 style='color:dodgerblue; font-size:18pt'>{file_name}</h3>", unsafe_allow_html=True)
    
            if df.empty:
                st.warning("No data points in this file.")
                continue
    
            # Aggregate intensity by m/z to do a line plot
            spectrum = df.groupby('mz')['intensity'].sum().reset_index()
    
            # 2-way slider for adjusting m/z range
            # Get mz min/max values from dict
            mz_min, mz_max = float(spectrum['mz'].min()), float(spectrum['mz'].max())
            mz_window = st.slider(
                                  f"Adjust m/z axis for {file_name}",
                                  min_value=mz_min,
                                  max_value=mz_max,
                                  value=(mz_min, mz_max),
                                  step=(mz_max - mz_min) / 100,
                                  key=f"mz_slider_{file_name}"
                                 )
    
            # Filter spectrum for selected x-axis window
            filtered_spectrum = spectrum[
                                        (spectrum['mz'] >= mz_window[0]) & (spectrum['mz'] <= mz_window[1])
                                        ]
    
            # Plot
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(filtered_spectrum['mz'], 
                    filtered_spectrum['intensity'], 
                    color='dodgerblue', linewidth=1)
            ax.set_xlabel("m/z")
            ax.set_ylabel("Intensity")
            ax.set_title(f"Mass Spectrum: {file_name}")
            ax.grid(True, linestyle='--', alpha=0.5)
            st.pyplot(fig)
            
            # Store the plot
            plot_dict[file_name] = fig
            
            # Create slider bars to adjust peak find parameters
            max_intensity = filtered_spectrum['intensity'].max()
           
            st.markdown("**Adjust Peak Detection Parameters:**")
            
            col1, col2 = st.columns([1,1])
            
            with col1:
                
                height_slider = st.slider(
                                         f"Minimum peak height.",
                                         min_value=0.0,
                                         max_value=float(max_intensity),
                                         value=float(max_intensity * 0.05),
                                         step=float(max_intensity / 100),
                                         key=f"height_slider_{file_name}"
                                         )
                
                prominence_slider = st.slider(
                                             f"Minimum peak prominence",
                                             min_value=0.0,
                                             max_value=float(max_intensity * 0.5),
                                             value=float(max_intensity * 0.02),
                                             step=float(max_intensity / 200),
                                             key=f"prominence_slider_{file_name}"
                                             )
            
            # Store these thresholds so peak_finder can use them later
            adjusted_spectra_dict[file_name] = {
                                                "spectrum": filtered_spectrum,
                                                "height": height_slider,
                                                "prominence": prominence_slider
                                               }
    
            # Update progress bar
            progress_bar.progress(idx / total_files)
    
        # ✅ Clear progress bar entirely when done
        # Stall in case progress bar zoom by too quickly to see
        time.sleep(1)
        progress_placeholder.empty()
        
        plt.close(fig)
        del fig, ax
        gc.collect()
        
        return adjusted_spectra_dict, plot_dict
    
    def peak_finder(self, adjusted_spectra, distance=5):
        """
        Finds peaks in mass spectra and returns peak information for each file.
        Uses height and prominence sliders stored in adjusted_spectra if available.
        """
    
        all_results_dict = {}
    
        for filename, data in adjusted_spectra.items():
            
            # Unpack spectrum + detection settings
            if isinstance(data, dict):
                spectrum = data["spectrum"]
                height = data.get("height", None)
                prominence = data.get("prominence", 0.01)
                
            else:
                spectrum = data
                height = None
                prominence = 0.01
    
            # Extract x (m/z) and y (intensity)
            if isinstance(spectrum, pd.DataFrame):
                x = spectrum.iloc[:, 0].values
                y = spectrum.iloc[:, 1].values
                
            else:
                x, y = np.array(spectrum)[:, 0], np.array(spectrum)[:, 1]
    
            # Auto-estimate height if not given by slider
            if height is None:
                
                height = np.max(y) * 0.05  
    
            # Find peaks
            peaks, props = find_peaks(
                                      y,
                                      height=height,
                                      prominence=prominence,
                                      distance=distance
                                     )
    
            peak_mz = x[peaks]
            peak_intensity = y[peaks]
    
            # Integrate around each peak with a small local window
            window_size = int(len(x) * 0.001)  
            window_size = max(3, window_size)
    
            peak_areas = []
            for p in peaks:
                left = max(0, p - window_size)
                right = min(len(x) - 1, p + window_size)
                area = simpson(y[left:right], x[left:right])
                peak_areas.append(area)
    
            peak_df = pd.DataFrame({
                                    "m/z": peak_mz,
                                    "peak_intensity": peak_intensity,
                                    "peak_area": peak_areas
                                   })
            
            # Sort from largest to smallest peak height
            peak_df = peak_df.sort_values(by="peak_intensity", 
                                          ascending=False).reset_index(drop=True)
            
            # Create column_config for all columns as "small"
            column_config = {
                            col: st.column_config.Column(col, width=200)
                            for col in peak_df.columns
                            }
            
            # st.markdown("### 📈 Peaks Detected ")
            st.markdown(f"<p style='color: Blue; \
                          font-size: 24px; \
                          margin: 0;'>📈 Peaks Detected</p>",
                          unsafe_allow_html=True
                        )
                
            # Number the entries
            peak_df.insert(0, "Peak #", np.arange(1, len(peak_df) + 1))
            
            current_file_results_table = st.data_editor(
                                                        peak_df,
                                                        hide_index=True,
                                                        use_container_width=False,
                                                        num_rows="static",
                                                        column_config=column_config,
                                                        key=filename
                                                        )
            
            st.divider()
            
    
            all_results_dict[filename] = peak_df.sort_values(by="m/z").reset_index(drop=True)
 
        
        return all_results_dict
    
    def combine_results(self, all_results):
        
        # Convert all_results to a new df with filenames as a new column
        combined_list = []

        for filename, df in all_results.items():
            
            # Add a new column for the filename
            df_with_filename = df.copy()
            df_with_filename.insert(0, 'filename', filename) 
            combined_list.append(df_with_filename)
    
        # Combine all DataFrames into one
        all_results_df = pd.concat(combined_list, ignore_index=True)
        
        return all_results_df
    
    def download_results(self, all_results_df, all_spectra):
        
        st.markdown(f"<p style='color: Blue; \
                      font-size: 24px; \
                      margin: 0;'>📈 Results for all Files Uploaded</p>",
                      unsafe_allow_html=True
                    )
        
        # Create column_config for all columns as "small"
        column_config = {
                        col: st.column_config.Column(col, width=200)
                        for col in all_results_df.columns
                        }
        
        all_results_table = st.data_editor(
                                           all_results_df,
                                           hide_index=True,
                                           use_container_width=False,
                                           num_rows="static",
                                           column_config=column_config,
                                           )
        
        # Convert DataFrame to CSV (in memory)
        csv_buffer = StringIO()
        all_results_df.to_csv(csv_buffer, index=False)
        csv_data = csv_buffer.getvalue()
    
        # Streamlit download button
        st.download_button(
                           label="📥 Download All Results",
                           data=csv_data,
                           file_name="all_results_data.csv",
                           mime="text/csv"
                          )
        
        pdf_buffer = BytesIO()
        with PdfPages(pdf_buffer) as pdf:
            for filename, fig in all_spectra.items():
                # Optional: add a title to each page
                # fig.suptitle(filename, fontsize=12)
                pdf.savefig(fig)
                plt.close(fig) 

        pdf_buffer.seek(0)  

        st.download_button(
                           label="📄 Download All Spectra",
                           data=pdf_buffer,
                           file_name="mass_spectra_plots.pdf",
                           mime="application/pdf"
                           )
        
        st.divider()
        
        
            


        
  
class Utilities():
    """Contains static methods that can be accessed from all other classes."""
    
    def __init__(self):
        
        pass
    
    @staticmethod
    def add_to_state(state_variables, overwrite=False):
        """
        Instantiates session state variables from a dictionary.
        If overwrite=True, existing keys will be updated.
        """
        for key, value in state_variables.items():
            if overwrite or key not in st.session_state:
                st.session_state[key] = value
                
    @staticmethod
    def progress_bar(message: str, current: int, total: int, bar=None):
        """
        Updates or creates a progress bar.
        Args:
            message (str): Label text above progress bar.
            current (int): Current step.
            total (int): Total number of steps.
            bar (st.progress): Existing progress bar (optional).
        Returns:
            Updated progress bar object.
        """
        if total == 0:
            return None
        
        progress = min(current / total, 1.0)
        if bar is None:
            st.caption(message)
            bar = st.progress(progress)
        else:
            bar.progress(progress)
        return bar




# Run 
if __name__ == '__main__':
    
    # Set up session_state variables
    
    # Delete directory when deploying to cloud
    directory = '/Users/danfeldheim/Documents/ms_app/'
        
    # Use this for cloud
    st.session_state['logo'] = 'logo-black.png'
    
    # Use this for local machine
    # st.session_state['logo'] = directory + 'logo-black.png'
        
    # Load image for favicon
    logo_img = Image.open(st.session_state['logo'])
        
    # Page config
    st.set_page_config(layout = "wide", 
                       page_title = 'Biota', 
                       page_icon = logo_img,
                       initial_sidebar_state="auto", 
                       menu_items = None)
    
    
    # Call Flow_Control class that makes all calls to other classes and methods
    obj1 = Flow_Control()
    all_calls = obj1.all_calls()
