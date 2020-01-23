import ConvertTDTToNWB
from datetime import datetime
from os import getcwd, path, listdir
import NWBHelperFunctions

# TODO: Write wrapper to use on cluster

# File to import
session_name = "Mouse1" # TODO: Fix the session name
tdt_block_dir = "C:\\Users\\tbuschma\\OneDrive\\Documents\\Python Scripts\\NWB conversion\\Example Data\\" \
                "Mouse_Ephys_Imaging\\Block-1\\" # TODO: Fix the block location for looping




#######################
## Initialize NWB file

# Information about the experiment
experiment = {'name': 'Ephys+MesoscaleDynamics', # TODO: Project Name
              'animal': 'Mouse1', # TODO: Need to fill in animal's name
              'experimenter': 'Camden MacDowell and Daniel Che',
              'lab': 'Buschman (www.timbuschman.com)',
              'institution': 'Princeton Neuroscience Institute'}
subject = {'subject_id': 'Mouse1, Tag #XXXX', 'age': 'XX', 'description': 'A very good mouse.',
           'genotype': 'N/A', 'sex': 'M', 'species': 'Thy1-GCaMP6', 'weight': 'XX g',
           'date_of_birth': datetime.strptime('12/12/2019', '%m/%d/%Y')} # TODO: Fill this out

notes = 'Animal was positioned on a running wheel...' #TODO: Add any notes you want

# Electrode groups should bring together of sites on the same electrode (i.e. for laminar arrays)
elec_type = 'PR'
electrode_group = {'name': 'default electrode group', # TODO: Add name of electrode group -- I would only create groups for those electrodes tied together into a single group (e.g. laminar probes, etc)
                   'description': 'Generic electrode group for ungrouped electrodes. Electrode type: %s' % (elec_type), # TODO: Add information about the electrode group, such as its type.
                   'location': 'Area'} # TODO: Import area from unit dataframe

# Information about electrodes -- can probably be extracted from ElectrodeTracker data
electrode_info = []  # TODO: Add code to import electrode information
electrode_info.append({'id': 1, #TODO: This is the channel number -- assuming that RAW1 is 1-16, RAW2 is 17-32, etc
                       'x': 0.0, 'y': 0.0, 'z': 0.0, # TODO: Location of electrode (A/P, M/L, relative to bregma
                       'impedance': float(-1), # TODO: We probably don't have this, so just leave as -1
                       'location': 'PAR', # TODO: region assignment
                       'group': None, # TODO: This should reference the above group -- use the name of the electrode group for readability
                       'filtering': ''}),  # TODO: Add any filtering -- not included in TDT files
electrode_info.append({'id': 2, 'x': 0.0, 'y': 0.0, 'z': 0.0, 'impedance': float(-1), 'location': 'PAR', 'group': None, 'filtering': ''})

# Initialize NWB file with information
nwb_file_name = ConvertTDTToNWB.InitializeNWBFromTDT(tdt_block_dir,
                                                     '',
                                                     experiment=experiment,
                                                     electrode_info=electrode_info,
                                                     electrode_group=electrode_group,
                                                     subject=subject,
                                                     notes=notes,
                                                     overwrite=True,
                                                     verbose=True)


## Add raw stimulus images to file -- if you have any stimuli you show to the animal that you want to save to file...
#stimulus_image_info = [{'file': path.dirname(blk_file_name) + '\\101.jpg', 'name': 'Stimulus 1 (101)', 'resolution': None, 'description': 'First stimulus in Group 1'}]
#NWBHelperFunctions.AddStimulusImageToNWB(nwb_file_name, stimulus_image_info, verbose=True)



#######################
## Add raw electrophysiology data
elec_ids=[elec['id'] for elec in electrode_info]
ConvertTDTToNWB.AddTDTRawDataToNWB(tdt_block_dir, nwb_file_name, elec_ids=elec_ids, verbose=True)


#######################
## Add camera signal data
signal_info = [{'label': 'CAME', 'name': 'Camera', 'comments': 'Analog signal recording camera captures.', 'elec_id': ['CAME']}]
ConvertTDTToNWB.AddTDTAnalogDataToNWB(tdt_block_dir, nwb_file_name, signal_info=signal_info, verbose=True)

#######################
## Process camera signal to extract onset times
nwb_file, nwb_io = NWBHelperFunctions.OpenNWBFile(nwb_file_name, mode='r')
# Grab analog signal for camera image acquisition
analog_signal = nwb_file.processing['behavior']['Analog signals (Camera)']['Camera']
smooth = {'type': 'boxcar', 'width': 5}
onset_times, offset_times = NWBHelperFunctions.ExtractTimestampsFromAnalogSignal(analog_signal=analog_signal, threshold=2.0, smooth=smooth, verbose=True)
nwb_io.close()


#######################
## Add mesoscale image data to NWB

# Get list of image files #TODO: you'll probably want to make sure that this is correct, rather than doing it programmatically
img_dir = "C:\\Users\\tbuschma\\OneDrive\\Documents\\Python Scripts\\NWB conversion\\Example Data\\" \
                "Mouse_Ephys_Imaging\\" # TODO: Fix the image location for looping
img_files = [(img_dir + f) for f in listdir(img_dir) if path.splitext(f)[1] == '.tif']

nwb_file, nwb_io = NWBHelperFunctions.OpenNWBFile(nwb_file_name, mode='a')
NWBHelperFunctions.AddMesoscaleImagingDataToNWB(nwb_file,   # This shows a different way of interacting with NWB -- rather than having function open and close it, you can often pass the file itself
                                                tiff_files=img_files,
                                                img_info=None,
                                                timestamps=onset_times[0:4123], #TODO: Make sure the timestamps are the same length as the image stack
                                                load_data=True,
                                                verbose=True)
# Write to disk
nwb_io.write(nwb_file)
nwb_io.close()